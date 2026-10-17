//! Host codon usage: parsing, thresholding and CAI weights.
//!
//! Two JSON shapes are accepted, because both exist in the wild and in this repo:
//! nested `{"A": {"GCA": 0.21, ...}, ...}` as produced by `python_codon_tables`, and
//! flat `{"GCA": 0.21, ...}`. Either way the frequencies are renormalised within each
//! amino acid, so a table whose single-codon amino acids do not sum to 1 is repaired
//! rather than trusted.

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::code::{
    aa_index, codon_aa, codon_string, parse_codon, AaIdx, CodonIdx, N_AA, N_CODONS,
};

#[derive(Debug, Clone)]
pub struct HostTable {
    freq: [f64; N_CODONS],
    weight: [f64; N_CODONS],
    by_aa: Vec<Vec<CodonIdx>>,
    threshold: f64,
}

impl HostTable {
    pub fn load(path: &Path, threshold: f64) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading host table {}", path.display()))?;
        Self::from_json_str(&text, threshold)
            .with_context(|| format!("parsing host table {}", path.display()))
    }

    pub fn from_json_str(text: &str, threshold: f64) -> Result<Self> {
        if !(0.0..1.0).contains(&threshold) {
            bail!("host threshold must be in [0, 1), got {threshold}");
        }
        let raw = parse_any_shape(text)?;
        let mut freq = [0.0f64; N_CODONS];
        for (codon, value) in raw {
            freq[codon as usize] = value;
        }

        let mut by_aa: Vec<Vec<CodonIdx>> = vec![Vec::new(); N_AA];
        for codon in 0..N_CODONS as CodonIdx {
            by_aa[codon_aa(codon) as usize].push(codon);
        }

        renormalise(&mut freq, &by_aa)?;
        for value in freq.iter_mut() {
            if *value < threshold {
                *value = 0.0;
            }
        }
        renormalise(&mut freq, &by_aa)?;

        let mut weight = [0.0f64; N_CODONS];
        for group in &by_aa {
            let max = group
                .iter()
                .map(|&c| freq[c as usize])
                .fold(0.0f64, f64::max);
            if max <= 0.0 {
                continue;
            }
            for &c in group {
                weight[c as usize] = freq[c as usize] / max;
            }
        }

        Ok(Self {
            freq,
            weight,
            by_aa,
            threshold,
        })
    }

    pub fn freq(&self, codon: CodonIdx) -> f64 {
        self.freq[codon as usize]
    }

    /// Relative adaptiveness, i.e. the codon's frequency over the most used frequency
    /// in its synonymous group. This is the per-codon term of the CAI.
    pub fn weight(&self, codon: CodonIdx) -> f64 {
        self.weight[codon as usize]
    }

    /// Every codon for this amino acid, including ones the threshold dropped.
    pub fn synonyms(&self, aa: AaIdx) -> &[CodonIdx] {
        &self.by_aa[aa as usize]
    }

    /// Synonyms the host actually uses. Empty only if an amino acid is missing from the
    /// table entirely, which `from_json_str` rejects.
    pub fn usable_synonyms(&self, aa: AaIdx) -> impl Iterator<Item = CodonIdx> + '_ {
        self.by_aa[aa as usize]
            .iter()
            .copied()
            .filter(move |&c| self.freq[c as usize] > 0.0)
    }

    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// Codons the threshold removed, for reporting.
    pub fn dropped(&self) -> Vec<CodonIdx> {
        (0..N_CODONS as CodonIdx)
            .filter(|&c| self.freq[c as usize] == 0.0)
            .collect()
    }
}

fn renormalise(freq: &mut [f64; N_CODONS], by_aa: &[Vec<CodonIdx>]) -> Result<()> {
    for (aa, group) in by_aa.iter().enumerate() {
        let total: f64 = group.iter().map(|&c| freq[c as usize]).sum();
        if total <= 0.0 {
            let letter = crate::code::aa_letter(aa as AaIdx) as char;
            bail!("host table has no usable codon for amino acid {letter}");
        }
        for &c in group {
            freq[c as usize] /= total;
        }
    }
    Ok(())
}

/// Accepts either the nested per-amino-acid shape or a flat codon map.
fn parse_any_shape(text: &str) -> Result<Vec<(CodonIdx, f64)>> {
    let value: serde_json::Value =
        serde_json::from_str(text).context("host table is not valid JSON")?;
    let object = value
        .as_object()
        .context("host table must be a JSON object")?;

    let mut out: BTreeMap<CodonIdx, f64> = BTreeMap::new();
    for (key, entry) in object {
        match entry {
            serde_json::Value::Object(inner) => {
                // nested: the key is an amino acid letter
                if key.len() != 1 || aa_index(key.as_bytes()[0]).is_none() {
                    bail!("unknown amino acid key {key:?} in host table");
                }
                for (codon_text, freq) in inner {
                    insert_codon(&mut out, codon_text, freq)?;
                }
            }
            _ => insert_codon(&mut out, key, entry)?,
        }
    }
    if out.is_empty() {
        bail!("host table has no codons");
    }
    Ok(out.into_iter().collect())
}

fn insert_codon(
    out: &mut BTreeMap<CodonIdx, f64>,
    codon_text: &str,
    value: &serde_json::Value,
) -> Result<()> {
    let codon = parse_codon(codon_text)
        .with_context(|| format!("{codon_text:?} is not a DNA or RNA codon"))?;
    let freq = value
        .as_f64()
        .with_context(|| format!("frequency for {codon_text} is not a number"))?;
    if !freq.is_finite() || freq < 0.0 {
        bail!("frequency for {codon_text} must be zero or more, got {freq}");
    }
    if out.insert(codon, freq).is_some() {
        bail!("codon {} appears twice in the host table", codon_string(codon));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::aa_letter;

    const ECOLI: &str = "tests/data/e_coli_316407.json";

    fn ecoli(threshold: f64) -> HostTable {
        HostTable::load(Path::new(ECOLI), threshold).expect("E. coli table loads")
    }

    fn freq_of(table: &HostTable, codon: &str) -> f64 {
        table.freq(parse_codon(codon).unwrap())
    }

    fn weight_of(table: &HostTable, codon: &str) -> f64 {
        table.weight(parse_codon(codon).unwrap())
    }

    #[test]
    fn parses_nested_table() {
        let table =
            HostTable::from_json_str(include_str!("../tests/data/e_coli_316407.json"), 0.0)
                .unwrap();
        assert!((freq_of(&table, "GCA") - 0.21).abs() < 1e-12);
    }

    #[test]
    fn parses_flat_table_identically() {
        let nested = ecoli(0.0);
        let flat: String = {
            let mut parts = Vec::new();
            for c in 0..N_CODONS as CodonIdx {
                parts.push(format!("\"{}\": {}", codon_string(c), nested.freq(c)));
            }
            format!("{{{}}}", parts.join(","))
        };
        let reparsed = HostTable::from_json_str(&flat, 0.0).unwrap();
        for c in 0..N_CODONS as CodonIdx {
            assert!((reparsed.freq(c) - nested.freq(c)).abs() < 1e-12);
        }
    }

    #[test]
    fn renormalises_per_amino_acid() {
        // Trp has a single codon, so whatever the file claims, it must come out at 1.0.
        let table = HostTable::from_json_str(r#"{"TGG": 0.48, "ATG": 0.5}"#, 0.0);
        // every other amino acid is missing, so this must be an error, not silence
        assert!(table.is_err(), "a table missing amino acids must be rejected");

        let table = ecoli(0.0);
        assert!((freq_of(&table, "TGG") - 1.0).abs() < 1e-12);
        assert!((freq_of(&table, "ATG") - 1.0).abs() < 1e-12);
    }

    #[test]
    fn every_amino_acid_group_sums_to_one() {
        for threshold in [0.0, 0.1, 0.2] {
            let table = ecoli(threshold);
            for aa in 0..N_AA as AaIdx {
                let total: f64 = table.synonyms(aa).iter().map(|&c| table.freq(c)).sum();
                assert!(
                    (total - 1.0).abs() < 1e-12,
                    "amino acid {} sums to {total} at threshold {threshold}",
                    aa_letter(aa) as char
                );
            }
        }
    }

    #[test]
    fn threshold_zeroes_and_renormalises_like_python() {
        // Values below are what the Python package produces for this table at 0.10;
        // they were read out of the verified viz-demo/pydata.json reference.
        let table = ecoli(0.10);
        for dropped in ["AGA", "AGG", "ATA", "CGA", "CTA", "TAG"] {
            assert_eq!(freq_of(&table, dropped), 0.0, "{dropped} should be dropped");
        }
        let expected = [
            ("CGT", 0.431_818_181_818_181_8),
            ("CGC", 0.454_545_454_545_454_5),
            ("CGG", 0.113_636_363_636_363_6),
            ("CTG", 0.520_833_333_333_333_3),
            ("TTA", 0.135_416_666_666_666_7),
            ("GCA", 0.21),
            ("GCG", 0.36),
        ];
        for (codon, want) in expected {
            assert!(
                (freq_of(&table, codon) - want).abs() < 1e-9,
                "{codon}: got {}, want {want}",
                freq_of(&table, codon)
            );
        }
    }

    #[test]
    fn cai_weights_are_relative_to_group_max() {
        let table = ecoli(0.10);
        let expected = [
            ("CGC", 1.0),
            ("CGT", 0.95),
            ("CGG", 0.25),
            ("CGA", 0.0),
            ("CTG", 1.0),
            ("TTA", 0.26),
            ("GCG", 1.0),
            ("GCA", 0.583_333_333_333_333_3),
            ("GCT", 0.444_444_444_444_444_4),
            ("TGG", 1.0),
        ];
        for (codon, want) in expected {
            assert!(
                (weight_of(&table, codon) - want).abs() < 1e-9,
                "{codon}: got {}, want {want}",
                weight_of(&table, codon)
            );
        }
        // the max weight in every group is exactly 1.0
        for aa in 0..N_AA as AaIdx {
            let max = table
                .synonyms(aa)
                .iter()
                .map(|&c| table.weight(c))
                .fold(0.0f64, f64::max);
            assert!((max - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn usable_synonyms_skip_dropped_codons() {
        let table = ecoli(0.10);
        let arg = aa_index(b'R').unwrap();
        let usable: Vec<String> = table.usable_synonyms(arg).map(codon_string).collect();
        assert_eq!(usable, vec!["CGC", "CGG", "CGT"]);
        assert_eq!(table.synonyms(arg).len(), 6);
    }

    #[test]
    fn rejects_bad_input() {
        assert!(HostTable::from_json_str("not json", 0.0).is_err());
        assert!(HostTable::from_json_str(r#"{"XYZ": 1.0}"#, 0.0).is_err());
        assert!(HostTable::from_json_str(r#"{"GCA": -1.0}"#, 0.0).is_err());
        assert!(HostTable::load(Path::new(ECOLI), 1.5).is_err());
    }
}
