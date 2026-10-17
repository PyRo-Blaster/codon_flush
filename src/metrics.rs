//! Numbers reported about a finished gene.

use crate::arrange::state::SeqState;
use crate::code::{codon_aa, is_gc, CodonIdx, N_CODONS};
use crate::composition;
use crate::host::HostTable;

pub use crate::composition::{deviation, gc_fraction};

/// Codon adaptation index: the geometric mean of the relative adaptiveness of every
/// codon, skipping Met and Trp, which have no alternatives.
///
/// Biopython divides the log sum by `len - 1` rather than `len`; this uses `len`, the
/// textbook definition. The difference is well under a percent at gene length.
pub fn cai(state: &SeqState, host: &HostTable) -> f64 {
    cai_of_counts(&state.counts, host)
}

/// CAI from codon counts alone. Summing over the 64 counts rather than walking the
/// sequence keeps the result bit-identical for any arrangement of the same multiset:
/// floating-point addition is not associative, so summation order would otherwise
/// show up in the last digit.
pub fn cai_of_counts(counts: &[u32; N_CODONS], host: &HostTable) -> f64 {
    let mut sum = 0.0;
    let mut counted = 0u32;
    for codon in 0..N_CODONS as CodonIdx {
        let n = counts[codon as usize];
        if n == 0 || host.synonyms(codon_aa(codon)).len() < 2 {
            continue;
        }
        let weight = host.weight(codon);
        if weight <= 0.0 {
            return 0.0; // a codon the host never uses collapses the whole index
        }
        sum += weight.ln() * n as f64;
        counted += n;
    }
    if counted == 0 {
        return 1.0;
    }
    (sum / counted as f64).exp()
}

pub fn counts_of(state: &SeqState) -> [u32; N_CODONS] {
    state.counts
}

pub fn deviation_of(state: &SeqState, host: &HostTable) -> f64 {
    composition::deviation(&state.counts, host)
}

pub fn gc_of(state: &SeqState) -> f64 {
    let gc = state.nt.iter().filter(|&&b| is_gc(b)).count();
    if state.nt.is_empty() {
        0.0
    } else {
        gc as f64 / state.nt.len() as f64
    }
}

/// Codons used that the host threshold removed. Should always be empty.
pub fn unused_codons(state: &SeqState, host: &HostTable) -> Vec<CodonIdx> {
    let mut out: Vec<CodonIdx> = (0..N_CODONS as CodonIdx)
        .filter(|&c| state.counts[c as usize] > 0 && host.freq(c) == 0.0)
        .collect();
    out.sort_unstable();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::{aa_index, parse_codon, AaIdx};
    use crate::testdata::ECOLI_TABLE;
    use std::path::Path;

    fn host() -> HostTable {
        HostTable::load(Path::new(ECOLI_TABLE), 0.10).unwrap()
    }

    fn state_of(dna: &str) -> SeqState {
        let codons: Vec<CodonIdx> = dna
            .as_bytes()
            .chunks(3)
            .map(|c| parse_codon(std::str::from_utf8(c).unwrap()).unwrap())
            .collect();
        let aa: Vec<AaIdx> = codons.iter().copied().map(codon_aa).collect();
        SeqState::new(aa, codons)
    }

    #[test]
    fn cai_of_all_preferred_codons_is_one() {
        // CGC and GCG are the most used Arg and Ala codons in this table
        let state = state_of("CGCGCGCGCGCG");
        assert!((cai(&state, &host()) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cai_collapses_on_a_codon_the_host_dropped() {
        let state = state_of("CGCAGACGCGCG"); // AGA is below the 10% threshold
        assert_eq!(cai(&state, &host()), 0.0);
    }

    #[test]
    fn met_and_trp_do_not_count_towards_cai() {
        let only_single = state_of("ATGTGGATGTGG");
        assert_eq!(cai(&only_single, &host()), 1.0);
        let _ = aa_index(b'M');
    }

    #[test]
    fn gc_is_counted_over_bases() {
        assert!((gc_of(&state_of("GGGCCC")) - 1.0).abs() < 1e-12);
        assert!((gc_of(&state_of("AAATTT")) - 0.0).abs() < 1e-12);
        assert!((gc_of(&state_of("GGGAAA")) - 0.5).abs() < 1e-12);
    }
}
