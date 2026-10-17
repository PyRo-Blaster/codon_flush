//! Minimal FASTA reading and writing.

use anyhow::{bail, Context, Result};

use crate::code::{aa_index, AaIdx};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub id: String,
    pub description: String,
    pub residues: Vec<AaIdx>,
}

impl Record {
    pub fn header(&self) -> String {
        if self.description.is_empty() {
            self.id.clone()
        } else {
            format!("{} {}", self.id, self.description)
        }
    }
}

/// Parse protein records. Lowercase is accepted, whitespace and digits are ignored, and
/// a trailing `*` is kept as a stop so the gene can carry its own terminator.
pub fn parse(text: &str) -> Result<Vec<Record>> {
    let mut records: Vec<Record> = Vec::new();
    let mut current: Option<Record> = None;
    for (number, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with(';') {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(record) = current.take() {
                records.push(record);
            }
            let header = header.trim();
            let (id, description) = match header.split_once(char::is_whitespace) {
                Some((id, rest)) => (id.to_string(), rest.trim().to_string()),
                None => (header.to_string(), String::new()),
            };
            current = Some(Record { id, description, residues: Vec::new() });
            continue;
        }
        let Some(record) = current.as_mut() else {
            bail!("line {}: sequence data before any > header", number + 1);
        };
        for letter in line.bytes() {
            if letter.is_ascii_whitespace() || letter.is_ascii_digit() || letter == b'-' {
                continue;
            }
            let residue = aa_index(letter).with_context(|| {
                format!(
                    "line {}: {:?} is not an amino acid; use the 20 standard letters or *",
                    number + 1,
                    letter as char
                )
            })?;
            record.residues.push(residue);
        }
    }
    if let Some(record) = current {
        records.push(record);
    }
    if records.is_empty() {
        bail!("no FASTA records found");
    }
    if let Some(empty) = records.iter().find(|r| r.residues.is_empty()) {
        bail!("record {:?} has no residues", empty.id);
    }
    Ok(records)
}

/// One record, wrapped at 60 bases like the Python wrote them.
pub fn write_record(header: &str, dna: &str) -> String {
    let mut out = format!(">{header}\n");
    for chunk in dna.as_bytes().chunks(60) {
        out.push_str(std::str::from_utf8(chunk).expect("ASCII"));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::aa_letter;

    fn letters(record: &Record) -> String {
        record.residues.iter().map(|&a| aa_letter(a) as char).collect()
    }

    #[test]
    fn reads_multi_line_records() {
        let text = ">one first protein\nMSKG\nEELF\n\n>two\nmskg\n";
        let records = parse(text).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "one");
        assert_eq!(records[0].description, "first protein");
        assert_eq!(letters(&records[0]), "MSKGEELF");
        assert_eq!(records[1].id, "two");
        assert_eq!(letters(&records[1]), "MSKG", "lowercase is accepted");
    }

    #[test]
    fn keeps_a_stop_codon() {
        let records = parse(">x\nMSK*\n").unwrap();
        assert_eq!(letters(&records[0]), "MSK*");
    }

    #[test]
    fn ignores_numbering_and_gaps() {
        let records = parse(">x\n1 MSK 4 G-G\n").unwrap();
        assert_eq!(letters(&records[0]), "MSKGG");
    }

    #[test]
    fn rejects_bad_input_clearly() {
        let err = parse(">x\nMSKZ\n").unwrap_err().to_string();
        assert!(err.contains('Z'), "{err}");
        assert!(parse("MSKG\n").unwrap_err().to_string().contains("before any"));
        assert!(parse("").unwrap_err().to_string().contains("no FASTA records"));
        assert!(parse(">empty\n").unwrap_err().to_string().contains("no residues"));
    }

    #[test]
    fn writes_wrapped_records() {
        let dna = "A".repeat(130);
        let out = write_record("gene one", &dna);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], ">gene one");
        assert_eq!(lines[1].len(), 60);
        assert_eq!(lines[3].len(), 10);
        assert_eq!(lines[1..].concat(), dna);
    }
}
