//! Codon and amino acid encoding.
//!
//! Bases are `A=0, C=1, G=2, T=3`, and a codon index is `b0 * 16 + b1 * 4 + b2`, so
//! codons occupy `0..64`. That layout makes a rolling k-mer a shift and a mask, which
//! the repeat and hairpin indices rely on.

pub type Base = u8;
pub type CodonIdx = u8;
pub type AaIdx = u8;

pub const N_CODONS: usize = 64;
/// 20 amino acids plus stop.
pub const N_AA: usize = 21;
pub const STOP: AaIdx = 20;

pub const BASE_LETTER: &[u8; 4] = b"ACGT";
pub const AA_LETTER: &[u8; N_AA] = b"ACDEFGHIKLMNPQRSTVWY*";

/// Standard genetic code, indexed by codon index. Generated from Biopython's
/// `CodonTable.unambiguous_dna_by_name["Standard"]`, stops written as `*`.
const CODON_AA_LETTER: &[u8; N_CODONS] =
    b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

/// Translation table from codon index to amino acid index.
pub fn codon_aa(codon: CodonIdx) -> AaIdx {
    aa_index(CODON_AA_LETTER[codon as usize]).expect("genetic code holds only valid letters")
}

pub fn base_index(letter: u8) -> Option<Base> {
    match letter.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

pub fn aa_index(letter: u8) -> Option<AaIdx> {
    AA_LETTER
        .iter()
        .position(|&l| l == letter.to_ascii_uppercase())
        .map(|p| p as AaIdx)
}

pub fn aa_letter(aa: AaIdx) -> u8 {
    AA_LETTER[aa as usize]
}

pub fn codon_bases(codon: CodonIdx) -> [Base; 3] {
    [codon >> 4, (codon >> 2) & 3, codon & 3]
}

pub fn codon_from_bases(b: [Base; 3]) -> CodonIdx {
    (b[0] << 4) | (b[1] << 2) | b[2]
}

pub fn codon_letters(codon: CodonIdx) -> [u8; 3] {
    let b = codon_bases(codon);
    [
        BASE_LETTER[b[0] as usize],
        BASE_LETTER[b[1] as usize],
        BASE_LETTER[b[2] as usize],
    ]
}

pub fn codon_string(codon: CodonIdx) -> String {
    String::from_utf8(codon_letters(codon).to_vec()).expect("ASCII")
}

pub fn parse_codon(text: &str) -> Option<CodonIdx> {
    let bytes = text.as_bytes();
    if bytes.len() != 3 {
        return None;
    }
    let b = [
        base_index(bytes[0])?,
        base_index(bytes[1])?,
        base_index(bytes[2])?,
    ];
    Some(codon_from_bases(b))
}

/// Number of G or C bases in a codon.
pub fn codon_gc(codon: CodonIdx) -> u32 {
    codon_bases(codon)
        .iter()
        .filter(|&&b| b == 1 || b == 2)
        .count() as u32
}

pub fn is_gc(base: Base) -> bool {
    base == 1 || base == 2
}

pub fn is_purine(base: Base) -> bool {
    base == 0 || base == 2
}

pub fn complement(base: Base) -> Base {
    3 - base
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codon_index_round_trips() {
        for c in 0..N_CODONS as CodonIdx {
            assert_eq!(codon_from_bases(codon_bases(c)), c);
            assert_eq!(parse_codon(&codon_string(c)), Some(c));
        }
    }

    #[test]
    fn genetic_code_matches_known_codons() {
        let cases = [
            ("ATG", b'M'),
            ("TGG", b'W'),
            ("TTT", b'F'),
            ("AAA", b'K'),
            ("CGT", b'R'),
            ("AGA", b'R'),
            ("GCA", b'A'),
            ("TAA", b'*'),
            ("TAG", b'*'),
            ("TGA", b'*'),
        ];
        for (text, letter) in cases {
            let codon = parse_codon(text).unwrap();
            assert_eq!(aa_letter(codon_aa(codon)), letter, "codon {text}");
        }
    }

    #[test]
    fn genetic_code_has_the_expected_group_sizes() {
        let mut sizes = [0usize; N_AA];
        for c in 0..N_CODONS as CodonIdx {
            sizes[codon_aa(c) as usize] += 1;
        }
        // Leu, Ser and Arg have six codons; Met and Trp have one; stops have three.
        assert_eq!(sizes[aa_index(b'L').unwrap() as usize], 6);
        assert_eq!(sizes[aa_index(b'S').unwrap() as usize], 6);
        assert_eq!(sizes[aa_index(b'R').unwrap() as usize], 6);
        assert_eq!(sizes[aa_index(b'M').unwrap() as usize], 1);
        assert_eq!(sizes[aa_index(b'W').unwrap() as usize], 1);
        assert_eq!(sizes[STOP as usize], 3);
        assert_eq!(sizes.iter().sum::<usize>(), N_CODONS);
    }

    #[test]
    fn gc_and_purine_helpers() {
        assert_eq!(codon_gc(parse_codon("GCC").unwrap()), 3);
        assert_eq!(codon_gc(parse_codon("ATT").unwrap()), 0);
        assert_eq!(codon_gc(parse_codon("ATG").unwrap()), 1);
        assert!(is_purine(base_index(b'A').unwrap()));
        assert!(is_purine(base_index(b'G').unwrap()));
        assert!(!is_purine(base_index(b'C').unwrap()));
        assert_eq!(complement(base_index(b'A').unwrap()), base_index(b'T').unwrap());
        assert_eq!(complement(base_index(b'G').unwrap()), base_index(b'C').unwrap());
    }
}
