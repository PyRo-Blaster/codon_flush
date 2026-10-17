//! RNA splice donor and acceptor motifs.
//!
//! The 20 donor and 10 acceptor patterns are the ones the Python carries, transcribed
//! as base-class masks instead of regular expressions. Every donor has `GT` fixed at
//! offset 2 and every acceptor has `CAG` at offset 6, which gives a cheap prefilter:
//! almost every position is rejected by one comparison.

use std::ops::Range;

use crate::arrange::state::{SeqState, View};
use crate::code::Base;

use super::{window_delta, Constraint, Kind, Violation};

/// Consensus donor is GGGTRAGT: the first GT is fixed and at least two other bases
/// match the consensus.
const DONORS: [&str; 20] = [
    "GGGT.AGT",
    ".GGT.[AT]GT",
    "G.GT.AGT",
    "..GT[AG]AGT",
    "GGGT[AG].[AG]T",
    "GGGT[AG].G.",
    "GGGT[AG]A[AG].",
    "GGGT[AG]AG.",
    "GGGT[AG]A[AG]T",
    "GGGT[AG].GT",
    ".GGT[AG][AT][ATG].",
    ".GGT[AG].GT",
    "G.GT[AG]A[AG].",
    "G.GT[AG]A[AG].",
    "G.GT[AG]AG.",
    "G.GT[AG]A[AG]T",
    "G.GT[AG].GT",
    "..GT[AG]AG.",
    "..GT[AG].GT",
    "..GT[AG].G.",
];

/// Consensus acceptor is YYYYYNCAGG: at least three of the five upstream bases are
/// pyrimidines.
const ACCEPTORS: [&str; 10] = [
    "[TC][TC][TC]..[ATCG]CAG.",
    "[TC][TC].[TC].[ATCG]CAG.",
    "[TC].[TC][TC].[ATCG]CAG.",
    ".[TC][TC][TC].[ATCG]CAG.",
    ".[TC][TC].[TC][ATCG]CAG.",
    ".[TC].[TC][TC][ATCG]CAG.",
    "..[TC][TC][TC][ATCG]CAG.",
    "[TC].[TC].[TC][ATCG]CAG.",
    "[TC]..[TC][TC][ATCG]CAG.",
    "[TC][TC]..[TC][ATCG]CAG.",
];

/// One bit per base: A=1, C=2, G=4, T=8.
type Mask = u8;

fn mask_of(base: Base) -> Mask {
    1 << base
}

/// Parse `GGT[AG].` into one mask per position.
fn parse_pattern(pattern: &str) -> Vec<Mask> {
    let bytes = pattern.as_bytes();
    let mut masks = Vec::new();
    let mut index = 0;
    while index < bytes.len() {
        match bytes[index] {
            b'.' => {
                masks.push(0b1111);
                index += 1;
            }
            b'[' => {
                let close = index + bytes[index..].iter().position(|&c| c == b']').expect("closed class");
                let mut mask = 0;
                for &letter in &bytes[index + 1..close] {
                    mask |= mask_of(crate::code::base_index(letter).expect("base letter"));
                }
                masks.push(mask);
                index = close + 1;
            }
            letter => {
                masks.push(mask_of(crate::code::base_index(letter).expect("base letter")));
                index += 1;
            }
        }
    }
    masks
}

/// The motif masks, for callers that want to avoid creating a splice site rather than
/// count the ones already there.
pub fn motif_masks() -> Vec<Vec<u8>> {
    DONORS
        .iter()
        .chain(ACCEPTORS.iter())
        .map(|p| parse_pattern(p))
        .collect()
}

struct Motif {
    masks: Vec<Mask>,
    donor: bool,
}

pub struct SpliceSites {
    motifs: Vec<Motif>,
    longest: usize,
    total: u32,
}

impl Default for SpliceSites {
    fn default() -> Self {
        Self::new()
    }
}

impl SpliceSites {
    pub fn new() -> Self {
        let mut motifs: Vec<Motif> = Vec::new();
        for pattern in DONORS {
            let masks = parse_pattern(pattern);
            debug_assert_eq!(masks[2], mask_of(2), "every donor has G at offset 2");
            debug_assert_eq!(masks[3], mask_of(3), "every donor has T at offset 3");
            motifs.push(Motif { masks, donor: true });
        }
        for pattern in ACCEPTORS {
            let masks = parse_pattern(pattern);
            debug_assert_eq!(masks[6], mask_of(1), "every acceptor has C at offset 6");
            debug_assert_eq!(masks[7], mask_of(0), "every acceptor has A at offset 7");
            debug_assert_eq!(masks[8], mask_of(2), "every acceptor has G at offset 8");
            motifs.push(Motif { masks, donor: false });
        }
        let longest = motifs.iter().map(|m| m.masks.len()).max().unwrap_or(0);
        Self { motifs, longest, total: 0 }
    }

    /// Reject positions that cannot match any motif, with one or three comparisons.
    fn prefilter(view: &View<'_>, at: usize, n: usize) -> (bool, bool) {
        let donor = at + 8 <= n && view.base(at + 2) == 2 && view.base(at + 3) == 3;
        let acceptor = at + 10 <= n
            && view.base(at + 6) == 1
            && view.base(at + 7) == 0
            && view.base(at + 8) == 2;
        (donor, acceptor)
    }

    fn matches_at(&self, view: &View<'_>, motif: &Motif, at: usize) -> bool {
        motif
            .masks
            .iter()
            .enumerate()
            .all(|(k, &mask)| mask & mask_of(view.base(at + k)) != 0)
    }

    fn count_in(&self, view: &View<'_>, range: Range<usize>) -> u32 {
        let n = view.n_bases();
        let mut total = 0;
        for at in range.start..range.end.min(n) {
            let (donor, acceptor) = Self::prefilter(view, at, n);
            if !donor && !acceptor {
                continue;
            }
            for motif in &self.motifs {
                if (motif.donor && !donor) || (!motif.donor && !acceptor) {
                    continue;
                }
                if at + motif.masks.len() <= n && self.matches_at(view, motif, at) {
                    total += 1;
                }
            }
        }
        total
    }
}

impl Constraint for SpliceSites {
    fn kind(&self) -> Kind {
        Kind::Splice
    }

    fn rebuild(&mut self, state: &SeqState) {
        self.total = self.count_in(&state.view(), 0..state.n_bases());
    }

    fn total(&self) -> u32 {
        self.total
    }

    fn delta(&self, before: &View<'_>, after: &View<'_>) -> i32 {
        window_delta(before, after, self.longest, |v, r| self.count_in(v, r))
    }

    fn apply(&mut self, before: &View<'_>, after: &View<'_>) {
        self.total = (self.total as i32 + self.delta(before, after)) as u32;
    }

    fn report(&self, state: &SeqState) -> Vec<Violation> {
        let view = state.view();
        let n = state.n_bases();
        let mut out = Vec::new();
        for at in 0..n {
            let (donor, acceptor) = Self::prefilter(&view, at, n);
            if !donor && !acceptor {
                continue;
            }
            for motif in &self.motifs {
                if (motif.donor && !donor) || (!motif.donor && !acceptor) {
                    continue;
                }
                if at + motif.masks.len() <= n && self.matches_at(&view, motif, at) {
                    let len = motif.masks.len();
                    let text: String = (0..len)
                        .map(|k| crate::code::BASE_LETTER[view.base(at + k) as usize] as char)
                        .collect();
                    out.push(Violation {
                        kind: Kind::Splice,
                        at: at as u32,
                        len: len as u32,
                        detail: format!(
                            "{} {text} at base {}",
                            if motif.donor { "donor" } else { "acceptor" },
                            at + 1
                        ),
                    });
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::constraints::harness;

    fn make() -> SpliceSites {
        SpliceSites::new()
    }

    /// Brute force: every pattern against every position, no prefilter.
    fn oracle(state: &SeqState) -> u32 {
        let view = state.view();
        let n = state.n_bases();
        let patterns: Vec<Vec<Mask>> = DONORS
            .iter()
            .chain(ACCEPTORS.iter())
            .map(|p| parse_pattern(p))
            .collect();
        let mut total = 0;
        for at in 0..n {
            for masks in &patterns {
                if at + masks.len() <= n
                    && masks
                        .iter()
                        .enumerate()
                        .all(|(k, &m)| m & mask_of(view.base(at + k)) != 0)
                {
                    total += 1;
                }
            }
        }
        total
    }

    #[test]
    fn matches_the_oracle() {
        harness::assert_matches_oracle(make, oracle, 200);
    }

    #[test]
    fn incremental_matches_rebuild() {
        harness::assert_incremental_matches_rebuild(make, 20_000);
    }

    #[test]
    fn agrees_with_the_python_patterns_on_known_strings() {
        // Verified against `seq_opt._splice_donors` / `_splice_acceptors` directly.
        // Note the last three: they appear in the Python's own test as "known" sites
        // but match no pattern, so that test asserted nothing about them.
        let cases: [(&str, bool, bool); 14] = [
            ("CGGTGAAT", true, false),
            ("AAGTACGG", true, false),
            ("CTGTGAGT", true, false),
            ("TAGTGTGT", true, false),
            ("TGGTCAGT", true, false),
            ("AGGTGAAA", true, false),
            ("GGGTGAGC", true, false),
            ("CTGCTACAGC", false, true),
            ("GACCCACAGG", false, true),
            ("CTGATTCAGG", false, true),
            ("CTCTTGCAGT", false, true),
            ("AGGTACCT", false, false),
            ("TGGTTAGA", false, false),
            ("TCTGGTTAGA", false, false),
        ];
        for (text, want_donor, want_acceptor) in cases {
            let bases: Vec<Base> = text
                .bytes()
                .map(|b| crate::code::base_index(b).unwrap())
                .collect();
            let hit = |patterns: &[&str]| {
                patterns.iter().map(|p| parse_pattern(p)).any(|masks| {
                    (0..bases.len().saturating_sub(masks.len()) + 1).any(|offset| {
                        masks.len() + offset <= bases.len()
                            && masks
                                .iter()
                                .enumerate()
                                .all(|(k, &m)| m & mask_of(bases[offset + k]) != 0)
                    })
                })
            };
            assert_eq!(hit(&DONORS), want_donor, "{text} as a donor");
            assert_eq!(hit(&ACCEPTORS), want_acceptor, "{text} as an acceptor");
        }
    }

    #[test]
    fn pattern_parsing_handles_classes() {
        assert_eq!(parse_pattern("A.[GT]"), vec![0b0001, 0b1111, 0b1100]);
    }
}
