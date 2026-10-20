//! RNA splice donor and acceptor motifs.
//!
//! The 20 donor and 10 acceptor patterns are the ones the Python carries, transcribed
//! as base-class masks instead of regular expressions. Every donor has `GT` fixed at
//! offset 2 and every acceptor has `CAG` at offset 6, which gives a cheap prefilter:
//! almost every position is rejected by one comparison.

use std::ops::Range;

use crate::arrange::state::{SeqState, View};
use crate::code::Base;

use super::{Constraint, Kind, Violation};

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

/// What counts as a splice problem.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum Rule {
    /// No donor and no acceptor anywhere. The default, and the original meaning.
    #[default]
    NoMotif,
    /// A donor alone, or an acceptor alone, is allowed; the two together are not.
    ///
    /// Splicing needs a pair: a 5' donor and a downstream 3' acceptor define an
    /// intron, so a lone half site cannot be spliced out *of this coding sequence*.
    /// That qualification matters. A transgene is expressed from a vector, and a lone
    /// donor inside the insert can still pair with an acceptor downstream in the
    /// vector's expression cassette, which is the usual mammalian arrangement. Under
    /// this rule the engine therefore reports every site it retains, with position
    /// and kind, so the gene can be checked against the vector it will go into.
    ///
    /// The rule is unordered: any donor together with any acceptor is a violation,
    /// whichever comes first. That is strictly more conservative than requiring the
    /// acceptor to be downstream, and it keeps the count maintainable incrementally.
    NoPair,
}

/// Sites found, counted once per position per kind however many patterns match.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub donors: u32,
    pub acceptors: u32,
}

impl Counts {
    fn shifted(self, donors: i32, acceptors: i32) -> Self {
        Self {
            donors: (self.donors as i32 + donors) as u32,
            acceptors: (self.acceptors as i32 + acceptors) as u32,
        }
    }
}

/// One motif site in a finished gene, for disclosure rather than repair.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Site {
    pub at: u32,
    pub donor: bool,
    pub text: String,
}

impl Site {
    pub fn kind(&self) -> &'static str {
        if self.donor {
            "donor"
        } else {
            "acceptor"
        }
    }
}

pub struct SpliceSites {
    motifs: Vec<Motif>,
    longest: usize,
    rule: Rule,
    found: Counts,
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
        Self { motifs, longest, rule: Rule::default(), found: Counts::default() }
    }

    pub fn with_rule(rule: Rule) -> Self {
        Self { rule, ..Self::new() }
    }

    pub fn rule(&self) -> Rule {
        self.rule
    }

    /// What the rule makes of a set of counts. This is the number the search drives
    /// to zero, so it is also what validity means.
    fn objective(rule: Rule, found: Counts) -> u32 {
        match rule {
            Rule::NoMotif => found.donors + found.acceptors,
            // Removing either side of the pair clears it, and the cheaper side is
            // the smaller count, so the distance to a legal gene is the minimum.
            Rule::NoPair => found.donors.min(found.acceptors),
        }
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

    /// Whether a donor and an acceptor start here, and how long the first matching
    /// pattern of each is. A position is one site per kind however many of the 20
    /// donor or 10 acceptor patterns match it: the patterns overlap heavily, and
    /// counting each hit made a single site read as several.
    fn site_at(&self, view: &View<'_>, at: usize, n: usize) -> (Option<usize>, Option<usize>) {
        let (donor, acceptor) = Self::prefilter(view, at, n);
        if !donor && !acceptor {
            return (None, None);
        }
        let (mut as_donor, mut as_acceptor) = (None, None);
        for motif in &self.motifs {
            if motif.donor {
                if !donor || as_donor.is_some() {
                    continue;
                }
            } else if !acceptor || as_acceptor.is_some() {
                continue;
            }
            if at + motif.masks.len() <= n && self.matches_at(view, motif, at) {
                if motif.donor {
                    as_donor = Some(motif.masks.len());
                } else {
                    as_acceptor = Some(motif.masks.len());
                }
            }
        }
        (as_donor, as_acceptor)
    }

    fn counts_in(&self, view: &View<'_>, range: Range<usize>) -> Counts {
        let n = view.n_bases();
        let mut found = Counts::default();
        for at in range.start..range.end.min(n) {
            let (donor, acceptor) = self.site_at(view, at, n);
            found.donors += u32::from(donor.is_some());
            found.acceptors += u32::from(acceptor.is_some());
        }
        found
    }

    /// Change in each count over the neighbourhood a move can reach.
    fn window_counts(&self, before: &View<'_>, after: &View<'_>) -> (i32, i32) {
        let (mut donors, mut acceptors) = (0i32, 0i32);
        for &(start, end) in after.affected(self.longest).iter() {
            let gained = self.counts_in(after, start..end);
            let lost = self.counts_in(before, start..end);
            donors += gained.donors as i32 - lost.donors as i32;
            acceptors += gained.acceptors as i32 - lost.acceptors as i32;
        }
        (donors, acceptors)
    }

    /// Every site in the gene, whether or not the rule objects to it. This is what
    /// the caller discloses for a gene that keeps a site by necessity.
    pub fn sites(&self, state: &SeqState) -> Vec<Site> {
        let view = state.view();
        let n = state.n_bases();
        let mut out = Vec::new();
        for at in 0..n {
            let (donor, acceptor) = self.site_at(&view, at, n);
            for (len, is_donor) in [(donor, true), (acceptor, false)] {
                if let Some(len) = len {
                    out.push(Site {
                        at: at as u32,
                        donor: is_donor,
                        text: (0..len)
                            .map(|k| crate::code::BASE_LETTER[view.base(at + k) as usize] as char)
                            .collect(),
                    });
                }
            }
        }
        out
    }

    pub fn counts(&self) -> Counts {
        self.found
    }
}

impl Constraint for SpliceSites {
    fn kind(&self) -> Kind {
        Kind::Splice
    }

    fn rebuild(&mut self, state: &SeqState) {
        self.found = self.counts_in(&state.view(), 0..state.n_bases());
    }

    fn total(&self) -> u32 {
        Self::objective(self.rule, self.found)
    }

    fn delta(&self, before: &View<'_>, after: &View<'_>) -> i32 {
        // Under the pair rule the objective is a function of the gene's two totals,
        // not of the window alone, so the window supplies the change in each count
        // and the running totals turn that into a change in the objective. Reading
        // `self.found` is sound here because both views are taken against the
        // pre-move state, which is what those totals describe.
        let (donors, acceptors) = self.window_counts(before, after);
        Self::objective(self.rule, self.found.shifted(donors, acceptors)) as i32
            - Self::objective(self.rule, self.found) as i32
    }

    fn apply(&mut self, before: &View<'_>, after: &View<'_>) {
        let (donors, acceptors) = self.window_counts(before, after);
        self.found = self.found.shifted(donors, acceptors);
    }

    fn report(&self, state: &SeqState) -> Vec<Violation> {
        let found = self.counts_in(&state.view(), 0..state.n_bases());
        // Which sites the search is being asked to remove. Under the pair rule a
        // lone half site is not a violation, and when both kinds are present only
        // the rarer one needs to go, so that is what is reported. Ties go to the
        // acceptors: no peptide forces one, while 30 tripeptides force a donor.
        let wanted = match self.rule {
            Rule::NoMotif => None,
            Rule::NoPair if found.donors == 0 || found.acceptors == 0 => return Vec::new(),
            Rule::NoPair => Some(found.donors < found.acceptors),
        };
        self.sites(state)
            .into_iter()
            .filter(|site| wanted.map_or(true, |donor| site.donor == donor))
            .map(|site| Violation {
                kind: Kind::Splice,
                at: site.at,
                len: site.text.len() as u32,
                detail: format!("{} {} at base {}", site.kind(), site.text, site.at + 1),
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::constraints::harness;

    fn make() -> SpliceSites {
        SpliceSites::new()
    }

    fn make_pair() -> SpliceSites {
        SpliceSites::with_rule(Rule::NoPair)
    }

    /// Brute force: every pattern against every position, no prefilter, counting a
    /// position once per kind however many patterns match it.
    fn oracle_counts(state: &SeqState) -> Counts {
        let view = state.view();
        let n = state.n_bases();
        let donors: Vec<Vec<Mask>> = DONORS.iter().map(|p| parse_pattern(p)).collect();
        let acceptors: Vec<Vec<Mask>> = ACCEPTORS.iter().map(|p| parse_pattern(p)).collect();
        let hit = |patterns: &[Vec<Mask>], at: usize| {
            patterns.iter().any(|masks| {
                at + masks.len() <= n
                    && masks
                        .iter()
                        .enumerate()
                        .all(|(k, &m)| m & mask_of(view.base(at + k)) != 0)
            })
        };
        let mut found = Counts::default();
        for at in 0..n {
            found.donors += u32::from(hit(&donors, at));
            found.acceptors += u32::from(hit(&acceptors, at));
        }
        found
    }

    fn oracle(state: &SeqState) -> u32 {
        let found = oracle_counts(state);
        found.donors + found.acceptors
    }

    fn oracle_pair(state: &SeqState) -> u32 {
        let found = oracle_counts(state);
        found.donors.min(found.acceptors)
    }

    #[test]
    fn matches_the_oracle() {
        harness::assert_matches_oracle(make, oracle, 200);
    }

    #[test]
    fn the_pair_rule_matches_its_own_oracle() {
        harness::assert_matches_oracle(make_pair, oracle_pair, 200);
    }

    #[test]
    fn incremental_matches_rebuild() {
        harness::assert_incremental_matches_rebuild(make, 20_000);
    }

    #[test]
    fn the_pair_rule_is_incremental_too() {
        // The objective is a function of the whole gene's two counts, so this is
        // the test that the window delta really does reconstruct it.
        harness::assert_incremental_matches_rebuild(make_pair, 20_000);
    }

    fn state_of(dna: &str) -> SeqState {
        use crate::code::{codon_aa, parse_codon, AaIdx, CodonIdx};
        let codons: Vec<CodonIdx> = dna
            .as_bytes()
            .chunks(3)
            .map(|c| parse_codon(std::str::from_utf8(c).unwrap()).unwrap())
            .collect();
        let aa: Vec<AaIdx> = codons.iter().copied().map(codon_aa).collect();
        SeqState::new(aa, codons)
    }

    fn scored(rule: Rule, dna: &str) -> (u32, Counts) {
        let state = state_of(dna);
        let mut sites = SpliceSites::with_rule(rule);
        sites.rebuild(&state);
        (sites.total(), sites.counts())
    }

    #[test]
    fn one_site_counts_once_however_many_patterns_match_it() {
        // TGGTATGT, the donor forced by the tripeptide WYV, matches five of the 20
        // donor patterns. It is one site, and a review read the old count of five
        // as five separate problems.
        let (total, found) = scored(Rule::NoMotif, "TGGTATGTA");
        assert_eq!(found, Counts { donors: 1, acceptors: 0 });
        assert_eq!(total, 1);
    }

    #[test]
    fn the_pair_rule_allows_a_lone_donor_but_not_a_pair() {
        let (lone_donor, counts) = scored(Rule::NoPair, "TGGTATGTA");
        assert!(counts.donors > 0 && counts.acceptors == 0);
        assert_eq!(lone_donor, 0, "a donor with no acceptor cannot splice anything out");
        // the same gene under the original rule is still a violation
        assert_eq!(scored(Rule::NoMotif, "TGGTATGTA").0, 1);
    }

    #[test]
    fn the_pair_rule_counts_the_cheaper_side_of_the_pair() {
        // CTGCTACAGC is an acceptor; put one in a gene that also has donors.
        let dna = "TGGTATGTACTGCTACAGCAAA";
        let dna = &dna[..dna.len() - dna.len() % 3];
        let (total, found) = scored(Rule::NoPair, dna);
        assert!(found.donors > 0, "{found:?}");
        assert!(found.acceptors > 0, "{found:?}");
        assert_eq!(
            total,
            found.donors.min(found.acceptors),
            "the distance to a legal gene is the side that is cheaper to clear"
        );
        assert!(total > 0, "a donor and an acceptor together is a violation");
    }

    #[test]
    fn the_pair_rule_reports_only_the_sites_it_wants_removed() {
        let dna = "TGGTATGTACTGCTACAGCAAA";
        let dna = &dna[..dna.len() - dna.len() % 3];
        let state = state_of(dna);
        let mut sites = SpliceSites::with_rule(Rule::NoPair);
        sites.rebuild(&state);
        let reported = sites.report(&state);
        assert!(!reported.is_empty());
        assert!(
            reported.iter().all(|v| v.detail.starts_with("acceptor")),
            "with more donors than acceptors the acceptors are the ones to clear: {reported:?}"
        );
        // but every site is still disclosed, including the ones left alone
        let all = sites.sites(&state);
        assert!(all.iter().any(|s| s.donor), "{all:?}");
        assert!(all.iter().any(|s| !s.donor), "{all:?}");
    }

    #[test]
    fn a_lone_donor_is_disclosed_even_though_it_is_not_a_violation() {
        let state = state_of("TGGTATGTA");
        let mut sites = SpliceSites::with_rule(Rule::NoPair);
        sites.rebuild(&state);
        assert_eq!(sites.total(), 0, "nothing for the search to fix");
        assert!(sites.report(&state).is_empty());
        let disclosed = sites.sites(&state);
        assert_eq!(disclosed.len(), 1, "{disclosed:?}");
        assert_eq!(disclosed[0].kind(), "donor");
        assert_eq!(disclosed[0].at, 0);
        assert_eq!(disclosed[0].text, "TGGTATGT");
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
