//! Alternate translation starts: a ribosome binding site sitting just upstream of a
//! start codon, which would let the ribosome begin translating inside the gene.
//!
//! Two definitions are available.
//!
//! [`Rule::PurineRun`] is what the tool has always used: five purines anywhere in the
//! 15 bases ending 3 before TTG, CTG or ATG. (The original listed 32 binding-site
//! strings; they are exactly the 32 five-base strings over A and G.)
//!
//! [`Rule::ShineDalgarno`] asks whether the upstream sequence could actually pair with
//! the 3' end of the 16S rRNA, which is what recruits the ribosome. The anti-SD is
//! 3'-AUUCCUCCA-5', so the mRNA consensus is AGGAGG, and a site needs four contiguous
//! bases of that core placed 4 to 14 bases ahead of the start codon. It also uses the
//! start codons E. coli really initiates at - ATG, GTG and TTG - rather than NCBI
//! table 1's list, which includes CTG and omits GTG.

use std::ops::Range;

use crate::arrange::state::{SeqState, View};
use crate::code::{is_purine, Base};

use super::{window_delta, Constraint, Kind, Violation};

/// Bases between the start of the search window and the xTG.
const UPSTREAM: usize = 18;
/// The three bases closest to the xTG are ignored, as before.
const GAP: usize = 3;
const RUN: usize = 5;
/// An edit this far ahead of an xTG can still change its upstream window.
const MARGIN: usize = UPSTREAM + GAP;

/// Which definition of an internal start site to enforce.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rule {
    /// Five purines before TTG, CTG or ATG.
    PurineRun,
    /// A Shine-Dalgarno core before ATG, GTG or TTG.
    ShineDalgarno,
}

/// Four contiguous bases of the AGGAGG consensus.
const SD_CORES: [[Base; 4]; 3] = [
    [0, 2, 2, 0], // AGGA
    [2, 2, 0, 2], // GGAG
    [2, 0, 2, 2], // GAGG
];

/// Bases allowed between the end of the core and the start codon.
const SPACER: (usize, usize) = (4, 14);

pub struct StartSites {
    rule: Rule,
    total: u32,
}

impl Default for StartSites {
    fn default() -> Self {
        Self::new()
    }
}

impl StartSites {
    pub fn new() -> Self {
        Self::with_rule(Rule::ShineDalgarno)
    }

    pub fn with_rule(rule: Rule) -> Self {
        Self { rule, total: 0 }
    }

    fn is_start_codon(&self, view: &View<'_>, at: usize) -> bool {
        // any offset, not just in frame: a ribosome does not know the reading frame
        let (b0, b1, b2) = (view.base(at), view.base(at + 1), view.base(at + 2));
        match self.rule {
            // TTG, CTG, ATG
            Rule::PurineRun => b1 == 3 && b2 == 2 && (b0 == 3 || b0 == 1 || b0 == 0),
            // ATG, GTG, TTG
            Rule::ShineDalgarno => b1 == 3 && b2 == 2 && (b0 == 0 || b0 == 2 || b0 == 3),
        }
    }

    /// Could the upstream sequence pair with the 16S rRNA well enough to recruit a
    /// ribosome, at a workable distance from the start codon?
    fn has_sd_core(view: &View<'_>, at: usize) -> bool {
        let last_core_start = at - SPACER.0 - 4;
        let first_core_start = at.saturating_sub(SPACER.1 + 4);
        for start in first_core_start..=last_core_start {
            for core in SD_CORES.iter() {
                if (0..4).all(|k| view.base(start + k) == core[k]) {
                    return true;
                }
            }
        }
        false
    }

    /// Does a run of five purines fit inside the window ending `GAP` before `at`?
    fn has_binding_site(view: &View<'_>, at: usize) -> bool {
        let window = at - UPSTREAM..at - GAP;
        let mut run = 0;
        for index in window {
            if is_purine(view.base(index)) {
                run += 1;
                if run >= RUN {
                    return true;
                }
            } else {
                run = 0;
            }
        }
        false
    }

    fn is_violation(&self, view: &View<'_>, at: usize) -> bool {
        if !self.is_start_codon(view, at) {
            return false;
        }
        match self.rule {
            Rule::PurineRun => Self::has_binding_site(view, at),
            Rule::ShineDalgarno => Self::has_sd_core(view, at),
        }
    }

    fn count_in(&self, view: &View<'_>, range: Range<usize>) -> u32 {
        let n = view.n_bases();
        let mut total = 0;
        for at in range.start.max(UPSTREAM)..range.end.min(n.saturating_sub(2)) {
            if self.is_violation(view, at) {
                total += 1;
            }
        }
        total
    }
}

impl Constraint for StartSites {
    fn kind(&self) -> Kind {
        Kind::StartSite
    }

    fn rebuild(&mut self, state: &SeqState) {
        self.total = self.count_in(&state.view(), 0..state.n_bases());
    }

    fn total(&self) -> u32 {
        self.total
    }

    fn delta(&self, before: &View<'_>, after: &View<'_>) -> i32 {
        window_delta(before, after, MARGIN, |v, r| self.count_in(v, r))
    }

    fn apply(&mut self, before: &View<'_>, after: &View<'_>) {
        self.total = (self.total as i32 + self.delta(before, after)) as u32;
    }

    fn report(&self, state: &SeqState) -> Vec<Violation> {
        let view = state.view();
        let n = state.n_bases();
        let mut out = Vec::new();
        for at in UPSTREAM..n.saturating_sub(2) {
            if self.is_violation(&view, at) {
                let codon: String = (0..3)
                    .map(|k| crate::code::BASE_LETTER[view.base(at + k) as usize] as char)
                    .collect();
                let what = match self.rule {
                    Rule::PurineRun => "purine run",
                    Rule::ShineDalgarno => "Shine-Dalgarno core",
                };
                out.push(Violation {
                    kind: Kind::StartSite,
                    at: (at - UPSTREAM) as u32,
                    len: (UPSTREAM + 3) as u32,
                    detail: format!("{what} upstream of {codon} at base {}", at + 1),
                });
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::constraints::harness;

    /// The original rule, kept for comparison and for anyone who needs the old output.
    fn make() -> StartSites {
        StartSites::with_rule(Rule::PurineRun)
    }

    fn oracle(state: &SeqState) -> u32 {
        let dna = state.dna();
        let b = dna.as_bytes();
        let mut total = 0;
        for at in UPSTREAM..b.len().saturating_sub(2) {
            let codon = &dna[at..at + 3];
            if !matches!(codon, "TTG" | "CTG" | "ATG") {
                continue;
            }
            let window = &dna[at - UPSTREAM..at - GAP];
            let purine_run = window
                .as_bytes()
                .split(|&c| c != b'A' && c != b'G')
                .map(|r| r.len())
                .max()
                .unwrap_or(0);
            if purine_run >= RUN {
                total += 1;
            }
        }
        let _ = b;
        total
    }

    fn make_sd() -> StartSites {
        StartSites::new()
    }

    fn sd_oracle(state: &SeqState) -> u32 {
        let dna = state.dna();
        let mut total = 0;
        for at in UPSTREAM..dna.len().saturating_sub(2) {
            if !matches!(&dna[at..at + 3], "ATG" | "GTG" | "TTG") {
                continue;
            }
            let window = &dna[at - SPACER.1 - 4..at - SPACER.0];
            if ["AGGA", "GGAG", "GAGG"].iter().any(|c| window.contains(c)) {
                total += 1;
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
    fn shine_dalgarno_matches_its_oracle() {
        harness::assert_matches_oracle(make_sd, sd_oracle, 200);
    }

    #[test]
    fn shine_dalgarno_incremental_matches_rebuild() {
        harness::assert_incremental_matches_rebuild(make_sd, 20_000);
    }

    #[test]
    fn the_two_rules_disagree_where_the_biology_says_they_should() {
        use crate::code::{codon_aa, parse_codon, AaIdx, CodonIdx};
        let state_of = |dna: &str| {
            let codons: Vec<CodonIdx> = dna
                .as_bytes()
                .chunks(3)
                .map(|c| parse_codon(std::str::from_utf8(c).unwrap()).unwrap())
                .collect();
            let aa: Vec<AaIdx> = codons.iter().copied().map(codon_aa).collect();
            SeqState::new(aa, codons)
        };

        // A real ribosome binding site: AGGAGG eight bases ahead of an ATG.
        let real = state_of("CCGCCGCCGAGGAGGCCGCCGATGCCG");
        let mut sd = make_sd();
        sd.rebuild(&real);
        assert_eq!(sd.total(), 1, "a genuine Shine-Dalgarno must be caught");

        // A purine run that cannot pair with the 16S rRNA: no GG anywhere.
        let not_sd = state_of("CCGCCGCCAAGAAAGACCGCCATGCCG");
        let mut sd = make_sd();
        sd.rebuild(&not_sd);
        let mut legacy = make();
        legacy.rebuild(&not_sd);
        assert_eq!(legacy.total(), 1, "the old rule flags any five purines");
        assert_eq!(sd.total(), 0, "but AAGAAAGA cannot recruit a ribosome");

        // GTG is the second most common start in E. coli and the old rule ignores it.
        let gtg = state_of("CCGCCGCCGAGGAGGCCGCCGGTGCCG");
        let mut sd = make_sd();
        sd.rebuild(&gtg);
        let mut legacy = make();
        legacy.rebuild(&gtg);
        assert_eq!(sd.total(), 1, "GTG starts translation too");
        assert_eq!(legacy.total(), 0, "the old rule never looks for GTG");
    }
}
