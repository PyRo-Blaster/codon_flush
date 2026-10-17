//! Homopolymers: runs of the same base.
//!
//! The original counted how often the commonest base appeared in each 6-nt window, so
//! `GCGCGC` counted as four Gs and Cs and tripped a limit of 4 without containing a run
//! at all. This counts actual runs, which is what the docstring always claimed.

use std::ops::Range;

use crate::arrange::state::{Ranges, SeqState, View};

use super::{Constraint, Kind, Violation};

pub struct Homopolymers {
    max_run: usize,
    total: u32,
}

impl Homopolymers {
    pub fn new(max_run: usize) -> Self {
        Self { max_run: max_run.max(1), total: 0 }
    }

    /// Grow a range until both ends sit on a run boundary, so every run that overlaps
    /// it is wholly inside it and can be measured.
    fn to_run_boundaries(view: &View<'_>, range: Range<usize>) -> Range<usize> {
        let n = view.n_bases();
        let mut start = range.start.min(n);
        let mut end = range.end.min(n);
        while start > 0 && view.base(start - 1) == view.base(start) {
            start -= 1;
        }
        while end < n && end > 0 && view.base(end) == view.base(end - 1) {
            end += 1;
        }
        start..end
    }

    /// Runs longer than the limit whose first base lies in `range`.
    fn count_in(&self, view: &View<'_>, range: Range<usize>) -> u32 {
        let n = view.n_bases();
        let mut total = 0;
        let mut index = range.start;
        while index < range.end.min(n) {
            let base = view.base(index);
            let mut end = index + 1;
            while end < n && view.base(end) == base {
                end += 1;
            }
            // only count a run anchored at its first base, and only if it starts here
            let starts_here = index == 0 || view.base(index - 1) != base;
            if starts_here && end - index > self.max_run {
                total += 1;
            }
            index = end;
        }
        total
    }

    /// The region a move can change, expanded to run boundaries in both views.
    fn regions(&self, before: &View<'_>, after: &View<'_>) -> Ranges {
        let mut ranges = Ranges::new();
        for &(start, end) in after.affected(self.max_run + 1).iter() {
            let a = Self::to_run_boundaries(before, start..end);
            let b = Self::to_run_boundaries(after, start..end);
            let merged = (a.start.min(b.start), a.end.max(b.end));
            match ranges.last_mut() {
                Some(last) if merged.0 <= last.1 => last.1 = last.1.max(merged.1),
                _ => ranges.push(merged),
            }
        }
        ranges
    }
}

impl Constraint for Homopolymers {
    fn kind(&self) -> Kind {
        Kind::Homopolymer
    }

    fn rebuild(&mut self, state: &SeqState) {
        let view = state.view();
        self.total = self.count_in(&view, 0..state.n_bases());
    }

    fn total(&self) -> u32 {
        self.total
    }

    fn delta(&self, before: &View<'_>, after: &View<'_>) -> i32 {
        let mut delta = 0i32;
        for &(start, end) in self.regions(before, after).iter() {
            delta += self.count_in(after, start..end) as i32;
            delta -= self.count_in(before, start..end) as i32;
        }
        delta
    }

    fn apply(&mut self, before: &View<'_>, after: &View<'_>) {
        self.total = (self.total as i32 + self.delta(before, after)) as u32;
    }

    fn report(&self, state: &SeqState) -> Vec<Violation> {
        let view = state.view();
        let n = state.n_bases();
        let mut out = Vec::new();
        let mut index = 0;
        while index < n {
            let base = view.base(index);
            let mut end = index + 1;
            while end < n && view.base(end) == base {
                end += 1;
            }
            if end - index > self.max_run {
                out.push(Violation {
                    kind: Kind::Homopolymer,
                    at: index as u32,
                    len: (end - index) as u32,
                    detail: format!(
                        "{} consecutive {}, limit is {}",
                        end - index,
                        crate::code::BASE_LETTER[base as usize] as char,
                        self.max_run
                    ),
                });
            }
            index = end;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::constraints::harness;
    use crate::code::{aa_index, codon_aa, parse_codon, AaIdx, CodonIdx};

    const LIMIT: usize = 4;

    fn make() -> Homopolymers {
        Homopolymers::new(LIMIT)
    }

    fn oracle(state: &SeqState) -> u32 {
        let dna = state.dna();
        let bytes = dna.as_bytes();
        let mut total = 0;
        let mut i = 0;
        while i < bytes.len() {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j] == bytes[i] {
                j += 1;
            }
            if j - i > LIMIT {
                total += 1;
            }
            i = j;
        }
        total
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
    fn matches_the_oracle() {
        harness::assert_matches_oracle(make, oracle, 200);
    }

    #[test]
    fn incremental_matches_rebuild() {
        harness::assert_incremental_matches_rebuild(make, 20_000);
    }

    #[test]
    fn counts_runs_not_window_composition() {
        // the case the old rule got wrong: alternating bases are not a homopolymer
        let mut c = make();
        c.rebuild(&state_of("GCGCGCGCG"));
        assert_eq!(c.total(), 0, "GCGCGC has no run at all");

        // five in a row is one violation, not several overlapping windows
        let mut c = make();
        let state = state_of("AAAAACCCG");
        c.rebuild(&state);
        assert_eq!(c.total(), 1);
        let report = c.report(&state);
        assert_eq!(report[0].at, 0);
        assert_eq!(report[0].len, 5);

        // exactly at the limit is allowed
        let mut c = make();
        c.rebuild(&state_of("AAAACCCGG"));
        assert_eq!(c.total(), 0);

        // one long run counts once, however long it is
        let mut c = make();
        c.rebuild(&state_of("AAAAAAAAAAAA"));
        assert_eq!(c.total(), 1);
        let _ = aa_index(b'K');
    }
}
