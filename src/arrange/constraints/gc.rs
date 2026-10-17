//! GC content in sliding windows.
//!
//! Unlike the original, which stepped a window of `3 x gene length` across the gene
//! (so the last window set was really the whole gene), these are true sliding windows
//! at every base offset. Whole-gene GC is settled in `composition`, where no
//! rearrangement can disturb it.

use crate::arrange::state::{Ranges, SeqState, View};
use crate::small::Small;

/// Base index and whether it gained (+1) or lost (-1) a G or C.
type Edits = Small<(usize, i32), 6>;

use super::{Constraint, Kind, Violation};

#[derive(Debug, Clone, Copy)]
pub struct WindowSet {
    pub name: &'static str,
    pub bases: usize,
    pub low: f64,
    pub high: f64,
}

/// The bands the tool has always used, minus the whole-gene one.
pub const DEFAULT_WINDOWS: [WindowSet; 3] = [
    WindowSet { name: "IDT", bases: 20, low: 0.15, high: 0.90 },
    WindowSet { name: "twist", bases: 50, low: 0.15, high: 0.80 },
    WindowSet { name: "IDT_long", bases: 100, low: 0.28, high: 0.68 },
];

pub struct GcWindows {
    sets: Vec<WindowSet>,
    /// GC count per window start, one vector per set.
    counts: Vec<Vec<u16>>,
    total: u32,
}

impl GcWindows {
    pub fn new(sets: Vec<WindowSet>) -> Self {
        Self { sets, counts: Vec::new(), total: 0 }
    }

    /// Window length actually used: a sequence shorter than the window gets one window
    /// covering all of it.
    fn window_len(set: &WindowSet, n_bases: usize) -> usize {
        set.bases.min(n_bases).max(1)
    }

    fn window_starts(set: &WindowSet, n_bases: usize) -> usize {
        n_bases.saturating_sub(Self::window_len(set, n_bases)) + 1
    }

    fn is_violation(set: &WindowSet, gc: u16, len: usize) -> bool {
        let frac = gc as f64 / len as f64;
        frac < set.low || frac > set.high
    }

    /// Window starts that could change when the bases in `touched` change.
    fn affected_starts(&self, set_index: usize, view: &View<'_>) -> Ranges {
        let n = view.n_bases();
        let len = Self::window_len(&self.sets[set_index], n);
        let starts = Self::window_starts(&self.sets[set_index], n);
        let mut ranges = Ranges::new();
        for &pos in view.touched() {
            let first_base = pos as usize * 3;
            let lo = (first_base + 1).saturating_sub(len);
            let hi = (first_base + 3).min(starts);
            if lo >= hi {
                continue;
            }
            match ranges.last_mut() {
                Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
                _ => ranges.push((lo, hi)),
            }
        }
        ranges
    }

    /// GC in a window read straight from a view, used for rebuilds.
    fn gc_of(view: &View<'_>, start: usize, len: usize) -> u16 {
        (start..start + len)
            .filter(|&i| crate::code::is_gc(view.base(i)))
            .count() as u16
    }

    /// Bases whose GC-ness the move changes, as (base index, +1 or -1). At most six,
    /// and often none at all when two synonymous codons share a GC pattern.
    fn edits(view: &View<'_>) -> Edits {
        let mut out = Edits::new();
        for &pos in view.touched() {
            let stored = crate::code::codon_bases(view.state().codons[pos as usize]);
            for (offset, &base) in stored.iter().enumerate() {
                let index = pos as usize * 3 + offset;
                let was = crate::code::is_gc(base) as i32;
                let now = crate::code::is_gc(view.base(index)) as i32;
                if was != now {
                    out.push((index, now - was));
                }
            }
        }
        out
    }

    fn shift_in_window(edits: &Edits, start: usize, len: usize) -> i32 {
        edits
            .iter()
            .filter(|(index, _)| *index >= start && *index < start + len)
            .map(|(_, delta)| delta)
            .sum()
    }
}

impl Constraint for GcWindows {
    fn kind(&self) -> Kind {
        Kind::Gc
    }

    fn rebuild(&mut self, state: &SeqState) {
        let view = state.view();
        let n = state.n_bases();
        self.counts.clear();
        self.total = 0;
        for set in &self.sets {
            let len = Self::window_len(set, n);
            let starts = Self::window_starts(set, n);
            let mut counts = Vec::with_capacity(starts);
            let mut running = Self::gc_of(&view, 0, len);
            counts.push(running);
            for start in 1..starts {
                if crate::code::is_gc(view.base(start - 1)) {
                    running -= 1;
                }
                if crate::code::is_gc(view.base(start + len - 1)) {
                    running += 1;
                }
                counts.push(running);
            }
            self.total += counts
                .iter()
                .filter(|&&gc| Self::is_violation(set, gc, len))
                .count() as u32;
            self.counts.push(counts);
        }
    }

    fn total(&self) -> u32 {
        self.total
    }

    fn delta(&self, _before: &View<'_>, after: &View<'_>) -> i32 {
        let edits = Self::edits(after);
        if edits.is_empty() {
            return 0; // the move does not move a single G or C
        }
        let n = after.n_bases();
        let mut delta = 0i32;
        for (index, set) in self.sets.iter().enumerate() {
            let len = Self::window_len(set, n);
            for &(from, to) in self.affected_starts(index, after).iter() {
                for start in from..to {
                    let shift = Self::shift_in_window(&edits, start, len);
                    if shift == 0 {
                        continue; // this window's GC is unchanged, so its status is too
                    }
                    let was = self.counts[index][start];
                    let now = (was as i32 + shift) as u16;
                    delta += Self::is_violation(set, now, len) as i32;
                    delta -= Self::is_violation(set, was, len) as i32;
                }
            }
        }
        delta
    }

    fn apply(&mut self, _before: &View<'_>, after: &View<'_>) {
        let edits = Self::edits(after);
        if edits.is_empty() {
            return;
        }
        let n = after.n_bases();
        for index in 0..self.sets.len() {
            let set = self.sets[index];
            let len = Self::window_len(&set, n);
            for &(from, to) in self.affected_starts(index, after).iter() {
                for start in from..to {
                    let shift = Self::shift_in_window(&edits, start, len);
                    if shift == 0 {
                        continue;
                    }
                    let was = self.counts[index][start];
                    let now = (was as i32 + shift) as u16;
                    self.total = (self.total as i32
                        + Self::is_violation(&set, now, len) as i32
                        - Self::is_violation(&set, was, len) as i32)
                        as u32;
                    self.counts[index][start] = now;
                }
            }
        }
    }

    fn report(&self, state: &SeqState) -> Vec<Violation> {
        let view = state.view();
        let n = state.n_bases();
        let mut out = Vec::new();
        for (index, set) in self.sets.iter().enumerate() {
            let len = Self::window_len(set, n);
            for start in 0..Self::window_starts(set, n) {
                let gc = if self.counts.is_empty() {
                    Self::gc_of(&view, start, len)
                } else {
                    self.counts[index][start]
                };
                if Self::is_violation(set, gc, len) {
                    out.push(Violation {
                        kind: Kind::Gc,
                        at: start as u32,
                        len: len as u32,
                        detail: format!(
                            "{}: {:.0}% GC over {} nt, outside {:.0}-{:.0}%",
                            set.name,
                            gc as f64 / len as f64 * 100.0,
                            len,
                            set.low * 100.0,
                            set.high * 100.0
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
    use crate::code::is_gc;

    fn oracle(state: &SeqState) -> u32 {
        let view = state.view();
        let n = state.n_bases();
        let mut total = 0;
        for set in DEFAULT_WINDOWS.iter() {
            let len = GcWindows::window_len(set, n);
            for start in 0..GcWindows::window_starts(set, n) {
                let gc = (start..start + len).filter(|&i| is_gc(view.base(i))).count();
                let frac = gc as f64 / len as f64;
                if frac < set.low || frac > set.high {
                    total += 1;
                }
            }
        }
        total
    }

    fn make() -> GcWindows {
        GcWindows::new(DEFAULT_WINDOWS.to_vec())
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
    fn band_edges_are_not_violations() {
        // 20-base window at exactly 15% GC: 3 of 20. The bound is strict, as before.
        let set = WindowSet { name: "t", bases: 20, low: 0.15, high: 0.90 };
        assert!(!GcWindows::is_violation(&set, 3, 20));
        assert!(GcWindows::is_violation(&set, 2, 20));
        assert!(!GcWindows::is_violation(&set, 18, 20));
        assert!(GcWindows::is_violation(&set, 19, 20));
    }
}
