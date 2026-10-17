//! Historical experiment module, not linked into the production engine. Only proposal
//! scoring was optimized; accepted updates still move sorted-vector elements. The
//! original module comments below describe the baseline's intended design, not a
//! verified complexity guarantee. See the accompanying technical design for bounds.
//!
//! Repeated stretches: three identical codons in a row, and any 9-mer that occurs
//! twice far enough apart to be two separate copies.
//!
//! The 9-mer rule is not local: changing one codon can resolve a repeat whose other
//! copy is thousands of bases away. So instead of scanning, this keeps an index from
//! k-mer to the positions holding it. A codon edit touches at most nine k-mers, so both
//! the delta and the update are O(1) however long the gene is.

use std::collections::HashMap;
use std::ops::Range;

use crate::arrange::state::{SeqState, View};

use super::{Constraint, Kind, Violation};

const K: usize = 9;

/// Two copies closer than this overlap, and the original's non-overlapping match count
/// would not have seen them as a repeat either.
const MIN_GAP: usize = K;

pub struct Repeats {
    /// Sorted positions per packed k-mer.
    index: HashMap<u32, Vec<u32>>,
    /// Number of positions that belong to a repeated k-mer. Counting every copy, not
    /// just the first, lets the search fix whichever copy is easiest to move.
    repeated_positions: u32,
    /// Number of runs of three identical codons.
    triples: u32,
    n_bases: usize,
}

impl Default for Repeats {
    fn default() -> Self {
        Self::new()
    }
}

impl Repeats {
    pub fn new() -> Self {
        Self {
            index: HashMap::new(),
            repeated_positions: 0,
            triples: 0,
            n_bases: 0,
        }
    }

    fn key_at(view: &View<'_>, at: usize) -> u32 {
        let mut key = 0u32;
        for k in 0..K {
            key = (key << 2) | view.base(at + k) as u32;
        }
        key
    }

    fn is_repeated(positions: &[u32]) -> bool {
        positions.len() > 1
            && (*positions.last().unwrap() - positions[0]) as usize >= MIN_GAP
    }

    /// How many violations a key contributes: one per copy once it counts as repeated.
    fn contribution(positions: &[u32]) -> u32 {
        if Self::is_repeated(positions) {
            positions.len() as u32
        } else {
            0
        }
    }

    /// k-mer start positions whose content the touched codons change.
    fn touched_kmers(&self, view: &View<'_>) -> Vec<usize> {
        let last_start = self.n_bases.saturating_sub(K);
        let mut out: Vec<usize> = Vec::new();
        for &pos in view.touched() {
            let first = pos as usize * 3;
            let lo = (first + 1).saturating_sub(K);
            let hi = (first + 3).min(last_start + 1);
            for at in lo..hi {
                out.push(at);
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    fn insert(&mut self, key: u32, at: u32) {
        let positions = self.index.entry(key).or_default();
        let was = Self::contribution(positions);
        let slot = positions.partition_point(|&p| p < at);
        positions.insert(slot, at);
        let now = Self::contribution(positions);
        self.repeated_positions = self.repeated_positions + now - was;
    }

    fn remove(&mut self, key: u32, at: u32) {
        let Some(positions) = self.index.get_mut(&key) else {
            return;
        };
        let was = Self::contribution(positions);
        if let Ok(slot) = positions.binary_search(&at) {
            positions.remove(slot);
        }
        let now = Self::contribution(positions);
        self.repeated_positions = self.repeated_positions + now - was;
        if positions.is_empty() {
            self.index.remove(&key);
        }
    }

    /// Change in repeated-key count if the touched k-mers took their new values.
    fn key_delta(&self, before: &View<'_>, after: &View<'_>) -> i32 {
        let touched = self.touched_kmers(after);
        if touched.is_empty() {
            return 0;
        }
        // Work out, per key, which positions leave and which arrive.
        let mut changes: HashMap<u32, Vec<(u32, bool)>> = HashMap::new();
        for at in touched {
            let old = Self::key_at(before, at);
            let new = Self::key_at(after, at);
            if old == new {
                continue;
            }
            changes.entry(old).or_default().push((at as u32, false));
            changes.entry(new).or_default().push((at as u32, true));
        }
        let mut delta = 0i32;
        for (key, edits) in changes {
            let empty: Vec<u32> = Vec::new();
            let current = self.index.get(&key).unwrap_or(&empty);
            let was = Self::contribution(current);
            // Only count and extrema determine the repeated-key contribution.
            // At most 22 starts are touched by a two-codon edit, regardless of N.
            let removed = |p: &u32| edits.iter().any(|&(e, arriving)| !arriving && e == *p);
            let leaving = edits.iter().filter(|(_, arriving)| !arriving).count();
            let arriving = edits.len() - leaving;
            let next_len = current.len() - leaving + arriving;
            let mut lo = current.iter().find(|p| !removed(p)).copied();
            let mut hi = current.iter().rev().find(|p| !removed(p)).copied();
            for &(at, arriving) in &edits {
                if arriving {
                    lo = Some(lo.map_or(at, |p| p.min(at)));
                    hi = Some(hi.map_or(at, |p| p.max(at)));
                }
            }
            let now = match (lo, hi) {
                (Some(lo), Some(hi)) if hi - lo >= MIN_GAP as u32 => next_len as u32,
                _ => 0,
            };
            delta += now as i32 - was as i32;
        }
        delta
    }

    fn triples_in(&self, view: &View<'_>, range: Range<usize>) -> u32 {
        let n = view.n_codons();
        let mut total = 0;
        for first in range.start..range.end.min(n.saturating_sub(2)) {
            let codon = view.codon(first as u32);
            if view.codon(first as u32 + 1) == codon && view.codon(first as u32 + 2) == codon {
                total += 1;
            }
        }
        total
    }

    /// Codon positions where a triple could start or stop existing.
    fn triple_ranges(view: &View<'_>) -> Vec<Range<usize>> {
        let mut ranges: Vec<Range<usize>> = Vec::new();
        for &pos in view.touched() {
            let lo = (pos as usize).saturating_sub(2);
            let hi = pos as usize + 1;
            match ranges.last_mut() {
                Some(last) if lo <= last.end => last.end = last.end.max(hi),
                _ => ranges.push(lo..hi),
            }
        }
        ranges
    }
}

impl Constraint for Repeats {
    fn kind(&self) -> Kind {
        Kind::Repeat
    }

    fn rebuild(&mut self, state: &SeqState) {
        let view = state.view();
        self.index.clear();
        self.repeated_positions = 0;
        self.triples = 0;
        self.n_bases = state.n_bases();
        if self.n_bases >= K {
            for at in 0..=self.n_bases - K {
                let key = Self::key_at(&view, at);
                self.insert(key, at as u32);
            }
        }
        self.triples = self.triples_in(&view, 0..state.n_codons());
    }

    fn total(&self) -> u32 {
        self.repeated_positions + self.triples
    }

    fn delta(&self, before: &View<'_>, after: &View<'_>) -> i32 {
        let mut delta = self.key_delta(before, after);
        for range in Self::triple_ranges(after) {
            delta += self.triples_in(after, range.clone()) as i32;
            delta -= self.triples_in(before, range) as i32;
        }
        delta
    }

    fn apply(&mut self, before: &View<'_>, after: &View<'_>) {
        for range in Self::triple_ranges(after) {
            let gained = self.triples_in(after, range.clone()) as i32;
            let lost = self.triples_in(before, range) as i32;
            self.triples = (self.triples as i32 + gained - lost) as u32;
        }
        for at in self.touched_kmers(after) {
            let old = Self::key_at(before, at);
            let new = Self::key_at(after, at);
            if old != new {
                self.remove(old, at as u32);
                self.insert(new, at as u32);
            }
        }
    }

    fn report(&self, state: &SeqState) -> Vec<Violation> {
        let view = state.view();
        let mut out = Vec::new();
        for first in 0..state.n_codons().saturating_sub(2) {
            let codon = view.codon(first as u32);
            if view.codon(first as u32 + 1) == codon && view.codon(first as u32 + 2) == codon {
                out.push(Violation {
                    kind: Kind::Repeat,
                    at: (first * 3) as u32,
                    len: 9,
                    detail: format!(
                        "codons {}-{} are all {}",
                        first + 1,
                        first + 3,
                        crate::code::codon_string(codon)
                    ),
                });
            }
        }
        // No sort here: `ConstraintSet::report` orders the combined list, and sorting
        // an index the size of the gene on every refresh was the search's main cost.
        for positions in self.index.values() {
            if !Self::is_repeated(positions) {
                continue;
            }
            let text: String = (0..K)
                .map(|k| {
                    crate::code::BASE_LETTER[view.base(positions[0] as usize + k) as usize] as char
                })
                .collect();
            let where_else = positions
                .iter()
                .map(|p| (p + 1).to_string())
                .collect::<Vec<_>>()
                .join(", ");
            for &at in positions {
                out.push(Violation {
                    kind: Kind::Repeat,
                    at,
                    len: K as u32,
                    detail: format!("{text} occurs at bases {where_else}"),
                });
            }
        }
        out
    }
}
