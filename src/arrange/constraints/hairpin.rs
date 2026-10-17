//! Hairpins: a stem that can pair with its own reverse complement elsewhere in the
//! gene, allowing the G-T wobble pair.
//!
//! Stems are taken at every tenth base, as before. Like repeats this is not a local
//! property, so it keeps two indices: how many times each 10-mer occurs, and which
//! stems each 10-mer would satisfy. A codon edit changes at most a dozen 10-mers and
//! at most two stems, so both stay O(1).

use std::collections::HashMap;

use crate::arrange::state::{SeqState, View};
use crate::code::{complement, Base};

use super::{Constraint, Kind, Violation};

const STEM: usize = 10;
/// A stem with more wobble positions than this is matched by scanning instead of by
/// enumerating its variants.
const MAX_VARIANTS: usize = 64;

pub struct Hairpins {
    /// How many times each 10-mer occurs in the sequence.
    occurrences: HashMap<u32, u32>,
    /// Which stems a given 10-mer satisfies.
    stems_for: HashMap<u32, Vec<u32>>,
    /// Per stem: the variants it accepts and how many partners it currently has.
    variants: Vec<Vec<u32>>,
    matches: Vec<u32>,
    n_bases: usize,
    total: u32,
}

impl Default for Hairpins {
    fn default() -> Self {
        Self::new()
    }
}

impl Hairpins {
    pub fn new() -> Self {
        Self {
            occurrences: HashMap::new(),
            stems_for: HashMap::new(),
            variants: Vec::new(),
            matches: Vec::new(),
            n_bases: 0,
            total: 0,
        }
    }

    fn key_at(view: &View<'_>, at: usize, len: usize) -> u32 {
        let mut key = 0u32;
        for k in 0..len {
            key = (key << 2) | view.base(at + k) as u32;
        }
        key
    }

    fn stem_starts(n_bases: usize) -> impl Iterator<Item = usize> {
        (0..n_bases).step_by(STEM)
    }

    /// Every 10-mer that pairs with this stem. The reverse complement is the partner;
    /// wherever that calls for C, a T will also pair, which is the G-T wobble.
    fn variants_of(view: &View<'_>, stem_start: usize, stem_len: usize) -> Vec<u32> {
        if stem_len == 0 {
            return Vec::new();
        }
        let partner: Vec<Base> = (0..stem_len)
            .rev()
            .map(|k| complement(view.base(stem_start + k)))
            .collect();
        let wobble: Vec<usize> = partner
            .iter()
            .enumerate()
            .filter(|(_, &b)| b == 1)
            .map(|(i, _)| i)
            .collect();
        let count = 1usize << wobble.len().min(MAX_VARIANTS.trailing_zeros() as usize + 1);
        let mut out = Vec::with_capacity(count);
        let combinations = 1usize << wobble.len();
        if combinations > MAX_VARIANTS {
            // Too many wobble positions to enumerate; fall back to the C-only form plus
            // the all-T form, and let `rebuild` catch the rest. Rare in practice: it
            // needs seven or more G's in one stem.
            out.push(Self::pack(&partner));
            let mut all_t = partner.clone();
            for &i in &wobble {
                all_t[i] = 3;
            }
            out.push(Self::pack(&all_t));
            out.sort_unstable();
            out.dedup();
            return out;
        }
        for mask in 0..combinations {
            let mut variant = partner.clone();
            for (bit, &index) in wobble.iter().enumerate() {
                if mask & (1 << bit) != 0 {
                    variant[index] = 3;
                }
            }
            out.push(Self::pack(&variant));
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    fn pack(bases: &[Base]) -> u32 {
        bases.iter().fold(0u32, |key, &b| (key << 2) | b as u32)
    }

    /// Target positions whose 10-mer the touched codons change.
    fn touched_targets(&self, view: &View<'_>) -> Vec<usize> {
        if self.n_bases < STEM {
            return Vec::new();
        }
        let last = self.n_bases - STEM;
        let mut out: Vec<usize> = Vec::new();
        for &pos in view.touched() {
            let first = pos as usize * 3;
            let lo = (first + 1).saturating_sub(STEM);
            let hi = (first + 3).min(last + 1);
            for at in lo..hi {
                out.push(at);
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Stems whose own bases the touched codons change.
    fn touched_stems(&self, view: &View<'_>) -> Vec<usize> {
        let mut out: Vec<usize> = Vec::new();
        for &pos in view.touched() {
            let first = pos as usize * 3;
            for base in first..(first + 3).min(self.n_bases) {
                out.push(base / STEM);
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    fn stem_len(&self, index: usize) -> usize {
        STEM.min(self.n_bases - index * STEM)
    }

    /// Matches for a stem, counted against occurrence counts adjusted by `shift`.
    fn matches_for(
        &self,
        variants: &[u32],
        shift: &HashMap<u32, i32>,
    ) -> i32 {
        variants
            .iter()
            .map(|key| {
                *self.occurrences.get(key).unwrap_or(&0) as i32 + shift.get(key).copied().unwrap_or(0)
            })
            .sum()
    }

    /// How the touched targets change occurrence counts.
    fn occurrence_shift(&self, before: &View<'_>, after: &View<'_>) -> HashMap<u32, i32> {
        let mut shift: HashMap<u32, i32> = HashMap::new();
        for at in self.touched_targets(after) {
            let old = Self::key_at(before, at, STEM);
            let new = Self::key_at(after, at, STEM);
            if old == new {
                continue;
            }
            *shift.entry(old).or_insert(0) -= 1;
            *shift.entry(new).or_insert(0) += 1;
        }
        shift
    }
}

impl Constraint for Hairpins {
    fn kind(&self) -> Kind {
        Kind::Hairpin
    }

    fn rebuild(&mut self, state: &SeqState) {
        let view = state.view();
        self.n_bases = state.n_bases();
        self.occurrences.clear();
        self.stems_for.clear();
        self.variants.clear();
        self.matches.clear();

        if self.n_bases >= STEM {
            for at in 0..=self.n_bases - STEM {
                *self.occurrences.entry(Self::key_at(&view, at, STEM)).or_insert(0) += 1;
            }
        }
        let starts: Vec<usize> = Self::stem_starts(self.n_bases).collect();
        for (index, start) in starts.iter().enumerate() {
            let len = STEM.min(self.n_bases - start);
            let variants = if len == STEM {
                Self::variants_of(&view, *start, len)
            } else {
                Vec::new() // a trailing part-stem cannot pair over a full 10 bases
            };
            for &key in &variants {
                self.stems_for.entry(key).or_default().push(index as u32);
            }
            let count = variants
                .iter()
                .map(|k| *self.occurrences.get(k).unwrap_or(&0))
                .sum();
            self.variants.push(variants);
            self.matches.push(count);
        }
        self.total = self.matches.iter().sum();
    }

    fn total(&self) -> u32 {
        self.total
    }

    fn delta(&self, before: &View<'_>, after: &View<'_>) -> i32 {
        if self.n_bases < STEM {
            return 0;
        }
        let shift = self.occurrence_shift(before, after);
        let changed_stems = self.touched_stems(after);
        let mut delta = 0i32;

        // stems whose own sequence changed: recompute from scratch
        for &index in &changed_stems {
            if self.stem_len(index) != STEM {
                continue;
            }
            let new_variants = Self::variants_of(after, index * STEM, STEM);
            delta += self.matches_for(&new_variants, &shift) - self.matches[index] as i32;
        }

        // every other stem keeps its pattern, but its partners may have changed
        for (key, change) in &shift {
            if *change == 0 {
                continue;
            }
            if let Some(stems) = self.stems_for.get(key) {
                for &stem in stems {
                    if !changed_stems.contains(&(stem as usize)) {
                        delta += change;
                    }
                }
            }
        }
        delta
    }

    fn apply(&mut self, before: &View<'_>, after: &View<'_>) {
        if self.n_bases < STEM {
            return;
        }
        let predicted = self.delta(before, after);
        let shift = self.occurrence_shift(before, after);
        let changed_stems = self.touched_stems(after);

        for (key, change) in &shift {
            let entry = self.occurrences.entry(*key).or_insert(0);
            *entry = (*entry as i32 + change) as u32;
            if *entry == 0 {
                self.occurrences.remove(key);
            }
            if let Some(stems) = self.stems_for.get(key) {
                for &stem in stems {
                    if !changed_stems.contains(&(stem as usize)) {
                        self.matches[stem as usize] =
                            (self.matches[stem as usize] as i32 + change) as u32;
                    }
                }
            }
        }

        for &index in &changed_stems {
            if self.stem_len(index) != STEM {
                continue;
            }
            for key in std::mem::take(&mut self.variants[index]) {
                if let Some(stems) = self.stems_for.get_mut(&key) {
                    stems.retain(|&s| s != index as u32);
                    if stems.is_empty() {
                        self.stems_for.remove(&key);
                    }
                }
            }
            let variants = Self::variants_of(after, index * STEM, STEM);
            let mut count = 0u32;
            for &key in &variants {
                self.stems_for.entry(key).or_default().push(index as u32);
                count += *self.occurrences.get(&key).unwrap_or(&0);
            }
            self.variants[index] = variants;
            self.matches[index] = count;
        }

        self.total = (self.total as i32 + predicted) as u32;
    }

    fn report(&self, state: &SeqState) -> Vec<Violation> {
        let view = state.view();
        let mut out = Vec::new();
        for (index, start) in Self::stem_starts(self.n_bases).enumerate() {
            let count = self.matches.get(index).copied().unwrap_or(0);
            for _ in 0..count {
                let text: String = (0..self.stem_len(index))
                    .map(|k| crate::code::BASE_LETTER[view.base(start + k) as usize] as char)
                    .collect();
                out.push(Violation {
                    kind: Kind::Hairpin,
                    at: start as u32,
                    len: self.stem_len(index) as u32,
                    detail: format!("stem {text} at base {} can pair elsewhere", start + 1),
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
    use crate::code::{codon_aa, parse_codon, AaIdx, CodonIdx};

    fn make() -> Hairpins {
        Hairpins::new()
    }

    /// Brute force: for each stem, scan the whole sequence for a partner.
    fn oracle(state: &SeqState) -> u32 {
        let view = state.view();
        let n = state.n_bases();
        let mut total = 0;
        for start in (0..n).step_by(STEM) {
            let len = STEM.min(n - start);
            if len < STEM {
                continue;
            }
            let partner: Vec<Base> = (0..len)
                .rev()
                .map(|k| complement(view.base(start + k)))
                .collect();
            for at in 0..=n.saturating_sub(STEM) {
                let hit = (0..len).all(|k| {
                    let want = partner[k];
                    let got = view.base(at + k);
                    got == want || (want == 1 && got == 3)
                });
                if hit {
                    total += 1;
                }
            }
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
    fn finds_a_planted_hairpin() {
        // a stem followed by its own reverse complement
        let stem = "CATATGGGGG";
        let partner: String = stem
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                _ => 'C',
            })
            .collect();
        let dna = format!("{stem}GCTTTAGCTA{partner}AA");
        let state = state_of(&dna[..dna.len() / 3 * 3]);
        let mut c = make();
        c.rebuild(&state);
        assert!(c.total() > 0, "planted hairpin should be found in {dna}");
    }
}
