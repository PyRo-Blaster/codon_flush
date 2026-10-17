//! The sequence being arranged, and the two kinds of move allowed on it.
//!
//! A [`Move::Swap`] exchanges the codons at two positions holding the same amino acid.
//! It cannot change the codon multiset, so it cannot change the deviation from the host
//! profile or the codon adaptation index. A [`Move::Subst`] replaces one codon with a
//! synonym, which does change the multiset and so is charged against the budget.
//!
//! Constraints evaluate a move through a [`View`], which layers an unapplied patch over
//! the current sequence. The same counting routine run against the patched and
//! unpatched view gives the move's exact effect, with no speculative mutation.

use crate::code::{codon_aa, codon_bases, AaIdx, Base, CodonIdx, N_AA, N_CODONS};
use crate::small::Small;

/// Base ranges, as (start, end) pairs.
pub type Ranges = Small<(usize, usize), 2>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Move {
    Swap { p: u32, q: u32 },
    Subst { p: u32, codon: CodonIdx },
}

/// The codon writes a move performs, without touching the sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Patch {
    pub first: (u32, CodonIdx),
    pub second: Option<(u32, CodonIdx)>,
}

impl Patch {
    pub fn positions(&self) -> impl Iterator<Item = u32> + '_ {
        std::iter::once(self.first.0).chain(self.second.map(|(p, _)| p))
    }

    fn codon_at(&self, pos: u32) -> Option<CodonIdx> {
        if self.first.0 == pos {
            return Some(self.first.1);
        }
        match self.second {
            Some((p, c)) if p == pos => Some(c),
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SeqState {
    pub aa: Vec<AaIdx>,
    pub codons: Vec<CodonIdx>,
    /// Bases, three per codon, kept in step with `codons`.
    pub nt: Vec<Base>,
    pub counts: [u32; N_CODONS],
    /// Positions grouped by amino acid, for finding swap partners.
    pub pos_by_aa: Vec<Vec<u32>>,
}

impl SeqState {
    pub fn new(aa: Vec<AaIdx>, codons: Vec<CodonIdx>) -> Self {
        assert_eq!(aa.len(), codons.len(), "one codon per residue");
        let mut nt = Vec::with_capacity(codons.len() * 3);
        let mut counts = [0u32; N_CODONS];
        let mut pos_by_aa = vec![Vec::new(); N_AA];
        for (i, &codon) in codons.iter().enumerate() {
            debug_assert_eq!(codon_aa(codon), aa[i], "codon must encode its residue");
            nt.extend_from_slice(&codon_bases(codon));
            counts[codon as usize] += 1;
            pos_by_aa[aa[i] as usize].push(i as u32);
        }
        Self {
            aa,
            codons,
            nt,
            counts,
            pos_by_aa,
        }
    }

    pub fn n_codons(&self) -> usize {
        self.codons.len()
    }

    pub fn n_bases(&self) -> usize {
        self.nt.len()
    }

    pub fn is_empty(&self) -> bool {
        self.codons.is_empty()
    }

    pub fn patch_of(&self, mv: Move) -> Patch {
        match mv {
            Move::Swap { p, q } => Patch {
                first: (p, self.codons[q as usize]),
                second: Some((q, self.codons[p as usize])),
            },
            Move::Subst { p, codon } => Patch {
                first: (p, codon),
                second: None,
            },
        }
    }

    pub fn apply(&mut self, mv: Move) {
        let patch = self.patch_of(mv);
        self.write(patch.first.0, patch.first.1);
        if let Some((pos, codon)) = patch.second {
            self.write(pos, codon);
        }
    }

    fn write(&mut self, pos: u32, codon: CodonIdx) {
        let old = self.codons[pos as usize];
        if old == codon {
            return;
        }
        debug_assert_eq!(codon_aa(codon), self.aa[pos as usize], "move changed the protein");
        self.counts[old as usize] -= 1;
        self.counts[codon as usize] += 1;
        self.codons[pos as usize] = codon;
        let bases = codon_bases(codon);
        let base0 = pos as usize * 3;
        self.nt[base0..base0 + 3].copy_from_slice(&bases);
    }

    /// The residues this sequence encodes, for checking that a run never changed them.
    pub fn translate(&self) -> Vec<AaIdx> {
        self.codons.iter().copied().map(codon_aa).collect()
    }

    pub fn dna(&self) -> String {
        self.nt
            .iter()
            .map(|&b| crate::code::BASE_LETTER[b as usize] as char)
            .collect()
    }

    pub fn view(&self) -> View<'_> {
        View {
            state: self,
            patch: None,
            touched: Vec::new(),
        }
    }

    /// A view of the sequence as the patch would leave it.
    pub fn patched(&self, patch: Patch) -> View<'_> {
        View {
            state: self,
            patch: Some(patch),
            touched: sorted_positions(patch),
        }
    }

    /// A view of the sequence as it is, scoped to the patch's positions, so that the
    /// same counting routine measures the same region before and after.
    pub fn scoped(&self, patch: Patch) -> View<'_> {
        View {
            state: self,
            patch: None,
            touched: sorted_positions(patch),
        }
    }
}

/// Touched positions in ascending order. Constraints merge adjacent regions as they
/// walk this list, so an unsorted pair would silently drop the lower region.
fn sorted_positions(patch: Patch) -> Vec<u32> {
    let mut positions: Vec<u32> = patch.positions().collect();
    positions.sort_unstable();
    positions.dedup();
    positions
}

pub struct View<'a> {
    state: &'a SeqState,
    patch: Option<Patch>,
    touched: Vec<u32>,
}

impl<'a> View<'a> {
    pub fn state(&self) -> &'a SeqState {
        self.state
    }

    pub fn n_codons(&self) -> usize {
        self.state.codons.len()
    }

    pub fn n_bases(&self) -> usize {
        self.state.nt.len()
    }

    pub fn codon(&self, pos: u32) -> CodonIdx {
        match self.patch.and_then(|p| p.codon_at(pos)) {
            Some(codon) => codon,
            None => self.state.codons[pos as usize],
        }
    }

    pub fn base(&self, index: usize) -> Base {
        match self.patch.and_then(|p| p.codon_at((index / 3) as u32)) {
            Some(codon) => codon_bases(codon)[index % 3],
            None => self.state.nt[index],
        }
    }

    /// Codon positions this view is scoped to. Empty for a whole-sequence view.
    pub fn touched(&self) -> &[u32] {
        &self.touched
    }

    /// The base ranges a constraint must re-examine, given how far its windows reach.
    /// `margin` is in bases. A move touches at most two codons, so there are at most
    /// two ranges and no allocation is needed.
    pub fn affected(&self, margin: usize) -> Ranges {
        let mut ranges = Ranges::new();
        for &pos in &self.touched {
            let start = (pos as usize * 3).saturating_sub(margin);
            let end = (pos as usize * 3 + 3 + margin).min(self.n_bases());
            match ranges.last_mut() {
                Some(last) if start <= last.1 => last.1 = last.1.max(end),
                _ => ranges.push((start, end)),
            }
        }
        ranges
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::{aa_index, parse_codon, N_CODONS};
    use crate::composition;
    use crate::host::HostTable;
    use crate::rng::Rng;
    use crate::testdata::{ECOLI_TABLE, GFP};
    use std::path::Path;

    fn gfp_state(seed: u64) -> (Vec<AaIdx>, SeqState) {
        let host = HostTable::load(Path::new(ECOLI_TABLE), 0.10).unwrap();
        let protein: Vec<AaIdx> = GFP.bytes().map(|b| aa_index(b).unwrap()).collect();
        let comp = composition::build(&protein, &host, 0.1, (0.30, 0.65)).unwrap();
        let mut rng = Rng::new(seed);
        let mut pools: Vec<Vec<CodonIdx>> = vec![Vec::new(); N_AA];
        for codon in 0..N_CODONS as CodonIdx {
            for _ in 0..comp.counts[codon as usize] {
                pools[codon_aa(codon) as usize].push(codon);
            }
        }
        for pool in pools.iter_mut() {
            rng.shuffle(pool);
        }
        let mut next = [0usize; N_AA];
        let codons: Vec<CodonIdx> = protein
            .iter()
            .map(|&aa| {
                let i = next[aa as usize];
                next[aa as usize] += 1;
                pools[aa as usize][i]
            })
            .collect();
        (protein.clone(), SeqState::new(protein, codons))
    }

    #[test]
    fn swaps_preserve_translation_and_multiset() {
        let (protein, mut state) = gfp_state(1);
        let counts_before = state.counts;
        let mut rng = Rng::new(9);
        for _ in 0..10_000 {
            let aa = state.aa[rng.below(state.n_codons())];
            let peers = &state.pos_by_aa[aa as usize];
            if peers.len() < 2 {
                continue;
            }
            let p = peers[rng.below(peers.len())];
            let q = peers[rng.below(peers.len())];
            state.apply(Move::Swap { p, q });
        }
        assert_eq!(state.translate(), protein, "translation must be untouched");
        assert_eq!(state.counts, counts_before, "a swap cannot change the multiset");
    }

    #[test]
    fn substitution_moves_exactly_one_count() {
        let (_, mut state) = gfp_state(2);
        let before = state.counts;
        let leu = aa_index(b'L').unwrap();
        let p = state.pos_by_aa[leu as usize][0];
        let old = state.codons[p as usize];
        let new = parse_codon("CTG").unwrap();
        let new = if new == old { parse_codon("CTT").unwrap() } else { new };
        state.apply(Move::Subst { p, codon: new });
        for (codon, &was) in before.iter().enumerate() {
            let expected =
                was - u32::from(codon == old as usize) + u32::from(codon == new as usize);
            assert_eq!(state.counts[codon], expected, "codon {codon}");
        }
    }

    #[test]
    fn bases_stay_in_step_with_codons() {
        let (_, mut state) = gfp_state(3);
        let mut rng = Rng::new(4);
        for _ in 0..2_000 {
            let aa = state.aa[rng.below(state.n_codons())];
            let peers = state.pos_by_aa[aa as usize].clone();
            if peers.len() < 2 {
                continue;
            }
            state.apply(Move::Swap {
                p: peers[rng.below(peers.len())],
                q: peers[rng.below(peers.len())],
            });
            let expected: Vec<Base> = state
                .codons
                .iter()
                .flat_map(|&c| codon_bases(c))
                .collect();
            assert_eq!(state.nt, expected);
        }
    }

    #[test]
    fn a_patched_view_reads_as_the_move_would_leave_it() {
        let (_, mut state) = gfp_state(5);
        let leu = aa_index(b'L').unwrap();
        let peers = state.pos_by_aa[leu as usize].clone();
        let (p, q) = (peers[0], peers[1]);
        if state.codons[p as usize] == state.codons[q as usize] {
            return; // nothing to observe
        }
        let mv = Move::Swap { p, q };
        let patch = state.patch_of(mv);
        let view = state.patched(patch);
        let (want_p, want_q) = (state.codons[q as usize], state.codons[p as usize]);
        assert_eq!(view.codon(p), want_p);
        assert_eq!(view.codon(q), want_q);
        for i in 0..3 {
            assert_eq!(view.base(p as usize * 3 + i), codon_bases(want_p)[i]);
        }
        // untouched positions read through unchanged
        let other = peers[2];
        assert_eq!(view.codon(other), state.codons[other as usize]);

        // and after applying, the plain view agrees with the patched one
        state.apply(mv);
        let plain = state.view();
        assert_eq!(plain.codon(p), want_p);
        assert_eq!(plain.codon(q), want_q);
    }

    #[test]
    fn affected_ranges_cover_the_margin_and_merge() {
        let (_, state) = gfp_state(6);
        let patch = Patch {
            first: (10, state.codons[10]),
            second: Some((11, state.codons[11])),
        };
        let view = state.scoped(patch);
        let ranges = view.affected(0);
        assert_eq!(ranges.as_slice(), &[(30, 36)], "adjacent codons merge into one range");

        let patch = Patch {
            first: (0, state.codons[0]),
            second: Some((100, state.codons[100])),
        };
        let ranges = state.scoped(patch).affected(9);
        assert_eq!(ranges.len(), 2, "distant codons stay separate");
        assert_eq!(ranges.as_slice()[0], (0, 12), "clamped at the start of the sequence");
        assert_eq!(ranges.as_slice()[1], (291, 312));
    }
}
