//! Constraints on the arrangement, and the machinery for evaluating a move against
//! them without mutating anything.
//!
//! Each constraint keeps a running violation total and knows how to answer two
//! questions cheaply: what would this move change (`delta`), and now that it has
//! happened, update yourself (`apply`). Both are given a *before* and an *after* view
//! of the same sequence, scoped to the positions the move touches, so a constraint only
//! ever examines a constant-sized neighbourhood.

use std::ops::Range;

use crate::arrange::state::{SeqState, View};

pub mod gc;
pub mod hairpin;
pub mod homopolymer;
pub mod repeat;
pub mod restriction;
pub mod splice;
pub mod start_site;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Kind {
    Gc,
    Homopolymer,
    Repeat,
    Hairpin,
    Restriction,
    StartSite,
    Splice,
}

impl Kind {
    pub fn label(self) -> &'static str {
        match self {
            Kind::Gc => "GC window",
            Kind::Homopolymer => "homopolymer",
            Kind::Repeat => "repeat",
            Kind::Hairpin => "hairpin",
            Kind::Restriction => "restriction site",
            Kind::StartSite => "start site",
            Kind::Splice => "splice motif",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Violation {
    pub kind: Kind,
    /// First base of the offending stretch.
    pub at: u32,
    /// Length in bases.
    pub len: u32,
    pub detail: String,
}

impl Violation {
    /// Codon positions the violation covers, which is what the search mutates.
    pub fn codons(&self) -> Range<u32> {
        let first = self.at / 3;
        let last = (self.at + self.len.max(1) - 1) / 3;
        first..last + 1
    }
}

pub trait Constraint: Send {
    fn kind(&self) -> Kind;
    /// Recount from scratch and reset any index.
    fn rebuild(&mut self, state: &SeqState);
    fn total(&self) -> u32;
    /// Change in the total this move would cause. Must not mutate.
    fn delta(&self, before: &View<'_>, after: &View<'_>) -> i32;
    /// Fold the move into the running total and any index. `before` and `after` are
    /// both taken against the pre-move state, as for `delta`.
    fn apply(&mut self, before: &View<'_>, after: &View<'_>);
    /// Every current violation, for the final report and for seeding the search queue.
    fn report(&self, state: &SeqState) -> Vec<Violation>;
}

/// Delta for a constraint whose violations are anchored within `margin` bases of any
/// change: count the affected neighbourhood before and after, and subtract.
pub fn window_delta<F>(before: &View<'_>, after: &View<'_>, margin: usize, count: F) -> i32
where
    F: Fn(&View<'_>, Range<usize>) -> u32,
{
    let mut delta = 0i32;
    for &(start, end) in after.affected(margin).iter() {
        delta += count(after, start..end) as i32;
        delta -= count(before, start..end) as i32;
    }
    delta
}

/// The whole constraint set, kept in one place so the search can score a move against
/// everything at once.
pub struct ConstraintSet {
    constraints: Vec<Box<dyn Constraint>>,
}

impl ConstraintSet {
    pub fn new(constraints: Vec<Box<dyn Constraint>>) -> Self {
        Self { constraints }
    }

    pub fn rebuild(&mut self, state: &SeqState) {
        for c in self.constraints.iter_mut() {
            c.rebuild(state);
        }
    }

    pub fn total(&self) -> u32 {
        self.constraints.iter().map(|c| c.total()).sum()
    }

    pub fn delta(&self, before: &View<'_>, after: &View<'_>) -> i32 {
        self.constraints.iter().map(|c| c.delta(before, after)).sum()
    }

    pub fn apply(&mut self, before: &View<'_>, after: &View<'_>) {
        for c in self.constraints.iter_mut() {
            c.apply(before, after);
        }
    }

    pub fn report(&self, state: &SeqState) -> Vec<Violation> {
        let mut all: Vec<Violation> = self
            .constraints
            .iter()
            .flat_map(|c| c.report(state))
            .collect();
        all.sort_by_key(|v| (v.at, v.kind));
        all
    }

    pub fn totals_by_kind(&self) -> Vec<(Kind, u32)> {
        self.constraints
            .iter()
            .map(|c| (c.kind(), c.total()))
            .filter(|(_, n)| *n > 0)
            .collect()
    }
}

#[cfg(test)]
pub mod harness {
    //! Shared differential tests. Every constraint is checked two ways: its count must
    //! match a brute-force oracle, and its incremental bookkeeping must match a full
    //! rebuild after every single move.

    use super::*;
    use crate::arrange::state::Move;
    use crate::code::{codon_aa, AaIdx, CodonIdx, N_AA};
    use crate::rng::Rng;

    /// A random in-frame sequence over the standard genetic code.
    pub fn random_state(rng: &mut Rng, n_codons: usize) -> SeqState {
        let codons: Vec<CodonIdx> = (0..n_codons).map(|_| rng.below(64) as CodonIdx).collect();
        let aa: Vec<AaIdx> = codons.iter().copied().map(codon_aa).collect();
        SeqState::new(aa, codons)
    }

    /// Random moves that keep the protein: swap two positions with the same residue, or
    /// substitute a synonymous codon.
    pub fn random_move(rng: &mut Rng, state: &SeqState) -> Option<Move> {
        let p = rng.below(state.n_codons()) as u32;
        let aa = state.aa[p as usize];
        let peers = &state.pos_by_aa[aa as usize];
        if peers.len() > 1 && rng.chance(0.7) {
            let q = peers[rng.below(peers.len())];
            if q != p {
                return Some(Move::Swap { p, q });
            }
        }
        let synonyms: Vec<CodonIdx> = (0..64u8).filter(|&c| codon_aa(c) == aa).collect();
        if synonyms.len() < 2 {
            return None;
        }
        let codon = synonyms[rng.below(synonyms.len())];
        (codon != state.codons[p as usize]).then_some(Move::Subst { p, codon })
    }

    /// The count must equal the brute-force oracle on fresh random sequences.
    pub fn assert_matches_oracle<C, F>(mut make: impl FnMut() -> C, oracle: F, sequences: usize)
    where
        C: Constraint,
        F: Fn(&SeqState) -> u32,
    {
        let mut rng = Rng::new(0x5EED);
        for round in 0..sequences {
            let n = rng.between(4, 120);
            let state = random_state(&mut rng, n);
            let mut c = make();
            c.rebuild(&state);
            assert_eq!(
                c.total(),
                oracle(&state),
                "round {round}: {:?} disagrees with the oracle on {}",
                c.kind(),
                state.dna()
            );
            assert_eq!(
                c.report(&state).len() as u32,
                c.total(),
                "round {round}: {:?} reports a different number than it counts",
                c.kind()
            );
        }
    }

    /// After every move, the incremental total must equal a full rebuild, and the delta
    /// predicted beforehand must match what actually happened.
    pub fn assert_incremental_matches_rebuild<C>(mut make: impl FnMut() -> C, moves: usize)
    where
        C: Constraint,
    {
        let mut rng = Rng::new(0xC0FFEE);
        let mut state = random_state(&mut rng, 90);
        let mut c = make();
        c.rebuild(&state);
        let mut fresh = make();

        let mut applied = 0;
        while applied < moves {
            let Some(mv) = random_move(&mut rng, &state) else {
                continue;
            };
            let patch = state.patch_of(mv);
            let before = state.scoped(patch);
            let after = state.patched(patch);

            let predicted = c.delta(&before, &after);
            let was = c.total();
            c.apply(&before, &after);
            state.apply(mv);

            fresh.rebuild(&state);
            assert_eq!(
                c.total(),
                fresh.total(),
                "{:?}: incremental total drifted after move {applied} ({mv:?}) on {}",
                c.kind(),
                state.dna()
            );
            assert_eq!(
                was as i32 + predicted,
                c.total() as i32,
                "{:?}: predicted delta {predicted} was wrong at move {applied} ({mv:?})",
                c.kind()
            );
            applied += 1;

            // occasionally start over, so the test sees many sequence shapes
            if applied % 500 == 0 {
                let n = rng.between(20, 140);
                state = random_state(&mut rng, n);
                c = make();
                c.rebuild(&state);
            }
        }
        let _ = N_AA;
    }
}
