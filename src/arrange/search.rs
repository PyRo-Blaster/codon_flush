//! Min-conflicts search over arrangements.
//!
//! The composition is already fixed, so the default move is a swap between two
//! positions holding the same amino acid: it cannot change the codon multiset, and
//! therefore cannot change the deviation from the host profile or the CAI. Only when no
//! swap can clear a violation does the search spend budget on a substitution, which
//! does change the multiset.

use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::arrange::constraints::{
    gc, hairpin, homopolymer, repeat, restriction, splice, start_site, Constraint,
    ConstraintSet, Kind, Violation,
};
use crate::arrange::state::{Move, SeqState};
use crate::code::{codon_aa, codon_bases, AaIdx, Base, CodonIdx, N_AA, N_CODONS};
use crate::composition::Composition;
use crate::host::HostTable;
use crate::metrics;
use crate::rng::Rng;

#[derive(Debug, Clone)]
pub struct Options {
    pub seed: u64,
    pub restarts: usize,
    /// Cap on accepted moves, as a multiple of the gene's codon count.
    pub moves_per_codon: usize,
    /// Give up on a restart after this many moves without an improvement.
    pub stall_limit: usize,
    /// Chance of taking a move that leaves the count unchanged, to cross plateaus.
    pub sideways: f64,
    /// Chance of taking a random move instead of the best one, to escape local minima.
    pub noise: f64,
    /// Swap partners sampled per position.
    pub peers: usize,
    /// Upper bound on the moves examined by the fallback sweep for one violation.
    pub exhaustive_cap: usize,
    pub homopolymer_limit: usize,
    pub enzymes: Vec<(String, String)>,
    pub remove_start_sites: bool,
    pub start_site_rule: start_site::Rule,
    pub remove_splice_sites: bool,
    /// What counts as a splice problem: no motif at all, or no donor-acceptor pair.
    pub splice_rule: splice::Rule,
    /// Whole-gene GC bounds. Swaps cannot move whole-gene GC, but substitutions
    /// change the multiset and so can walk it out of the band; this keeps them
    /// inside it rather than leaving acceptance to throw the arrangement away.
    pub gc_band: (f64, f64),
}

impl Default for Options {
    fn default() -> Self {
        Self {
            seed: 1,
            restarts: 4,
            moves_per_codon: 200,
            stall_limit: 2000,
            sideways: 0.3,
            noise: 0.05,
            peers: 8,
            exhaustive_cap: 128,
            homopolymer_limit: 4,
            enzymes: restriction::DEFAULT_ENZYMES
                .iter()
                .map(|(n, s)| (n.to_string(), s.to_string()))
                .collect(),
            remove_start_sites: true,
            start_site_rule: start_site::Rule::ShineDalgarno,
            remove_splice_sites: true,
            splice_rule: splice::Rule::NoMotif,
            gc_band: (0.30, 0.65),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Outcome {
    pub state: SeqState,
    pub violations: u32,
    pub by_kind: Vec<(Kind, u32)>,
    pub remaining: Vec<Violation>,
    pub deviation: f64,
    pub cai: f64,
    pub gc: f64,
    pub moves: usize,
    pub substitutions: usize,
    pub restart: usize,
    /// Splice motifs the finished gene still contains, whether or not the rule
    /// objects to them. Under the pair rule a gene can be valid and still carry a
    /// lone donor, and the caller has to be able to say where it is.
    pub retained_splice_sites: Vec<splice::Site>,
}

fn build_constraints(options: &Options) -> ConstraintSet {
    let mut constraints: Vec<Box<dyn Constraint>> = vec![
        Box::new(gc::GcWindows::new(gc::DEFAULT_WINDOWS.to_vec())),
        Box::new(homopolymer::Homopolymers::new(options.homopolymer_limit)),
        Box::new(repeat::Repeats::new()),
        Box::new(hairpin::Hairpins::new()),
    ];
    if !options.enzymes.is_empty() {
        constraints.push(Box::new(restriction::RestrictionSites::new(&options.enzymes)));
    }
    if options.remove_start_sites {
        constraints.push(Box::new(start_site::StartSites::with_rule(options.start_site_rule)));
    }
    if options.remove_splice_sites {
        constraints.push(Box::new(splice::SpliceSites::with_rule(options.splice_rule)));
    }
    ConstraintSet::new(constraints)
}

/// Lay the chosen multiset out over the protein, left to right, taking at each
/// position whichever remaining codon introduces the fewest problems in the bases just
/// written. Starting from a sensible arrangement leaves the search far less to undo
/// than starting from a shuffle.
fn seed_arrangement(
    protein: &[AaIdx],
    composition: &Composition,
    options: &Options,
    rng: &mut Rng,
) -> SeqState {
    // Counts rather than a pool of individual codons: the choices at any position are
    // the handful of synonyms still available, found in constant time. Cloning and
    // sorting a pool that grows with the gene made seeding quadratic.
    let mut remaining = composition.counts;
    let mut synonyms: Vec<Vec<CodonIdx>> = vec![Vec::new(); N_AA];
    for codon in 0..N_CODONS as CodonIdx {
        synonyms[codon_aa(codon) as usize].push(codon);
    }
    for list in synonyms.iter_mut() {
        rng.shuffle(list);
    }

    let sites: Vec<Vec<Base>> = options
        .enzymes
        .iter()
        .map(|(_, site)| site.bytes().filter_map(crate::code::base_index).collect())
        .collect();
    let motifs = splice::motif_masks();

    let mut nt: Vec<Base> = Vec::with_capacity(protein.len() * 3);
    let mut codons: Vec<CodonIdx> = Vec::with_capacity(protein.len());
    let mut seen: HashSet<u32> = HashSet::new();

    for &aa in protein {
        let mut best: Option<(u32, CodonIdx)> = None;
        for &candidate in &synonyms[aa as usize] {
            if remaining[candidate as usize] == 0 {
                continue;
            }
            let penalty = local_penalty(
                &nt, &codons, &seen, candidate, &sites, &motifs, options,
            );
            if best.map_or(true, |(score, _)| penalty < score) {
                best = Some((penalty, candidate));
            }
            if penalty == 0 {
                break; // nothing to gain from looking further
            }
        }
        let chosen = best.map(|(_, c)| c).expect("the composition supplies every residue");
        remaining[chosen as usize] -= 1;

        codons.push(chosen);
        nt.extend_from_slice(&codon_bases(chosen));
        if nt.len() >= 9 {
            let mut key = 0u32;
            for &base in &nt[nt.len() - 9..] {
                key = (key << 2) | base as u32;
            }
            seen.insert(key);
        }
    }
    SeqState::new(protein.to_vec(), codons)
}

/// Problems a candidate codon would create in the bases it adds. Covers the local
/// rules; whole-window GC and hairpins are left to the search, which sees the finished
/// sequence.
fn local_penalty(
    nt: &[Base],
    codons: &[CodonIdx],
    seen: &HashSet<u32>,
    candidate: CodonIdx,
    sites: &[Vec<Base>],
    motifs: &[Vec<u8>],
    options: &Options,
) -> u32 {
    let mut buffer = [0u8; 24];
    let keep = nt.len().min(21);
    buffer[..keep].copy_from_slice(&nt[nt.len() - keep..]);
    buffer[keep..keep + 3].copy_from_slice(&codon_bases(candidate));
    let tail = &buffer[..keep + 3];
    let added_at = keep;
    let mut penalty = 0u32;

    // a run of identical bases reaching past the limit
    let mut run = 1usize;
    let mut index = tail.len() - 1;
    while index > 0 && tail[index - 1] == tail[index] {
        run += 1;
        index -= 1;
    }
    if run > options.homopolymer_limit {
        penalty += 1;
    }

    // three identical codons in a row
    if codons.len() >= 2
        && codons[codons.len() - 1] == candidate
        && codons[codons.len() - 2] == candidate
    {
        penalty += 1;
    }

    // a restriction site or splice motif ending in the new bases
    for site in sites {
        for start in tail.len().saturating_sub(site.len() + 2)..=tail.len().saturating_sub(site.len()) {
            if start + site.len() <= tail.len()
                && start + site.len() > added_at
                && (0..site.len()).all(|k| tail[start + k] == site[k])
            {
                penalty += 1;
            }
        }
    }
    if options.remove_splice_sites {
        for masks in motifs {
            for start in tail.len().saturating_sub(masks.len() + 2)..=tail.len().saturating_sub(masks.len()) {
                if start + masks.len() <= tail.len()
                    && start + masks.len() > added_at
                    && masks.iter().enumerate().all(|(k, &m)| m & (1 << tail[start + k]) != 0)
                {
                    penalty += 1;
                }
            }
        }
    }

    // a 9-mer that has already been used
    if tail.len() >= 9 && nt.len() + 3 >= 9 {
        for end in (added_at + 1).max(9)..=tail.len() {
            let mut key = 0u32;
            for &base in &tail[end - 9..end] {
                key = (key << 2) | base as u32;
            }
            if seen.contains(&key) {
                penalty += 1;
            }
        }
    }
    penalty
}

struct Candidate {
    mv: Move,
    delta: i32,
    deviation_cost: f64,
}

/// Run one restart to completion.
fn run_once(
    protein: &[AaIdx],
    composition: &Composition,
    host: &HostTable,
    options: &Options,
    restart: usize,
) -> Outcome {
    let mut rng = Rng::new(options.seed ^ (restart as u64).wrapping_mul(0x9E37_79B9));
    let mut state = seed_arrangement(protein, composition, options, &mut rng);
    let mut constraints = build_constraints(options);
    constraints.rebuild(&state);

    let max_moves = options.moves_per_codon * state.n_codons().max(1);
    let mut queue: Vec<Violation> = Vec::new();
    // Violations that no single move can improve. Some are genuinely unsatisfiable:
    // Glu is GAA/GAG and Lys is AAA/AAG, so "EK" is six purines whatever you choose,
    // and if an in-frame ATG follows, no codon assignment removes that start site.
    // Violation span -> how many times a full sweep of it has come up empty.
    let mut blocked: HashMap<(u32, u32), u32> = HashMap::new();
    let mut moves = 0usize;
    let mut substitutions = 0usize;
    // A whole-run allowance for the expensive fallback sweep, so its cost is linear in
    // the gene's length however awkward the protein turns out to be.
    let mut sweep_budget = 40 * state.n_codons();
    let mut best_seen = constraints.total();
    let mut since_improvement = 0usize;

    while constraints.total() > 0 && moves < max_moves && since_improvement < options.stall_limit {
        if queue.is_empty() {
            queue = constraints
                .report(&state)
                .into_iter()
                .filter(|v| !blocked.contains_key(&(v.at, v.len)))
                .collect();
            rng.shuffle(&mut queue);
            if queue.is_empty() {
                break; // everything left is known to be unfixable
            }
        }
        let violation = queue.pop().expect("queue is not empty");
        let span = violation.codons();
        let positions: Vec<u32> = (span.start..span.end.min(state.n_codons() as u32)).collect();
        if positions.is_empty() {
            continue;
        }

        // Min-conflicts: one position from the offending window, a sample of partners.
        let p = positions[rng.below(positions.len())];
        let mut best = sample_swaps(&state, &constraints, &mut rng, p, options.peers);

        let mut chosen = match &best {
            Some(c) if c.delta < 0 => Some(c.mv),
            Some(c) if c.delta == 0 && rng.chance(options.sideways) => Some(c.mv),
            _ if rng.chance(options.noise) => best.as_ref().map(|c| c.mv),
            _ => None,
        };

        if chosen.is_none() {
            // Nothing in the sample helped. Look at every partner for every position
            // before concluding this violation is stuck.
            let cap = options.exhaustive_cap.min(sweep_budget);
            if cap > 0 {
                let (candidate, examined) =
                    exhaustive_swaps(&state, &constraints, &positions, cap);
                sweep_budget -= examined;
                best = candidate;
            }
            chosen = match &best {
                Some(c) if c.delta < 0 => Some(c.mv),
                _ => None,
            };
            if chosen.is_none() {
                if let Some(candidate) =
                    best_substitution(
                        &state, &constraints, host, &positions, composition,
                        options.gc_band, 0,
                    )
                {
                    substitutions += 1;
                    chosen = Some(candidate.mv);
                } else if sweep_budget > 0
                    && blocked.get(&(violation.at, violation.len)).copied().unwrap_or(0) == 0
                    && try_paired_substitution(
                    &mut state,
                    &mut constraints,
                    host,
                    &positions,
                    composition,
                    options.gc_band,
                    options.exhaustive_cap,
                ) {
                    // two coordinated substitutions were needed, and they paid off
                    substitutions += 2;
                    moves += 2;
                } else {
                    *blocked.entry((violation.at, violation.len)).or_insert(0) += 1;
                }
            }
        }

        if let Some(mv) = chosen {
            let touched = touched_bases(&state, mv);
            apply(&mut state, &mut constraints, mv);
            moves += 1;
            // Only violations near the change get another look. Clearing the whole set
            // after every improvement made the search re-sweep every stuck violation in
            // the gene, which is quadratic in the gene's length.
            unblock_near(&mut blocked, touched);
        }

        let total = constraints.total();
        if total < best_seen {
            best_seen = total;
            since_improvement = 0;
        } else {
            since_improvement += 1;
        }
    }

    let remaining = constraints.report(&state);
    // Listed from a fresh detector rather than the set, which hands out only the
    // violations the rule objects to. A retained lone site is not a violation and
    // still has to be disclosed.
    let retained_splice_sites = if options.remove_splice_sites {
        splice::SpliceSites::with_rule(options.splice_rule).sites(&state)
    } else {
        Vec::new()
    };
    Outcome {
        violations: constraints.total(),
        by_kind: constraints.totals_by_kind(),
        remaining,
        retained_splice_sites,
        deviation: metrics::deviation_of(&state, host),
        cai: metrics::cai(&state, host),
        gc: metrics::gc_of(&state),
        moves,
        substitutions,
        restart,
        state,
    }
}

/// A sample of swap partners for one position.
fn sample_swaps(
    state: &SeqState,
    constraints: &ConstraintSet,
    rng: &mut Rng,
    p: u32,
    peers_to_try: usize,
) -> Option<Candidate> {
    let peers = &state.pos_by_aa[state.aa[p as usize] as usize];
    if peers.len() < 2 {
        return None;
    }
    let mut best: Option<Candidate> = None;
    for _ in 0..peers_to_try.min(peers.len()) {
        let q = peers[rng.below(peers.len())];
        if q == p || state.codons[q as usize] == state.codons[p as usize] {
            continue;
        }
        let mv = Move::Swap { p, q };
        let patch = state.patch_of(mv);
        let delta = constraints.delta(&state.scoped(patch), &state.patched(patch));
        if best.as_ref().map_or(true, |b| delta < b.delta) {
            best = Some(Candidate { mv, delta, deviation_cost: 0.0 });
        }
    }
    best
}

/// Every partner for every position in the window, used only when the sample fails.
fn exhaustive_swaps(
    state: &SeqState,
    constraints: &ConstraintSet,
    positions: &[u32],
    cap: usize,
) -> (Option<Candidate>, usize) {
    let mut best: Option<Candidate> = None;
    let mut examined = 0usize;
    for &p in positions {
        for &q in &state.pos_by_aa[state.aa[p as usize] as usize] {
            if examined >= cap {
                return (best, examined);
            }
            examined += 1;
            if q == p || state.codons[q as usize] == state.codons[p as usize] {
                continue;
            }
            let mv = Move::Swap { p, q };
            let patch = state.patch_of(mv);
            let delta = constraints.delta(&state.scoped(patch), &state.patched(patch));
            if best.as_ref().map_or(true, |b| delta < b.delta) {
                best = Some(Candidate { mv, delta, deviation_cost: 0.0 });
            }
            if delta < 0 {
                return (best, examined);
            }
        }
    }
    (best, examined)
}

/// Some violations need two substitutions before anything improves: four lysines with
/// three AAA in the multiset always give six A's in a row, and only swapping two of
/// them for AAG fixes it. Try one substitution that holds the line, look for a move
/// that now improves, and undo if none appears, so budget is only ever spent on a fix.
fn try_paired_substitution(
    state: &mut SeqState,
    constraints: &mut ConstraintSet,
    host: &HostTable,
    positions: &[u32],
    composition: &Composition,
    gc_band: (f64, f64),
    cap: usize,
) -> bool {
    let Some(holding) =
        best_substitution(state, constraints, host, positions, composition, gc_band, 1)
    else {
        return false;
    };
    let Move::Subst { p, .. } = holding.mv else {
        return false;
    };
    let original = state.codons[p as usize];
    apply(state, constraints, holding.mv);

    let follow = exhaustive_swaps(state, constraints, positions, cap)
        .0
        .filter(|c| c.delta < 0)
        .or_else(|| {
            best_substitution(state, constraints, host, positions, composition, gc_band, 0)
        });
    match follow {
        Some(candidate) => {
            apply(state, constraints, candidate.mv);
            true
        }
        None => {
            apply(state, constraints, Move::Subst { p, codon: original });
            false
        }
    }
}

/// Bases a move rewrites.
fn touched_bases(state: &SeqState, mv: Move) -> (u32, u32) {
    let patch = state.patch_of(mv);
    let positions: Vec<u32> = patch.positions().collect();
    let lo = positions.iter().copied().min().unwrap_or(0) * 3;
    let hi = positions.iter().copied().max().unwrap_or(0) * 3 + 3;
    (lo, hi)
}

/// After this many empty sweeps, a violation is left alone for the rest of the run.
/// Without it, a violation nothing can fix is swept again every time a move lands near
/// it, which is what made long genes slow.
const GIVE_UP_AFTER: u32 = 2;

/// Reconsider a stuck violation only when a move actually rewrote part of it.
fn unblock_near(blocked: &mut HashMap<(u32, u32), u32>, touched: (u32, u32)) {
    blocked.retain(|&(at, len), failures| {
        let overlaps = at < touched.1 && touched.0 < at + len;
        !overlaps || *failures >= GIVE_UP_AFTER
    });
}

fn apply(state: &mut SeqState, constraints: &mut ConstraintSet, mv: Move) {
    let patch = state.patch_of(mv);
    let before = state.scoped(patch);
    let after = state.patched(patch);
    constraints.apply(&before, &after);
    state.apply(mv);
}

/// The cheapest substitution that clears violations without overrunning the budget.
#[allow(clippy::too_many_arguments)]
fn best_substitution(
    state: &SeqState,
    constraints: &ConstraintSet,
    host: &HostTable,
    positions: &[u32],
    composition: &Composition,
    gc_band: (f64, f64),
    allow_above: i32,
) -> Option<Candidate> {
    if composition.budget_left <= 0.0 {
        return None; // no budget means the multiset is frozen
    }
    let budget = composition.deviation + composition.budget_left;
    let mut best: Option<Candidate> = None;
    for &p in positions {
        let aa = state.aa[p as usize];
        let current = state.codons[p as usize];
        for codon in host.usable_synonyms(aa) {
            if codon == current {
                continue;
            }
            let mut counts = state.counts;
            counts[current as usize] -= 1;
            counts[codon as usize] += 1;
            let deviation = crate::composition::deviation(&counts, host);
            if deviation > budget {
                continue;
            }
            // The band is the caller's requirement, not a preference, so a
            // substitution that leaves it is not available however much it would
            // buy elsewhere.
            let gc = crate::composition::gc_fraction(&counts);
            if gc < gc_band.0 || gc > gc_band.1 {
                continue;
            }
            let mv = Move::Subst { p, codon };
            let patch = state.patch_of(mv);
            let delta = constraints.delta(&state.scoped(patch), &state.patched(patch));
            if delta >= allow_above {
                continue;
            }
            let better = best.as_ref().map_or(true, |b| {
                (delta, deviation) < (b.delta, b.deviation_cost)
            });
            if better {
                best = Some(Candidate { mv, delta, deviation_cost: deviation });
            }
        }
    }
    best
}

/// Optimise one protein. Restarts run in parallel and the best is chosen
/// deterministically, so the result does not depend on thread scheduling.
/// Every restart's arrangement, best first.
///
/// Restarts of one composition share its multiset and so share its deviation and CAI,
/// but they are different sequences, which is what a candidate set needs. Reducing to
/// one outcome throws that away, so the reduction lives in `optimise` instead.
pub fn optimise_all(
    protein: &[AaIdx],
    composition: &Composition,
    host: &HostTable,
    options: &Options,
) -> Vec<Outcome> {
    let restarts = options.restarts.max(1);
    let mut outcomes: Vec<Outcome> = (0..restarts)
        .into_par_iter()
        .map(|restart| run_once(protein, composition, host, options, restart))
        .collect();
    // Sorted by index-derived fields only, so the order never depends on which thread
    // finished first.
    outcomes.sort_by(|a, b| {
        a.violations
            .cmp(&b.violations)
            .then(a.deviation.partial_cmp(&b.deviation).unwrap())
            .then(b.cai.partial_cmp(&a.cai).unwrap())
            .then(a.restart.cmp(&b.restart))
    });
    outcomes
}

/// The single best arrangement. The original one-outcome API, kept for callers that
/// want one gene rather than a candidate set.
pub fn optimise(
    protein: &[AaIdx],
    composition: &Composition,
    host: &HostTable,
    options: &Options,
) -> Outcome {
    let mut outcomes = optimise_all(protein, composition, host, options);
    outcomes.swap_remove(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::aa_index;
    use crate::composition;
    use crate::testdata::{ECOLI_TABLE, GFP, STRESS};
    use std::path::Path;

    fn host() -> HostTable {
        HostTable::load(Path::new(ECOLI_TABLE), 0.10).unwrap()
    }

    fn protein(seq: &str) -> Vec<AaIdx> {
        seq.bytes().map(|b| aa_index(b).unwrap()).collect()
    }

    fn run(seq: &str, options: Options) -> (Outcome, Composition) {
        let host = host();
        let protein = protein(seq);
        let comp = composition::build(&protein, &host, 0.1, (0.30, 0.65)).unwrap();
        (optimise(&protein, &comp, &host, &options), comp)
    }

    #[test]
    fn clears_gfp_completely() {
        let (outcome, comp) = run(GFP, Options { seed: 1, ..Default::default() });
        assert_eq!(outcome.state.translate(), protein(GFP), "protein must survive");
        assert_eq!(outcome.violations, 0, "left {:?}", outcome.remaining);
        assert_eq!(outcome.deviation, comp.deviation, "no budget needed to be spent");
        assert!(outcome.deviation < 0.01, "deviation {}", outcome.deviation);
        assert!(outcome.cai > 0.70, "cai {}", outcome.cai);
    }

    #[test]
    fn the_old_purine_rule_flags_sites_no_move_can_fix() {
        // Ala-Glu-Lys-Ala-Ala-Ala-Met puts GA[AG]AA[AG] in the window before an
        // in-frame ATG. Glu and Lys are purine-only in every encoding and Met has no
        // alternative codon, so nothing can clear it. This is why the default rule now
        // asks whether the upstream sequence could pair with the 16S rRNA instead.
        let host = host();
        let protein = protein("AEKAAAM");
        let comp = composition::build(&protein, &host, 0.1, (0.0, 1.0)).unwrap();
        let legacy = Options {
            seed: 5,
            restarts: 2,
            start_site_rule: start_site::Rule::PurineRun,
            ..Default::default()
        };
        let outcome = optimise(&protein, &comp, &host, &legacy);
        let left = outcome
            .remaining
            .iter()
            .find(|v| v.kind == Kind::StartSite)
            .expect("the old rule should still report it");
        assert!(left.detail.contains("ATG"), "{}", left.detail);

        // and it really is unsatisfiable under that rule: no single move improves it
        let mut set = ConstraintSet::new(vec![Box::new(start_site::StartSites::with_rule(
            start_site::Rule::PurineRun,
        ))]);
        set.rebuild(&outcome.state);
        let n = outcome.state.n_codons() as u32;
        for p in 0..n {
            let aa = outcome.state.aa[p as usize];
            for q in 0..n {
                if q == p || outcome.state.aa[q as usize] != aa {
                    continue;
                }
                let patch = outcome.state.patch_of(Move::Swap { p, q });
                let delta = set.delta(&outcome.state.scoped(patch), &outcome.state.patched(patch));
                assert!(delta >= 0, "a swap at {p}/{q} would have fixed it");
            }
            for codon in host.usable_synonyms(aa) {
                let patch = outcome.state.patch_of(Move::Subst { p, codon });
                let delta = set.delta(&outcome.state.scoped(patch), &outcome.state.patched(patch));
                assert!(delta >= 0, "a substitution at {p} would have fixed it");
            }
        }

        // the default rule does not flag it: GAAAAG has no adjacent G's, so it cannot
        // recruit a ribosome however many purines it has
        let outcome = optimise(&protein, &comp, &host, &Options { seed: 5, restarts: 2, ..Default::default() });
        assert!(outcome.remaining.iter().all(|v| v.kind != Kind::StartSite));
    }

    #[test]
    fn swaps_alone_never_move_deviation_or_cai() {
        let host = host();
        let protein = protein(GFP);
        let comp = composition::build(&protein, &host, 0.0, (0.30, 0.65)).unwrap();
        // max_relax 0 leaves no budget, so substitutions cannot be used at all
        let mut seen = Vec::new();
        for seed in 1..=5 {
            let outcome = optimise(
                &protein,
                &comp,
                &host,
                &Options { seed, restarts: 2, ..Default::default() },
            );
            assert_eq!(outcome.substitutions, 0, "no budget means no substitutions");
            assert_eq!(
                outcome.deviation, comp.deviation,
                "a swap cannot change the deviation"
            );
            seen.push(outcome.cai);
        }
        for pair in seen.windows(2) {
            assert_eq!(pair[0], pair[1], "CAI is a property of the multiset alone");
        }
    }

    #[test]
    fn a_seed_reproduces_its_run() {
        let a = run(GFP, Options { seed: 7, ..Default::default() }).0;
        let b = run(GFP, Options { seed: 7, ..Default::default() }).0;
        assert_eq!(a.state.dna(), b.state.dna());
        assert_eq!(a.moves, b.moves);

        let single = run(GFP, Options { seed: 7, restarts: 1, ..Default::default() }).0;
        let many = run(GFP, Options { seed: 7, restarts: 8, ..Default::default() }).0;
        // more restarts may find a better arrangement, but each restart is reproducible
        assert!(many.violations <= single.violations);
    }

    #[test]
    fn handles_a_protein_with_no_choices() {
        // Trp and Met have one codon each, so no move exists at all
        let (outcome, _) = run("MWMWMWMWMWMW", Options { seed: 2, ..Default::default() });
        assert_eq!(outcome.state.translate(), protein("MWMWMWMWMWMW"));
        assert_eq!(outcome.moves, 0, "there is nothing to try");
    }

    #[test]
    fn reports_what_it_cannot_fix() {
        let host = host();
        let protein = protein("WWWWWWWWWWWWWWWWWWWW");
        // TGG pins GC at two thirds, so this one needs a band it can actually meet
        let comp = composition::build(&protein, &host, 0.1, (0.0, 1.0)).unwrap();
        let outcome = optimise(
            &protein,
            &comp,
            &host,
            &Options { seed: 3, restarts: 1, ..Default::default() },
        );
        // poly-Trp is TGGTGG..., which repeats and cannot be rearranged
        assert!(outcome.violations > 0);
        assert_eq!(outcome.remaining.len(), outcome.violations as usize);
        assert!(outcome.remaining.iter().all(|v| !v.detail.is_empty()));
    }

    #[test]
    fn clears_the_stress_protein_or_says_why() {
        let (outcome, comp) = run(STRESS, Options { seed: 4, ..Default::default() });
        assert_eq!(outcome.state.translate(), protein(STRESS));
        assert!(outcome.deviation <= comp.deviation + comp.budget_left + 1e-9);
        if outcome.violations > 0 {
            assert_eq!(outcome.remaining.len(), outcome.violations as usize);
        }
    }
}
