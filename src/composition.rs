//! Composition: how many of each codon the gene uses.
//!
//! This half of the problem decides the deviation from the host profile, the codon
//! adaptation index and the whole-gene GC content, and nothing else. It is solved once,
//! exactly, before any arrangement work begins.
//!
//! Largest-remainder (Hamilton) apportionment minimises `sum |k_c - f_c * n_a|` subject
//! to the counts summing to the number of positions, so it is the L1-optimal integer
//! allocation of a synonymous group.

use std::collections::HashSet;

use crate::code::{codon_aa, codon_gc, codon_string, AaIdx, CodonIdx, N_AA, N_CODONS};
use crate::host::HostTable;

#[derive(Debug, Clone)]
pub struct Composition {
    /// How many times each codon is used.
    pub counts: [u32; N_CODONS],
    /// Deviation from the host profile, using the metric the Python reports.
    pub deviation: f64,
    /// Whole-gene GC fraction implied by these counts.
    pub gc: f64,
    /// Codon adaptation index implied by these counts. Arrangement cannot change it.
    pub cai: f64,
    /// Deviation budget still available to the arrangement search.
    pub budget_left: f64,
    /// Codon moves spent pulling GC into its band.
    pub gc_moves: u32,
    /// Whether whole-gene GC ended up inside its band. It can be false when correcting
    /// it would cost more deviation than `max_relax` allows; the gene is still built,
    /// and the shortfall is reported rather than paid for out of budget the caller
    /// said was not available.
    pub gc_in_band: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CompositionError {
    /// GC cannot reach the band by any synonymous reallocation.
    GcUnreachable { best: f64, band: (f64, f64) },
    /// Even the L1-optimal allocation deviates more than a positive limit allows.
    /// Integer rounding alone costs this much, so no multiset can do better.
    DeviationFloor { best: f64, limit: f64 },
    /// GC can reach the band, but not without deviating more than the limit allows.
    GcOutsideBudget { best: f64, band: (f64, f64), limit: f64 },
    /// The protein has no residues.
    Empty,
}

impl std::fmt::Display for CompositionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::GcUnreachable { best, band } => write!(
                f,
                "whole-gene GC cannot reach {:.0}-{:.0}% for this protein; the closest \
                 any synonymous choice gets is {:.1}%",
                band.0 * 100.0,
                band.1 * 100.0,
                best * 100.0
            ),
            Self::DeviationFloor { best, limit } => write!(
                f,
                "the smallest deviation any codon multiset can have for this protein is \
                 {best:.6}, above the requested limit of {limit:.6}; integer rounding of \
                 the host profile alone costs that much, so raise --max-relax to at \
                 least {best:.6}"
            ),
            Self::GcOutsideBudget { best, band, limit } => write!(
                f,
                "whole-gene GC cannot reach {:.0}-{:.0}% without deviating from the host \
                 profile by more than the requested {limit:.6}; the closest any \
                 composition within that limit gets is {:.1}%",
                band.0 * 100.0,
                band.1 * 100.0,
                best * 100.0
            ),
            Self::Empty => write!(f, "the protein is empty"),
        }
    }
}

impl std::error::Error for CompositionError {}

/// Python's `round()`: halves go to the even neighbour. The deviation metric is defined
/// in terms of it, so matching it keeps the numbers comparable with the old tool.
fn round_half_even(x: f64) -> f64 {
    let floor = x.floor();
    let frac = x - floor;
    if frac < 0.5 {
        floor
    } else if frac > 0.5 {
        floor + 1.0
    } else if (floor as i64) % 2 == 0 {
        floor
    } else {
        floor + 1.0
    }
}

/// Deviation from the host profile, as `compare_profiles` in the Python defines it:
/// half the total absolute difference from the per-codon ideal, over the total ideal.
pub fn deviation(counts: &[u32; N_CODONS], host: &HostTable) -> f64 {
    let mut group_total = [0u32; N_AA];
    for codon in 0..N_CODONS {
        group_total[codon_aa(codon as CodonIdx) as usize] += counts[codon];
    }
    let mut diff_total = 0.0;
    let mut residues = 0.0;
    for codon in 0..N_CODONS as CodonIdx {
        let n_a = group_total[codon_aa(codon) as usize] as f64;
        let ideal = round_half_even(host.freq(codon) * n_a);
        diff_total += (counts[codon as usize] as f64 - ideal).abs();
        residues += ideal;
    }
    if residues == 0.0 {
        return 0.0;
    }
    (diff_total / 2.0) / residues
}

pub fn gc_fraction(counts: &[u32; N_CODONS]) -> f64 {
    let mut gc = 0u32;
    let mut total = 0u32;
    for codon in 0..N_CODONS as CodonIdx {
        gc += codon_gc(codon) * counts[codon as usize];
        total += 3 * counts[codon as usize];
    }
    if total == 0 {
        0.0
    } else {
        gc as f64 / total as f64
    }
}

/// Largest-remainder apportionment of one synonymous group.
fn apportion(host: &HostTable, aa: AaIdx, n: u32, counts: &mut [u32; N_CODONS]) {
    let usable: Vec<CodonIdx> = host.usable_synonyms(aa).collect();
    if usable.is_empty() || n == 0 {
        return;
    }
    let mut assigned = 0u32;
    let mut remainders: Vec<(f64, f64, CodonIdx)> = Vec::with_capacity(usable.len());
    for &codon in &usable {
        let quota = host.freq(codon) * n as f64;
        let floor = quota.floor();
        counts[codon as usize] = floor as u32;
        assigned += floor as u32;
        remainders.push((quota - floor, host.freq(codon), codon));
    }
    // Largest fractional part wins; ties go to the more frequent codon, then to the
    // lower index, so the result never depends on iteration order.
    remainders.sort_by(|a, b| {
        b.0.partial_cmp(&a.0)
            .unwrap()
            .then(b.1.partial_cmp(&a.1).unwrap())
            .then(a.2.cmp(&b.2))
    });
    for i in 0..(n - assigned) as usize {
        counts[remainders[i % remainders.len()].2 as usize] += 1;
    }
}

/// Build the codon multiset for a protein.
///
/// `gc_band` is the whole-gene GC requirement; the sliding-window GC constraints are
/// the arrangement's business, but the gene-wide fraction is fixed here because no
/// rearrangement can change it.
pub fn build(
    protein: &[AaIdx],
    host: &HostTable,
    max_relax: f64,
    gc_band: (f64, f64),
) -> Result<Composition, CompositionError> {
    if protein.is_empty() {
        return Err(CompositionError::Empty);
    }
    let mut group_total = [0u32; N_AA];
    for &aa in protein {
        group_total[aa as usize] += 1;
    }

    let mut counts = [0u32; N_CODONS];
    for aa in 0..N_AA as AaIdx {
        apportion(host, aa, group_total[aa as usize], &mut counts);
    }

    let mut gc_moves = 0;
    let mut gc = gc_fraction(&counts);
    let mut gc_in_band = true;
    while gc < gc_band.0 || gc > gc_band.1 {
        let want_more_gc = gc < gc_band.0;
        let Some(mv) = best_gc_move(&counts, host, want_more_gc) else {
            return Err(CompositionError::GcUnreachable {
                best: gc,
                band: gc_band,
            });
        };
        // Spending here is still spending: a move that would push the gene past
        // `max_relax` is not affordable, whatever it buys. Without this check the
        // composition alone could blow the limit before the search even started.
        let mut trial = counts;
        trial[mv.from as usize] -= 1;
        trial[mv.to as usize] += 1;
        if deviation(&trial, host) > max_relax {
            gc_in_band = false;
            break;
        }
        counts = trial;
        gc_moves += 1;
        gc = gc_fraction(&counts);
    }

    let deviation = deviation(&counts, host);
    Ok(Composition {
        counts,
        deviation,
        gc,
        cai: crate::metrics::cai_of_counts(&counts, host),
        budget_left: (max_relax - deviation).max(0.0),
        gc_moves,
        gc_in_band,
    })
}

struct GcMove {
    from: CodonIdx,
    to: CodonIdx,
}

/// Per-codon ideal counts and their total, which a within-group unit move never
/// changes. Precomputing them turns a move's deviation cost into O(1).
struct Ideals {
    ideal: [f64; N_CODONS],
    residues: f64,
}

impl Ideals {
    fn new(counts: &[u32; N_CODONS], host: &HostTable) -> Self {
        let mut group_total = [0u32; N_AA];
        for codon in 0..N_CODONS {
            group_total[codon_aa(codon as CodonIdx) as usize] += counts[codon];
        }
        let mut ideal = [0.0f64; N_CODONS];
        let mut residues = 0.0;
        for codon in 0..N_CODONS as CodonIdx {
            let n_a = group_total[codon_aa(codon) as usize] as f64;
            let want = round_half_even(host.freq(codon) * n_a);
            ideal[codon as usize] = want;
            residues += want;
        }
        Self { ideal, residues }
    }

    /// Change in deviation from moving one count from `from` to `to`.
    fn cost(&self, counts: &[u32; N_CODONS], from: CodonIdx, to: CodonIdx) -> f64 {
        if self.residues == 0.0 {
            return 0.0;
        }
        let at = |c: CodonIdx, delta: f64| {
            (counts[c as usize] as f64 + delta - self.ideal[c as usize]).abs()
                - (counts[c as usize] as f64 - self.ideal[c as usize]).abs()
        };
        (at(from, -1.0) + at(to, 1.0)) / (2.0 * self.residues)
    }
}

/// The synonymous unit move that buys the most GC change per unit of deviation.
fn best_gc_move(
    counts: &[u32; N_CODONS],
    host: &HostTable,
    want_more_gc: bool,
) -> Option<GcMove> {
    let ideals = Ideals::new(counts, host);
    let mut best: Option<(f64, GcMove)> = None;
    for from in 0..N_CODONS as CodonIdx {
        if counts[from as usize] == 0 {
            continue;
        }
        let aa = codon_aa(from);
        for to in host.usable_synonyms(aa) {
            if to == from {
                continue;
            }
            let gain = codon_gc(to) as i32 - codon_gc(from) as i32;
            if (want_more_gc && gain <= 0) || (!want_more_gc && gain >= 0) {
                continue;
            }
            let cost = ideals.cost(counts, from, to).max(0.0);
            // Prefer the biggest GC shift per unit of profile damage. The epsilon keeps
            // free moves (cost 0) from dividing by zero while still ranking first.
            let score = gain.abs() as f64 / (cost + 1e-9);
            if best.as_ref().map_or(true, |(s, _)| score > *s) {
                best = Some((score, GcMove { from, to }));
            }
        }
    }
    best.map(|(_, mv)| mv)
}

/// Limits on the bounded composition search. They are a budget, not a guarantee: a
/// search that runs out reports `exhausted` rather than claiming no other multiset
/// exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SearchLimits {
    /// States carried from one round to the next.
    pub beam_width: usize,
    /// Total states expanded. This is the only depth limit: the search keeps
    /// transferring counts until the budget is spent or the frontier empties.
    pub max_expansions: usize,
    /// Valid compositions handed on to arrangement.
    pub max_compositions: usize,
}

impl Default for SearchLimits {
    fn default() -> Self {
        Self { beam_width: 12, max_expansions: 2000, max_compositions: 8 }
    }
}

/// What a bounded composition search found.
#[derive(Debug, Clone)]
pub struct Exploration {
    /// Valid compositions, best first. Never empty: the baseline is always included.
    pub compositions: Vec<Composition>,
    /// States expanded, for reporting against the budget.
    pub expansions: usize,
    /// The budget ran out while the frontier still held states. This is search
    /// exhaustion, not proof that nothing better exists.
    pub exhausted: bool,
    /// The whole-gene GC band that was asked for.
    pub band: (f64, f64),
    /// The deviation a composition, and any arrangement of it, must stay within.
    /// This is `max_relax` whenever that is positive: a positive limit is a hard
    /// bound, and a protein whose rounding floor sits above it is an error rather
    /// than a quiet enlargement. For a frozen composition it is the rounding floor
    /// itself, which is the deviation the frozen multiset already has.
    pub dev_limit: f64,
}

/// How far a GC fraction sits outside its band; zero inside it.
pub fn gc_distance(gc: f64, band: (f64, f64)) -> f64 {
    if gc < band.0 {
        band.0 - gc
    } else if gc > band.1 {
        gc - band.1
    } else {
        0.0
    }
}

/// Comparisons against a hard limit allow this much floating-point slack. It is far
/// below any difference a single codon can make to either metric at gene length.
pub const LIMIT_TOLERANCE: f64 = 1e-9;

/// Several codon multisets for one protein, to be arranged and then ranked by CAI.
///
/// Every arrangement of one multiset has the same CAI, so searching arrangements alone
/// cannot optimise it. This walks synonymous one-count transfers from the host-profile
/// allocation and from a CAI-greedy allocation, keeping the compositions that stay
/// inside the deviation limit and the GC band.
///
/// Both limits are hard. A positive `max_relax` bounds every composition, including
/// the apportioned baseline: largest-remainder apportionment is L1-optimal, so when
/// even it deviates more, `DeviationFloor` says so rather than quietly raising the
/// bound to fit. `--max-relax 0` keeps its own meaning, freezing the baseline and
/// reporting the integer-rounding deviation it already has. The GC band is hard
/// throughout: a composition outside it is not a composition, and a protein that
/// cannot reach it within the deviation limit is `GcOutsideBudget`, never a widened
/// band. Only the CAI expectation is ever relaxed.
pub fn explore(
    protein: &[AaIdx],
    host: &HostTable,
    max_relax: f64,
    gc_band: (f64, f64),
    limits: &SearchLimits,
) -> Result<Exploration, CompositionError> {
    let baseline = build(protein, host, max_relax, gc_band)?;
    let frozen = max_relax <= 0.0;
    // No multiset deviates less than the baseline, so a positive limit below its
    // deviation is unsatisfiable, and saying so beats returning sequences that
    // overrun the budget the caller set.
    if !frozen && baseline.deviation > max_relax + LIMIT_TOLERANCE {
        return Err(CompositionError::DeviationFloor {
            best: baseline.deviation,
            limit: max_relax,
        });
    }
    let dev_limit = if frozen { baseline.deviation } else { max_relax };
    let max_compositions = limits.max_compositions.max(1);
    let in_band = |gc: f64| gc_distance(gc, gc_band) <= LIMIT_TOLERANCE;

    // `--max-relax 0` freezes the multiset. Arrangement still explores, but there is
    // exactly one composition to arrange, and nothing can move its GC.
    if frozen {
        if !in_band(baseline.gc) {
            return Err(CompositionError::GcOutsideBudget {
                best: baseline.gc,
                band: gc_band,
                limit: max_relax,
            });
        }
        return Ok(Exploration {
            compositions: vec![baseline],
            expansions: 0,
            exhausted: false,
            band: gc_band,
            dev_limit,
        });
    }

    let valid = |counts: &[u32; N_CODONS]| -> Option<Composition> {
        let dev = deviation(counts, host);
        let gc = gc_fraction(counts);
        if !dev.is_finite() || !gc.is_finite() {
            return None;
        }
        if dev > dev_limit + LIMIT_TOLERANCE || !in_band(gc) {
            return None;
        }
        Some(Composition {
            counts: *counts,
            deviation: dev,
            gc,
            cai: crate::metrics::cai_of_counts(counts, host),
            budget_left: (max_relax - dev).max(0.0),
            gc_moves: 0,
            gc_in_band: gc_distance(gc, gc_band) == 0.0,
        })
    };

    // Rank: first how far a state is from being usable at all, then CAI, then profile
    // deviation. The count vector breaks the last tie so the order never depends on
    // discovery order.
    let score = |counts: &[u32; N_CODONS]| -> (f64, f64, f64) {
        let dev = deviation(counts, host);
        let gc = gc_fraction(counts);
        let slack = (dev - dev_limit).max(0.0) + gc_distance(gc, gc_band);
        (slack, -crate::metrics::cai_of_counts(counts, host), dev)
    };
    let better = |a: &(f64, f64, f64, [u32; N_CODONS]), b: &(f64, f64, f64, [u32; N_CODONS])| {
        a.0.partial_cmp(&b.0)
            .unwrap()
            .then(a.1.partial_cmp(&b.1).unwrap())
            .then(a.2.partial_cmp(&b.2).unwrap())
            .then(a.3.cmp(&b.3))
    };

    let mut seen: HashSet<[u32; N_CODONS]> = HashSet::new();
    let mut found: Vec<Composition> = Vec::new();
    // The GC closest to the band among multisets inside the deviation limit. It is
    // what the error reports when nothing reaches the band, so the caller learns how
    // far off the request was rather than only that it failed.
    let mut closest_gc: Option<f64> = None;
    let consider = |counts: &[u32; N_CODONS],
                    found: &mut Vec<Composition>,
                    closest: &mut Option<f64>| {
        if deviation(counts, host) > dev_limit + LIMIT_TOLERANCE {
            return;
        }
        let gc = gc_fraction(counts);
        let nearer = match *closest {
            Some(best) => gc_distance(gc, gc_band) < gc_distance(best, gc_band),
            None => true,
        };
        if nearer {
            *closest = Some(gc);
        }
        if let Some(composition) = valid(counts) {
            found.push(composition);
        }
    };

    let mut frontier: Vec<[u32; N_CODONS]> = Vec::new();
    for start in [baseline.counts, cai_greedy(protein, host)] {
        if seen.insert(start) {
            consider(&start, &mut found, &mut closest_gc);
            frontier.push(start);
        }
    }

    let mut expansions = 0usize;
    let mut exhausted = false;
    while !frontier.is_empty() && expansions < limits.max_expansions {
        let mut next: Vec<(f64, f64, f64, [u32; N_CODONS])> = Vec::new();
        for counts in std::mem::take(&mut frontier) {
            if expansions >= limits.max_expansions {
                exhausted = true;
                break;
            }
            expansions += 1;
            for child in transfers(&counts, host) {
                if !seen.insert(child) {
                    continue;
                }
                consider(&child, &mut found, &mut closest_gc);
                let s = score(&child);
                next.push((s.0, s.1, s.2, child));
                // A transfer that buys CAI often pushes GC out of band. Offer the
                // paired state that pays the cheapest synonymous move back, so the
                // beam is not stopped by a barrier one transfer wide.
                if !in_band(gc_fraction(&child)) {
                    let want_more_gc = gc_fraction(&child) < gc_band.0;
                    if let Some(mv) = best_gc_move(&child, host, want_more_gc) {
                        let mut paired = child;
                        paired[mv.from as usize] -= 1;
                        paired[mv.to as usize] += 1;
                        if seen.insert(paired) {
                            consider(&paired, &mut found, &mut closest_gc);
                            let s = score(&paired);
                            next.push((s.0, s.1, s.2, paired));
                        }
                    }
                }
            }
        }
        if expansions >= limits.max_expansions && !next.is_empty() {
            exhausted = true;
        }
        next.sort_by(better);
        next.truncate(limits.beam_width.max(1));
        frontier = next.into_iter().map(|(_, _, _, counts)| counts).collect();
        // Retained states are capped too: with a generous deviation limit most of the
        // frontier's children are valid, and keeping every one of them would grow the
        // memory with the budget rather than with what is returned.
        if found.len() > 4 * max_compositions {
            prune(&mut found, max_compositions);
        }
        if exhausted {
            break;
        }
    }

    prune(&mut found, max_compositions);
    // The baseline is kept when it is valid, since arrangement can fail for the
    // high-CAI multiset and succeed for it. It is not forced in when it is not:
    // an invalid composition is not a fallback, it is a limit violation.
    let baseline_gc = baseline.gc;
    if in_band(baseline.gc)
        && baseline.deviation <= dev_limit + LIMIT_TOLERANCE
        && !found.iter().any(|c| c.counts == baseline.counts)
    {
        found.truncate(max_compositions.saturating_sub(1));
        found.push(baseline);
    }
    if found.is_empty() {
        // Deviation alone is satisfiable, since the baseline cleared it above, so
        // the band is what nothing inside the limit could reach.
        return Err(CompositionError::GcOutsideBudget {
            best: closest_gc.unwrap_or(baseline_gc),
            band: gc_band,
            limit: max_relax,
        });
    }

    Ok(Exploration { compositions: found, expansions, exhausted, band: gc_band, dev_limit })
}

/// Rank by descending CAI, then ascending deviation, then the count vector, drop
/// repeats, and keep the best `keep`. Every key is a number or the counts themselves,
/// so the order does not depend on the order states were discovered in.
fn prune(found: &mut Vec<Composition>, keep: usize) {
    found.sort_by(|a, b| {
        b.cai
            .partial_cmp(&a.cai)
            .unwrap()
            .then(a.deviation.partial_cmp(&b.deviation).unwrap())
            .then(a.counts.cmp(&b.counts))
    });
    found.dedup_by(|a, b| a.counts == b.counts);
    found.truncate(keep);
}

/// Every multiset one synonymous count transfer away.
fn transfers(counts: &[u32; N_CODONS], host: &HostTable) -> Vec<[u32; N_CODONS]> {
    let mut out = Vec::new();
    for from in 0..N_CODONS as CodonIdx {
        if counts[from as usize] == 0 {
            continue;
        }
        for to in host.usable_synonyms(codon_aa(from)) {
            if to == from {
                continue;
            }
            let mut child = *counts;
            child[from as usize] -= 1;
            child[to as usize] += 1;
            out.push(child);
        }
    }
    out
}

/// Every residue on its host-preferred codon: the CAI ceiling, usually far outside the
/// deviation limit, and useful only as the other end of the search.
fn cai_greedy(protein: &[AaIdx], host: &HostTable) -> [u32; N_CODONS] {
    let mut best: [Option<CodonIdx>; N_AA] = [None; N_AA];
    for aa in 0..N_AA as AaIdx {
        // Ties go to the lower codon index, so the allocation never depends on the
        // order `usable_synonyms` happens to yield.
        let mut pick = (-1.0f64, 0 as CodonIdx, false);
        for codon in host.usable_synonyms(aa) {
            let weight = host.weight(codon);
            if !pick.2 || weight > pick.0 || (weight == pick.0 && codon < pick.1) {
                pick = (weight, codon, true);
            }
        }
        best[aa as usize] = pick.2.then_some(pick.1);
    }
    let mut counts = [0u32; N_CODONS];
    for &aa in protein {
        if let Some(codon) = best[aa as usize] {
            counts[codon as usize] += 1;
        }
    }
    counts
}

/// Human-readable summary of which codons the composition uses, for reporting.
pub fn summary(counts: &[u32; N_CODONS]) -> String {
    let mut parts: Vec<String> = Vec::new();
    for codon in 0..N_CODONS as CodonIdx {
        if counts[codon as usize] > 0 {
            parts.push(format!("{}:{}", codon_string(codon), counts[codon as usize]));
        }
    }
    parts.join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::aa_index;
    use crate::host::HostTable;
    use crate::rng::Rng;
    use crate::testdata::{ECOLI_TABLE, GFP};
    use std::path::Path;

    fn ecoli() -> HostTable {
        HostTable::load(Path::new(ECOLI_TABLE), 0.10).unwrap()
    }

    fn protein(seq: &str) -> Vec<AaIdx> {
        seq.bytes().map(|b| aa_index(b).unwrap()).collect()
    }

    const WIDE_BAND: (f64, f64) = (0.0, 1.0);

    #[test]
    fn counts_sum_to_group_size_for_random_proteins() {
        let host = ecoli();
        let mut rng = Rng::new(1);
        let letters = b"ACDEFGHIKLMNPQRSTVWY";
        for n in [1usize, 2, 7, 23, 200, 999] {
            let seq: Vec<AaIdx> = (0..n)
                .map(|_| aa_index(letters[rng.below(letters.len())]).unwrap())
                .collect();
            let comp = build(&seq, &host, 0.1, WIDE_BAND).unwrap();
            let mut group_total = [0u32; N_AA];
            for &aa in &seq {
                group_total[aa as usize] += 1;
            }
            for aa in 0..N_AA as AaIdx {
                let used: u32 = host
                    .synonyms(aa)
                    .iter()
                    .map(|&c| comp.counts[c as usize])
                    .sum();
                assert_eq!(
                    used, group_total[aa as usize],
                    "amino acid {} with n={n}",
                    crate::code::aa_letter(aa) as char
                );
            }
            // and no dropped codon is ever used
            for codon in 0..N_CODONS as CodonIdx {
                if host.freq(codon) == 0.0 {
                    assert_eq!(comp.counts[codon as usize], 0);
                }
            }
        }
    }

    /// What the Python effectively builds: round each codon, then add or remove units
    /// from the most- and least-used codons until the group total is right.
    fn python_style_allocation(host: &HostTable, aa: AaIdx, n: u32) -> Vec<(CodonIdx, i64)> {
        let usable: Vec<CodonIdx> = host.usable_synonyms(aa).collect();
        let mut counts: Vec<(CodonIdx, i64)> = usable
            .iter()
            .map(|&c| (c, round_half_even(host.freq(c) * n as f64) as i64))
            .collect();
        let total = |cs: &Vec<(CodonIdx, i64)>| cs.iter().map(|(_, k)| *k).sum::<i64>();
        while total(&counts) > n as i64 {
            let pick = counts
                .iter_mut()
                .filter(|(_, k)| *k > 0)
                .min_by(|a, b| host.freq(a.0).partial_cmp(&host.freq(b.0)).unwrap())
                .unwrap();
            pick.1 -= 1;
        }
        while total(&counts) < n as i64 {
            let pick = counts
                .iter_mut()
                .max_by(|a, b| host.freq(a.0).partial_cmp(&host.freq(b.0)).unwrap())
                .unwrap();
            pick.1 += 1;
        }
        counts
    }

    #[test]
    fn no_unit_move_improves_the_apportionment() {
        let host = ecoli();
        let mut rng = Rng::new(5);
        let letters = b"ACDEFGHIKLMNPQRSTVWY";
        for _ in 0..30 {
            let n = rng.between(10, 300);
            let seq: Vec<AaIdx> = (0..n)
                .map(|_| aa_index(letters[rng.below(letters.len())]).unwrap())
                .collect();
            let comp = build(&seq, &host, 0.1, WIDE_BAND).unwrap();
            let mut group_total = [0u32; N_AA];
            for &aa in &seq {
                group_total[aa as usize] += 1;
            }
            let err = |c: CodonIdx, k: f64| {
                let n_a = group_total[codon_aa(c) as usize] as f64;
                (k - host.freq(c) * n_a).abs()
            };
            // Largest remainder is L1-optimal, so no transfer of one unit between two
            // synonyms can reduce the total distance to the ideal.
            for aa in 0..N_AA as AaIdx {
                let usable: Vec<CodonIdx> = host.usable_synonyms(aa).collect();
                for &x in &usable {
                    if comp.counts[x as usize] == 0 {
                        continue;
                    }
                    for &y in &usable {
                        if x == y {
                            continue;
                        }
                        let (kx, ky) =
                            (comp.counts[x as usize] as f64, comp.counts[y as usize] as f64);
                        let before = err(x, kx) + err(y, ky);
                        let after = err(x, kx - 1.0) + err(y, ky + 1.0);
                        assert!(
                            after >= before - 1e-9,
                            "moving one unit {} -> {} would improve L1",
                            codon_string(x),
                            codon_string(y)
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn apportionment_beats_rounding_on_l1_distance() {
        let host = ecoli();
        let mut rng = Rng::new(2);
        let letters = b"ACDEFGHIKLMNPQRSTVWY";
        for _ in 0..40 {
            let n = rng.between(10, 300);
            let seq: Vec<AaIdx> = (0..n)
                .map(|_| aa_index(letters[rng.below(letters.len())]).unwrap())
                .collect();
            let comp = build(&seq, &host, 0.1, WIDE_BAND).unwrap();
            let mut group_total = [0u32; N_AA];
            for &aa in &seq {
                group_total[aa as usize] += 1;
            }
            let (mut ours, mut theirs) = (0.0f64, 0.0f64);
            for aa in 0..N_AA as AaIdx {
                let n_a = group_total[aa as usize];
                if n_a == 0 {
                    continue;
                }
                for &codon in &host.usable_synonyms(aa).collect::<Vec<_>>() {
                    let target = host.freq(codon) * n_a as f64;
                    ours += (comp.counts[codon as usize] as f64 - target).abs();
                }
                for (codon, k) in python_style_allocation(&host, aa, n_a) {
                    theirs += (k as f64 - host.freq(codon) * n_a as f64).abs();
                }
            }
            assert!(
                ours <= theirs + 1e-9,
                "largest remainder {ours} should beat the round-and-fix allocation {theirs}"
            );
        }
    }

    #[test]
    fn gfp_deviation_is_under_one_percent() {
        let host = ecoli();
        let comp = build(&protein(GFP), &host, 0.1, (0.30, 0.65)).unwrap();
        assert!(
            comp.deviation < 0.01,
            "deviation {} should be under 1%",
            comp.deviation
        );
        // the prototype measured 0.0083 with the same table and protein
        assert!((comp.deviation - 0.0083).abs() < 0.002, "{}", comp.deviation);
        assert!(comp.budget_left > 0.09);
        assert_eq!(comp.gc_moves, 0, "GFP needs no GC correction");
    }

    #[test]
    fn gc_adjustment_pulls_into_the_band() {
        let host = ecoli();
        let seq = protein(GFP);
        let natural = build(&seq, &host, 1.0, WIDE_BAND).unwrap();
        assert!(
            (0.45..0.58).contains(&natural.gc),
            "setup expects GFP's natural GC near 0.5, got {}",
            natural.gc
        );

        // squeeze from above and from below; both must land inside the band
        for band in [(0.30, 0.45), (0.58, 0.70)] {
            let fixed = build(&seq, &host, 1.0, band).unwrap();
            assert!(
                fixed.gc >= band.0 && fixed.gc <= band.1,
                "GC {} should be inside {band:?}",
                fixed.gc
            );
            assert!(fixed.gc_moves > 0);
            assert!(
                fixed.deviation > natural.deviation,
                "pulling GC into {band:?} costs profile accuracy"
            );
        }
        let fixed = build(&seq, &host, 1.0, (0.30, 0.45)).unwrap();
        // the protein is still encodable: every group still sums correctly
        let mut group_total = [0u32; N_AA];
        for &aa in &seq {
            group_total[aa as usize] += 1;
        }
        for aa in 0..N_AA as AaIdx {
            let used: u32 = host
                .synonyms(aa)
                .iter()
                .map(|&c| fixed.counts[c as usize])
                .sum();
            assert_eq!(used, group_total[aa as usize]);
        }
    }

    #[test]
    fn gc_adjustment_never_spends_past_max_relax() {
        let host = ecoli();
        let seq = protein(GFP);
        // A band GFP's natural composition misses by a long way, so the correction
        // wants to spend far more than the limit allows.
        // The apportionment itself has a floor: the best integer allocation of this
        // protein is already 0.0083 away from the host profile, so a limit below that
        // cannot be met by anything. What must hold is that the GC correction never
        // spends beyond the limit, and spends nothing at all when the floor is above it.
        let floor = build(&seq, &host, 1.0, WIDE_BAND).unwrap().deviation;
        for band in [(0.30, 0.40), (0.60, 0.70)] {
            for max_relax in [0.0, 0.02, 0.10] {
                let comp = build(&seq, &host, max_relax, band).unwrap();
                assert!(
                    comp.deviation <= floor.max(max_relax) + 1e-9,
                    "band {band:?}, max_relax {max_relax}: deviation {} exceeds the limit",
                    comp.deviation
                );
                if floor >= max_relax {
                    assert_eq!(comp.gc_moves, 0, "nothing was affordable, so nothing should be spent");
                }
                assert!(!comp.gc_in_band, "this band is meant to be unaffordable");
            }
        }
        // and with budget to spare it still reaches the band
        let comp = build(&seq, &host, 1.0, (0.30, 0.40)).unwrap();
        assert!(comp.gc_in_band);
        assert!(comp.gc <= 0.40 && comp.gc >= 0.30);
    }

    #[test]
    fn gc_adjustment_reports_what_it_cannot_reach() {
        let host = ecoli();
        // Trp has exactly one codon, TGG, so GC is pinned at 2/3 with no way to move it.
        let seq = protein("WWWWWWWWWW");
        let err = build(&seq, &host, 1.0, (0.30, 0.40)).unwrap_err();
        match err {
            CompositionError::GcUnreachable { best, .. } => {
                assert!((best - 2.0 / 3.0).abs() < 1e-9, "best was {best}")
            }
            other => panic!("expected GcUnreachable, got {other:?}"),
        }
    }

    #[test]
    fn empty_protein_is_an_error() {
        assert!(matches!(
            build(&[], &ecoli(), 0.1, WIDE_BAND),
            Err(CompositionError::Empty)
        ));
    }

    #[test]
    fn deviation_is_zero_when_counts_match_the_ideal_exactly() {
        let host = ecoli();
        // 100 alanines: E. coli fractions .16/.27/.21/.36 are exact hundredths here.
        let comp = build(&protein(&"A".repeat(100)), &host, 0.1, WIDE_BAND).unwrap();
        assert_eq!(comp.deviation, 0.0);
        assert_eq!(comp.counts[crate::code::parse_codon("GCG").unwrap() as usize], 36);
        assert_eq!(comp.counts[crate::code::parse_codon("GCC").unwrap() as usize], 27);
        assert_eq!(comp.counts[crate::code::parse_codon("GCA").unwrap() as usize], 21);
        assert_eq!(comp.counts[crate::code::parse_codon("GCT").unwrap() as usize], 16);
    }

    #[test]
    fn stop_codons_are_handled_like_any_other_amino_acid() {
        let host = ecoli();
        let seq = protein(&format!("{GFP}*"));
        let comp = build(&seq, &host, 0.1, (0.30, 0.65)).unwrap();
        let stops: u32 = host
            .synonyms(crate::code::STOP)
            .iter()
            .map(|&c| comp.counts[c as usize])
            .sum();
        assert_eq!(stops, 1);
    }

    fn group_totals(counts: &[u32; N_CODONS]) -> [u32; N_AA] {
        let mut totals = [0u32; N_AA];
        for codon in 0..N_CODONS as CodonIdx {
            totals[codon_aa(codon) as usize] += counts[codon as usize];
        }
        totals
    }

    #[test]
    fn every_explored_composition_encodes_the_same_protein() {
        let host = ecoli();
        let seq = protein(GFP);
        let found = explore(&seq, &host, 0.1, (0.30, 0.65), &SearchLimits::default()).unwrap();
        let baseline = build(&seq, &host, 0.1, (0.30, 0.65)).unwrap();
        assert!(found.compositions.len() > 1, "the search should offer alternatives");
        for composition in &found.compositions {
            assert_eq!(
                group_totals(&composition.counts),
                group_totals(&baseline.counts),
                "a transfer may only move a count within its synonymous group"
            );
            assert!(composition.deviation <= found.dev_limit + LIMIT_TOLERANCE);
            assert_eq!(
                gc_distance(composition.gc, found.band),
                0.0,
                "the requested band is hard"
            );
            for codon in 0..N_CODONS as CodonIdx {
                if composition.counts[codon as usize] > 0 {
                    assert!(host.freq(codon) > 0.0, "only codons the host uses");
                }
            }
        }
    }

    #[test]
    fn the_baseline_is_always_offered_and_the_set_is_ranked_by_cai() {
        let host = ecoli();
        let seq = protein(GFP);
        let baseline = build(&seq, &host, 0.1, (0.30, 0.65)).unwrap();
        let found = explore(&seq, &host, 0.1, (0.30, 0.65), &SearchLimits::default()).unwrap();
        assert!(
            found.compositions.iter().any(|c| c.counts == baseline.counts),
            "arrangement can fail for the high-CAI multiset and succeed for the baseline"
        );
        for pair in found.compositions.windows(2) {
            if pair[1].counts == baseline.counts {
                continue; // the baseline is kept whatever its rank
            }
            assert!(pair[0].cai >= pair[1].cai - LIMIT_TOLERANCE);
        }
        assert!(
            found.compositions[0].cai > baseline.cai,
            "spending deviation should buy CAI: {} vs {}",
            found.compositions[0].cai,
            baseline.cai
        );
    }

    #[test]
    fn a_zero_limit_freezes_the_composition() {
        let host = ecoli();
        let seq = protein(GFP);
        let found = explore(&seq, &host, 0.0, (0.30, 0.65), &SearchLimits::default()).unwrap();
        assert_eq!(found.compositions.len(), 1);
        assert_eq!(found.expansions, 0);
        assert!(!found.exhausted);
        assert_eq!(found.compositions[0].counts, build(&seq, &host, 0.0, (0.30, 0.65)).unwrap().counts);
    }

    #[test]
    fn a_larger_limit_reaches_a_higher_cai() {
        let host = ecoli();
        let seq = protein(GFP);
        let tight = explore(&seq, &host, 0.02, (0.30, 0.65), &SearchLimits::default()).unwrap();
        let loose = explore(&seq, &host, 0.20, (0.30, 0.65), &SearchLimits::default()).unwrap();
        assert!(
            loose.compositions[0].cai > tight.compositions[0].cai,
            "{} should beat {}",
            loose.compositions[0].cai,
            tight.compositions[0].cai
        );
    }

    #[test]
    fn a_spent_budget_is_reported_as_exhaustion_not_impossibility() {
        let host = ecoli();
        let seq = protein(GFP);
        let limits = SearchLimits { max_expansions: 5, ..SearchLimits::default() };
        let found = explore(&seq, &host, 0.2, (0.30, 0.65), &limits).unwrap();
        assert!(found.exhausted, "the budget was spent with states still to expand");
        assert!(found.expansions <= 5);
        assert!(!found.compositions.is_empty(), "what it did find is still returned");
    }

    #[test]
    fn exploration_is_reproducible() {
        let host = ecoli();
        let seq = protein(GFP);
        let once = explore(&seq, &host, 0.1, (0.30, 0.65), &SearchLimits::default()).unwrap();
        let twice = explore(&seq, &host, 0.1, (0.30, 0.65), &SearchLimits::default()).unwrap();
        let counts = |e: &Exploration| e.compositions.iter().map(|c| c.counts).collect::<Vec<_>>();
        assert_eq!(counts(&once), counts(&twice));
    }

    #[test]
    fn a_positive_limit_below_the_rounding_floor_is_refused() {
        let host = ecoli();
        let seq = protein(GFP);
        let floor = build(&seq, &host, 0.1, (0.30, 0.65)).unwrap().deviation;
        assert!(floor > 0.0, "integer rounding costs something");
        let err = explore(&seq, &host, floor / 2.0, (0.30, 0.65), &SearchLimits::default())
            .unwrap_err();
        match err {
            CompositionError::DeviationFloor { best, limit } => {
                assert!((best - floor).abs() < 1e-12);
                assert!((limit - floor / 2.0).abs() < 1e-12);
            }
            other => panic!("expected the floor to be reported, got {other:?}"),
        }
    }

    #[test]
    fn a_band_out_of_reach_within_the_limit_is_refused() {
        let host = ecoli();
        let seq = protein(GFP);
        let err = explore(&seq, &host, 0.01, (0.30, 0.40), &SearchLimits::default())
            .unwrap_err();
        match err {
            CompositionError::GcOutsideBudget { best, band, limit } => {
                assert_eq!(band, (0.30, 0.40));
                assert!((limit - 0.01).abs() < 1e-12);
                assert!(best > 0.40, "it should say how close it got: {best}");
            }
            other => panic!("expected the band to be reported, got {other:?}"),
        }
    }

    #[test]
    fn a_frozen_composition_outside_the_band_is_refused_too() {
        let host = ecoli();
        let seq = protein(GFP);
        // nothing can move a frozen multiset into a band it does not already meet
        let err =
            explore(&seq, &host, 0.0, (0.30, 0.40), &SearchLimits::default()).unwrap_err();
        assert!(matches!(err, CompositionError::GcOutsideBudget { .. }), "{err:?}");
    }

    #[test]
    fn a_frozen_composition_reports_its_rounding_deviation_as_the_bound() {
        let host = ecoli();
        let seq = protein(GFP);
        let found = explore(&seq, &host, 0.0, (0.30, 0.65), &SearchLimits::default()).unwrap();
        assert_eq!(found.compositions.len(), 1);
        assert_eq!(found.dev_limit, found.compositions[0].deviation);
        assert!(found.dev_limit > 0.0, "the rounding cost is reported, not hidden");
    }
}
