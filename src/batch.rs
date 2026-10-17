//! Per-record orchestration: several candidates, ranked by CAI, with one retry.
//!
//! The engine optimises one protein into one gene. This layer asks it for a set of
//! distinct genes instead, ranks them by codon adaptation index, and, when none of
//! them clears the quality threshold, runs the whole thing a second time with more
//! deviation allowed and more search effort, then returns the best of both runs.
//!
//! The threshold is a retry trigger and a quality warning, never a filter: a complete
//! set whose best CAI is low is still returned, with the shortfall reported.

use crate::arrange::search::{optimise_all, Options, Outcome};
use crate::code::AaIdx;
use crate::composition::{self, gc_distance, Exploration, SearchLimits, LIMIT_TOLERANCE};
use crate::host::HostTable;

/// Settings the retry run changes. Everything not named here is the first run's.
#[derive(Debug, Clone)]
pub struct RetryPolicy {
    /// Deviation from the host profile the retry may spend. Must be at least the
    /// first run's limit; zero keeps the composition frozen.
    pub max_relax: f64,
    /// Arrangements tried per composition.
    pub restarts: usize,
    /// States the composition search may expand.
    pub expansions: usize,
}

impl RetryPolicy {
    /// The defaults the CLI documents: ten more points of deviation for a limit that
    /// is already positive, twice the restarts, twice the composition budget. A
    /// frozen composition stays frozen; only arrangement effort grows.
    pub fn doubling(max_relax: f64, restarts: usize, expansions: usize) -> Self {
        let relax = if max_relax > 0.0 { (max_relax + 0.10).min(1.0) } else { 0.0 };
        Self {
            max_relax: relax,
            restarts: restarts.saturating_mul(2).max(1),
            expansions: expansions.saturating_mul(2).max(1),
        }
    }
}

#[derive(Debug, Clone)]
pub struct BatchOptions {
    /// Candidates wanted per protein.
    pub wanted: usize,
    /// CAI a candidate must strictly exceed for the first run to stand. Equality does
    /// not pass.
    pub cai_threshold: f64,
    pub max_relax: f64,
    pub gc_band: (f64, f64),
    pub limits: SearchLimits,
    pub retry: RetryPolicy,
    /// Arrangement settings. `seed` is the base of every derived stream; `restarts`
    /// is the first run's count.
    pub arrange: Options,
}

impl Default for BatchOptions {
    fn default() -> Self {
        let arrange = Options::default();
        let limits = SearchLimits::default();
        Self {
            wanted: 5,
            cai_threshold: 0.80,
            max_relax: 0.10,
            gc_band: (0.30, 0.65),
            retry: RetryPolicy::doubling(0.10, arrange.restarts, limits.max_expansions),
            limits,
            arrange,
        }
    }
}

/// One returned gene, with enough provenance to say where it came from and which
/// bound it was built under.
#[derive(Debug, Clone)]
pub struct Candidate {
    pub dna: String,
    pub cai: f64,
    pub deviation: f64,
    pub gc: f64,
    /// 1 for the first run, 2 for the retry.
    pub run: u8,
    pub composition: usize,
    pub restart: usize,
    pub seed: u64,
    /// The deviation the run asked to be allowed. A retry candidate may legitimately
    /// exceed the first run's.
    pub max_relax: f64,
    /// The deviation bound actually applied. It is `max_relax`, except that the
    /// largest-remainder baseline is always admitted: no multiset deviates less, so a
    /// short protein whose integer rounding alone exceeds the limit is still built,
    /// against this bound rather than a silently enlarged one.
    pub dev_limit: f64,
    pub moves: usize,
    pub substitutions: usize,
}

/// What one run cost and found.
#[derive(Debug, Clone)]
pub struct Attempt {
    pub run: u8,
    pub max_relax: f64,
    pub restarts: usize,
    /// The deviation bound actually applied; see `Candidate::dev_limit`.
    pub dev_limit: f64,
    pub expansion_budget: usize,
    pub compositions: usize,
    pub expansions: usize,
    /// The composition search hit its budget. Not proof that no better multiset or
    /// further distinct encoding exists.
    pub search_exhausted: bool,
    pub arrangements: usize,
    pub valid: usize,
    pub best_cai: Option<f64>,
    /// Why the arrangements that were turned away were turned away. Without this a
    /// record with no candidate says only that it has none.
    pub rejected: Rejections,
}

/// Counts and one example of each reason an arrangement was not a candidate.
#[derive(Debug, Clone, Default)]
pub struct Rejections {
    /// Arrangements that still held sequence constraint violations, and the fewest
    /// any of them left, with the kinds involved.
    pub unresolved: usize,
    pub fewest_violations: Option<u32>,
    pub violation_detail: Option<String>,
    /// Arrangements whose profile deviation or whole-gene GC left the hard bound.
    pub over_limit: usize,
    /// Arrangements that did not translate back, or whose metrics were not finite.
    /// Either is an engine bug rather than a hard case.
    pub unsound: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// The full set, with at least one candidate above the threshold.
    Ok,
    /// The full set, but nothing cleared the threshold. A successful fallback.
    LowCai,
    /// Fewer distinct valid genes than asked for.
    InsufficientCandidates,
    /// No valid gene at all.
    NoValidCandidates,
    /// The protein or its composition could not be processed.
    Error,
}

impl Status {
    pub fn label(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::LowCai => "low_cai",
            Self::InsufficientCandidates => "insufficient_candidates",
            Self::NoValidCandidates => "no_valid_candidates",
            Self::Error => "error",
        }
    }

    /// Whether this status alone should make the process exit non-zero. A complete
    /// set with low CAI is a warning, not a failure.
    pub fn is_failure(self) -> bool {
        matches!(
            self,
            Self::InsufficientCandidates | Self::NoValidCandidates | Self::Error
        )
    }
}

#[derive(Debug, Clone)]
pub struct RecordReport {
    pub index: usize,
    pub status: Status,
    pub candidates: Vec<Candidate>,
    pub attempts: Vec<Attempt>,
    pub retry_reason: Option<String>,
    /// True when the best candidate did not clear the threshold, whatever else is
    /// also wrong. Kept separate from `status` so several conditions stay visible.
    pub low_cai: bool,
    /// True when any run's composition search ran out of budget.
    pub search_exhausted: bool,
    pub error: Option<String>,
}

impl RecordReport {
    pub fn best_cai(&self) -> Option<f64> {
        self.candidates.first().map(|c| c.cai)
    }
}

/// Mix stable indices into a stream seed. Index arithmetic only: no hashing of names,
/// no thread identity, so the same seed reproduces a run whatever the thread count.
pub fn stream_seed(base: u64, record: usize, run: u8, composition: usize) -> u64 {
    let mut z = base
        ^ (record as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ (run as u64).wrapping_mul(0xD1B5_4A32_D192_ED03)
        ^ (composition as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Optimise one protein into a ranked candidate set, retrying at most once.
pub fn run_record(
    index: usize,
    protein: &[AaIdx],
    host: &HostTable,
    options: &BatchOptions,
) -> RecordReport {
    let wanted = options.wanted.max(1);
    let mut report = RecordReport {
        index,
        status: Status::NoValidCandidates,
        candidates: Vec::new(),
        attempts: Vec::new(),
        retry_reason: None,
        low_cai: false,
        search_exhausted: false,
        error: None,
    };

    let mut pool: Vec<Candidate> = Vec::new();
    for run in 1u8..=2 {
        let (max_relax, restarts, expansions) = if run == 1 {
            (options.max_relax, options.arrange.restarts.max(1), options.limits.max_expansions)
        } else {
            (options.retry.max_relax, options.retry.restarts.max(1), options.retry.expansions)
        };
        let limits = SearchLimits { max_expansions: expansions, ..options.limits };

        let exploration =
            match composition::explore(protein, host, max_relax, options.gc_band, &limits) {
                Ok(exploration) => exploration,
                Err(error) => {
                    report.status = Status::Error;
                    report.error = Some(error.to_string());
                    return report;
                }
            };

        let attempt = collect(
            index, protein, host, options, &exploration, run, max_relax, restarts, expansions,
            &mut pool,
        );
        report.search_exhausted |= attempt.search_exhausted;
        report.attempts.push(attempt);

        rank(&mut pool);
        let qualifies = pool.first().is_some_and(|c| c.cai > options.cai_threshold);
        if qualifies {
            break;
        }
        if run == 1 {
            report.retry_reason = Some(match pool.first() {
                Some(best) => format!(
                    "best CAI {:.6} does not exceed the threshold {:.6}",
                    best.cai, options.cai_threshold
                ),
                None => "the first run produced no valid candidate".to_string(),
            });
        }
    }

    rank(&mut pool);
    pool.truncate(wanted);
    report.low_cai = !pool.first().is_some_and(|c| c.cai > options.cai_threshold);
    report.status = if pool.is_empty() {
        Status::NoValidCandidates
    } else if pool.len() < wanted {
        Status::InsufficientCandidates
    } else if report.low_cai {
        Status::LowCai
    } else {
        Status::Ok
    };
    report.candidates = pool;
    report
}

/// Arrange every composition the search offered and add the valid genes to the pool.
#[allow(clippy::too_many_arguments)]
fn collect(
    index: usize,
    protein: &[AaIdx],
    host: &HostTable,
    options: &BatchOptions,
    exploration: &Exploration,
    run: u8,
    max_relax: f64,
    restarts: usize,
    expansion_budget: usize,
    pool: &mut Vec<Candidate>,
) -> Attempt {
    let mut attempt = Attempt {
        run,
        max_relax,
        restarts,
        dev_limit: exploration.dev_limit,
        expansion_budget,
        compositions: exploration.compositions.len(),
        expansions: exploration.expansions,
        search_exhausted: exploration.exhausted,
        arrangements: 0,
        valid: 0,
        best_cai: None,
        rejected: Rejections::default(),
    };

    // Every composition the search offered is arranged, so the candidate set depends
    // on the search budget and not on the threshold: a low threshold must not buy a
    // worse set of genes than a high one.
    for (composition_index, composition) in exploration.compositions.iter().enumerate() {
        let seed = stream_seed(options.arrange.seed, index, run, composition_index);
        let arrange = Options { seed, restarts, ..options.arrange.clone() };
        for outcome in optimise_all(protein, composition, host, &arrange) {
            attempt.arrangements += 1;
            let candidate = match accept(
                protein, &outcome, exploration, run, composition_index, seed, max_relax,
            ) {
                Ok(candidate) => candidate,
                Err(reason) => {
                    note(&mut attempt.rejected, reason, &outcome);
                    continue;
                }
            };
            attempt.valid += 1;
            attempt.best_cai = Some(match attempt.best_cai {
                Some(best) if best >= candidate.cai => best,
                _ => candidate.cai,
            });
            pool.push(candidate);
        }
    }
    attempt
}

/// Why an arrangement was not admitted.
enum Reject {
    Unresolved,
    OverLimit,
    Unsound,
}

fn note(rejected: &mut Rejections, reason: Reject, outcome: &Outcome) {
    match reason {
        Reject::Unresolved => {
            rejected.unresolved += 1;
            let closest = match rejected.fewest_violations {
                Some(fewest) => outcome.violations < fewest,
                None => true,
            };
            if closest {
                rejected.fewest_violations = Some(outcome.violations);
                rejected.violation_detail = Some(
                    outcome
                        .by_kind
                        .iter()
                        .map(|(kind, n)| format!("{} x{n}", kind.label()))
                        .collect::<Vec<_>>()
                        .join(", "),
                );
            }
        }
        Reject::OverLimit => rejected.over_limit += 1,
        Reject::Unsound => rejected.unsound += 1,
    }
}

/// Re-derive everything from the finished sequence and admit it only if every hard
/// requirement holds. CAI is not one of them: a low-CAI gene is still valid, and is
/// what the fallback returns.
fn accept(
    protein: &[AaIdx],
    outcome: &Outcome,
    exploration: &Exploration,
    run: u8,
    composition: usize,
    seed: u64,
    max_relax: f64,
) -> Result<Candidate, Reject> {
    if outcome.violations > 0 {
        return Err(Reject::Unresolved);
    }
    if outcome.state.translate() != protein {
        return Err(Reject::Unsound);
    }
    let (cai, deviation, gc) = (outcome.cai, outcome.deviation, outcome.gc);
    if !cai.is_finite() || !deviation.is_finite() || !gc.is_finite() {
        return Err(Reject::Unsound);
    }
    if deviation > exploration.dev_limit + LIMIT_TOLERANCE
        || gc_distance(gc, exploration.band) > exploration.gc_limit + LIMIT_TOLERANCE
    {
        return Err(Reject::OverLimit);
    }
    Ok(Candidate {
        dna: outcome.state.dna(),
        cai,
        deviation,
        gc,
        run,
        composition,
        restart: outcome.restart,
        seed,
        max_relax,
        dev_limit: exploration.dev_limit,
        moves: outcome.moves,
        substitutions: outcome.substitutions,
    })
}

/// Sort by descending CAI, then ascending deviation, then the sequence itself, and
/// drop repeats of the same gene. Every key is a number or the DNA, so the order is
/// the same whatever order the arrangements finished in.
fn rank(pool: &mut Vec<Candidate>) {
    pool.sort_by(|a, b| {
        b.cai
            .partial_cmp(&a.cai)
            .unwrap()
            .then(a.deviation.partial_cmp(&b.deviation).unwrap())
            .then(a.dna.cmp(&b.dna))
            .then(a.run.cmp(&b.run))
            .then(a.composition.cmp(&b.composition))
            .then(a.restart.cmp(&b.restart))
    });
    pool.dedup_by(|a, b| a.dna == b.dna);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::aa_index;
    use crate::testdata::{ECOLI_TABLE, GFP};
    use std::collections::HashSet;
    use std::path::Path;

    fn host() -> HostTable {
        HostTable::load(Path::new(ECOLI_TABLE), 0.10).unwrap()
    }

    fn protein(seq: &str) -> Vec<AaIdx> {
        seq.bytes().map(|b| aa_index(b).unwrap()).collect()
    }

    /// Long enough to have real synonymous freedom, short enough for a debug build.
    fn short() -> Vec<AaIdx> {
        protein(&GFP[..60])
    }

    fn options() -> BatchOptions {
        BatchOptions::default()
    }

    fn candidate(dna: &str, cai: f64, deviation: f64, run: u8) -> Candidate {
        Candidate {
            dna: dna.to_string(),
            cai,
            deviation,
            gc: 0.5,
            run,
            composition: 0,
            restart: 0,
            seed: 0,
            max_relax: 0.1,
            dev_limit: 0.1,
            moves: 0,
            substitutions: 0,
        }
    }

    #[test]
    fn a_full_set_is_distinct_valid_and_ranked_by_cai() {
        let host = host();
        let protein = short();
        let report = run_record(0, &protein, &host, &options());
        assert_eq!(report.status, Status::Ok, "{report:?}");
        assert_eq!(report.candidates.len(), 5);
        let unique: HashSet<&str> = report.candidates.iter().map(|c| c.dna.as_str()).collect();
        assert_eq!(unique.len(), 5, "candidates must be distinct sequences");
        for pair in report.candidates.windows(2) {
            assert!(pair[0].cai >= pair[1].cai, "ranked by descending CAI");
        }
        for candidate in &report.candidates {
            assert_eq!(candidate.dna.len(), protein.len() * 3);
            assert!(candidate.deviation <= candidate.dev_limit + LIMIT_TOLERANCE);
            assert!(candidate.dev_limit >= options().max_relax);
            assert!(candidate.gc >= 0.30 && candidate.gc <= 0.65);
        }
    }

    #[test]
    fn a_qualifying_first_run_is_not_retried() {
        let host = host();
        let options = BatchOptions { cai_threshold: 0.70, ..options() };
        let report = run_record(0, &short(), &host, &options);
        assert!(report.best_cai().unwrap() > options.cai_threshold);
        assert_eq!(report.attempts.len(), 1, "no retry was needed");
        assert!(report.retry_reason.is_none());
        assert!(!report.low_cai);
    }

    #[test]
    fn an_unreachable_threshold_retries_once_and_still_returns_the_set() {
        let host = host();
        let options = BatchOptions { cai_threshold: 0.999, ..options() };
        let report = run_record(0, &short(), &host, &options);
        assert_eq!(report.attempts.len(), 2, "exactly one retry, never more");
        assert_eq!(report.attempts[0].run, 1);
        assert_eq!(report.attempts[1].run, 2);
        assert_eq!(report.status, Status::LowCai);
        assert!(report.low_cai);
        assert_eq!(report.candidates.len(), 5, "a low-CAI set is still returned");
        assert!(report.retry_reason.unwrap().contains("does not exceed"));
    }

    #[test]
    fn a_threshold_equal_to_the_best_cai_does_not_pass() {
        let host = host();
        let low = BatchOptions { cai_threshold: 0.70, ..options() };
        let first = run_record(0, &short(), &host, &low);
        assert_eq!(first.attempts.len(), 1);
        let best = first.best_cai().unwrap();
        let exact = BatchOptions { cai_threshold: best, ..low };
        let report = run_record(0, &short(), &host, &exact);
        assert_eq!(report.attempts.len(), 2, "equality at the threshold still retries");
        assert_eq!(
            report.attempts[0].best_cai.unwrap(),
            best,
            "the first run matched the threshold exactly and still did not stand"
        );
    }

    #[test]
    fn a_retry_never_loses_the_better_first_run_candidates() {
        let host = host();
        let plain = run_record(0, &short(), &host, &options());
        let forced = run_record(0, &short(), &host, &BatchOptions {
            cai_threshold: 0.999,
            ..options()
        });
        assert!(
            forced.best_cai().unwrap() >= plain.best_cai().unwrap() - LIMIT_TOLERANCE,
            "pooling both runs cannot make the best candidate worse"
        );
    }

    #[test]
    fn a_met_trp_protein_reports_the_shortage_rather_than_padding() {
        let host = host();
        let report = run_record(0, &protein("MWMWM"), &host, &options());
        assert_eq!(report.status, Status::InsufficientCandidates);
        assert_eq!(report.candidates.len(), 1, "MWMWM has exactly one encoding");
        assert!(report.status.is_failure(), "a shortage must not exit zero");
        let unique: HashSet<&str> = report.candidates.iter().map(|c| c.dna.as_str()).collect();
        assert_eq!(unique.len(), report.candidates.len(), "never duplicate to pad");
    }

    #[test]
    fn a_frozen_composition_stays_frozen_through_the_retry() {
        let host = host();
        let frozen = BatchOptions {
            max_relax: 0.0,
            cai_threshold: 0.999, // force the retry
            retry: RetryPolicy::doubling(0.0, 4, 2000),
            ..options()
        };
        assert_eq!(frozen.retry.max_relax, 0.0, "the default retry keeps zero frozen");
        let report = run_record(0, &short(), &host, &frozen);
        assert_eq!(report.attempts.len(), 2);
        for attempt in &report.attempts {
            assert_eq!(attempt.compositions, 1, "one frozen multiset per run");
            assert_eq!(attempt.expansions, 0, "a frozen search expands nothing");
            assert_eq!(attempt.dev_limit, report.attempts[0].dev_limit);
        }
        assert!(
            report.attempts[1].restarts > report.attempts[0].restarts,
            "the retry still buys more arrangement effort"
        );
        let deviations: HashSet<String> =
            report.candidates.iter().map(|c| format!("{:.12}", c.deviation)).collect();
        assert!(deviations.len() <= 1, "every candidate shares the frozen deviation");
    }

    #[test]
    fn an_explicit_positive_retry_limit_unfreezes_the_composition() {
        let host = host();
        let options = BatchOptions {
            max_relax: 0.0,
            cai_threshold: 0.999,
            retry: RetryPolicy { max_relax: 0.15, restarts: 8, expansions: 4000 },
            ..options()
        };
        let report = run_record(0, &short(), &host, &options);
        assert_eq!(report.attempts.len(), 2);
        assert_eq!(report.attempts[0].compositions, 1, "the first run stays frozen");
        assert!(
            report.attempts[1].compositions > 1,
            "an explicit positive retry limit authorises other multisets"
        );
    }

    #[test]
    fn the_same_seed_reproduces_the_same_candidates() {
        let host = host();
        let a = run_record(3, &short(), &host, &options());
        let b = run_record(3, &short(), &host, &options());
        let dna = |r: &RecordReport| r.candidates.iter().map(|c| c.dna.clone()).collect::<Vec<_>>();
        assert_eq!(dna(&a), dna(&b));
        let other = BatchOptions {
            arrange: Options { seed: 99, ..options().arrange },
            ..options()
        };
        assert_ne!(dna(&a), dna(&run_record(3, &short(), &host, &other)));
    }

    #[test]
    fn streams_differ_by_record_run_and_composition() {
        let seeds: HashSet<u64> = [
            stream_seed(1, 0, 1, 0),
            stream_seed(1, 1, 1, 0),
            stream_seed(1, 0, 2, 0),
            stream_seed(1, 0, 1, 1),
            stream_seed(2, 0, 1, 0),
        ]
        .into_iter()
        .collect();
        assert_eq!(seeds.len(), 5, "each index must move the stream");
        assert_eq!(stream_seed(7, 2, 2, 3), stream_seed(7, 2, 2, 3), "and be reproducible");
    }

    #[test]
    fn ranking_sorts_by_cai_then_deviation_then_sequence_and_drops_repeats() {
        let mut pool = vec![
            candidate("TTT", 0.5, 0.01, 1),
            candidate("AAA", 0.9, 0.05, 1),
            candidate("CCC", 0.9, 0.01, 2),
            candidate("GGG", 0.9, 0.01, 1),
            candidate("AAA", 0.9, 0.05, 2), // same gene from the retry
        ];
        rank(&mut pool);
        let order: Vec<&str> = pool.iter().map(|c| c.dna.as_str()).collect();
        assert_eq!(order, ["CCC", "GGG", "AAA", "TTT"]);
        assert_eq!(pool[2].run, 1, "the earlier run keeps the provenance of a repeat");
    }

    #[test]
    fn a_status_is_a_failure_only_when_something_is_actually_missing() {
        assert!(!Status::Ok.is_failure());
        assert!(!Status::LowCai.is_failure(), "a complete low-CAI set is a fallback");
        assert!(Status::InsufficientCandidates.is_failure());
        assert!(Status::NoValidCandidates.is_failure());
        assert!(Status::Error.is_failure());
    }
}
