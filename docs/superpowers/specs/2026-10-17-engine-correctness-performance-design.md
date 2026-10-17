# Rust engine: correctness, bounded memory, and search performance

Date: 2026-10-17

Status: partly implemented; see the status note below

Reviewed baseline: `06061f40d4822b35795743900b5438235abbf161`

> **Status.** The two correctness failures in section 2.1 were fixed in commit
> `3402d01`, together with the quadratic report named in section 2.2. Each now has a
> regression test in `tests/engine.rs`: `g_rich_stems_keep_their_mixed_wobble_matches`,
> `a_positive_deviation_limit_is_never_exceeded` and
> `the_repeat_report_stays_linear_in_the_number_of_copies`. Sections 2.1 and 2.2 are
> kept as the record of what was wrong and how it was reproduced, so the reproducers
> stay readable next to the tests. Everything else here — occurrence lists that still
> grow with multiplicity, the separate scheduler, coordinated repair, and the phased
> delivery in sections 3 to 6 — remains a proposal.

## 1. Problem and decision

The composition-first Rust rewrite is the right foundation: synonymous swaps preserve
codon counts, deviation and CAI, while incremental constraint scores guide arrangement.
However, several operations still grow with occurrence multiplicity or sequence length,
and two correctness failures undermine comparisons at tight deviation limits.

Retain this architecture and implement two complementary algorithms:

1. An exact incremental global search with compact violation scheduling, aggregated
   hairpin counts, and efficient repeat and partner indices.
2. Bounded neighborhood repair that jointly changes several codons when individual
   swaps cannot improve the sequence, preserving quotas whenever possible.

Fix the validity contract before making performance claims. Optimize time to a validated
solution, success rate within a budget, and peak memory together. Fast termination with
unresolved constraints is not successful optimization.

This proposal supersedes the unconditional O(1)-per-move and O(N)-successful-search
claims in the original design. Largest-remainder allocation optimizes its composition
objective; it does not solve the joint sequence-constrained problem exactly. Whole-gene
GC is invariant under swaps, but must also be checked when substitutions change counts.

Non-goals: a wholesale rewrite, a new biological model of RNA folding, an exponential
whole-gene state table, or a promise that every feasible input is solved in linear time.

## 2. Evidence and corrections to the existing design

### 2.1 Reproduced correctness failures

**Positive deviation limits can be exceeded before search.**
`composition::build` clamps `budget_left` to zero when the starting deviation exceeds
`max_relax`, rather than rejecting or repairing that composition. For GFP with the
repository E. coli table and host threshold 0.10, a requested 0.001 limit returns
deviation 0.008333333333333333. A CLI run with the repeat-only experimental patch exited
successfully with deviation 0.0083 and zero reported violations. The patch does not
change composition or validity logic; the unmodified baseline composition probe also
reproduces the excess.

**Mixed wobble matches can be missed.** `Hairpins::variants_of` retains only the all-C
and all-T alternatives when more than 64 variants exist. Contrary to the comment,
rebuild does not perform another exact scan. For DNA `GGGGGGGAAATTTCTCCCCCAC`, the stem
`GGGGGGGAAA` accepts target `TTTCTCCCCC` under the detector's own directional G-to-C/T
rule. The baseline reports zero hairpin matches; a brute-force oracle counts one.

Both failures above are fixed as of commit `3402d01`; the suite that missed them now
covers them directly. The wobble cap is gone, so all 2^g variants of a stem are
enumerated, bounded at 2^10 because a stem is ten bases. The whole-gene GC correction
now consults `max_relax` before each move and reports `gc_in_band: false` rather than
spending past the limit.

`scripts/validate.py` still independently recomputes translation, deviation and CAI
only: it does not check all sequence constraints, so it remains an incomplete validity
oracle, which is the gap section 3 addresses.

### 2.2 Complexity and search limitations found in the code

| Current operation | Consequence | Design correction |
|---|---|---|
| Repeat delta copies, filters and sorts occurrence buckets | Proposal work grows with bucket size | Derive the contribution from count and extrema |
| Accepted repeat edits insert/remove in sorted vectors | Potential O(N) element movement | Keep compact vectors initially; measure large buckets before adopting ordered indices |
| Hairpin delta visits every stem linked to a changed target | Work grows with repeated stem multiplicity | Aggregate source and target sequence classes |
| Hairpin report creates a record per matching pair | O(N²) records on some inputs | One compact witness per source/anchor, with multiplicity |
| Repeat report duplicates its full position list per occurrence | At least quadratic diagnostic bytes on repetitive inputs | Shared/grouped bounded summaries |
| Homopolymer scoring expands to maximal run boundaries | A long run can cause O(N) rescanning | Ordered run boundaries with O(log N) edits |
| Queue refresh invokes formatted, sorted reports | Allocation and formatting in the search loop | Separate scheduling from presentation |
| Block invalidation scans a map and uses a distant swap's bounding interval | Revisits unrelated regions; misses pattern dependencies as a model | Two actual spans plus pattern versions |
| A capped sweep repeatedly starts at the first positions | Wasted/no-op proposals and incomplete coverage | Alternative-codon pools and rotating or seeded sampling |
| Failed spans can be permanently blocked | Heuristic exhaustion is described as unfixability | Versioned retry status and explicit termination reasons |
| Only the best score is saved | Final sequence can be worse than an earlier state | Keep and return the best validated sequence |
| Every restart outcome and batch record is retained | Working memory can grow with restarts and batch size | Deterministic reduction and bounded batch streaming |

### 2.3 Measured experiment, and what it does not establish

A minimal repeat-delta prototype preserves the original score, but uses only the new
count, minimum and maximum occurrence positions. It passed all 80 tests and 20,000
additional baseline/prototype move comparisons. The persistent index is unchanged.

Mean time per proposal across 10,000 evaluations, on repeated `GCT` codons with a
synonymous `GCT -> GCC` edit:

| Length | Baseline | Prototype | Ratio |
|---:|---:|---:|---:|
| 100 aa | 4.418 µs | 1.313 µs | 3.36× |
| 400 aa | 5.710 µs | 1.052 µs | 5.43× |
| 1,600 aa | 15.797 µs | 1.087 µs | 14.53× |
| 6,400 aa | 55.937 µs | 1.052 µs | 53.16× |

These are kernel timings, excluding accepted-move application. Complete-run timings
did not establish a meaningful speedup. In the second comparison round, mean times
across seeds 1/2/3 were 9.741 vs 9.697 ms for GFP, 297.000 vs 296.882 ms for uniform
random 1,600-aa proteins, and 2,824.680 vs 2,840.623 ms for 6,400-aa proteins. Apparent
3–4% improvements on larger inputs in the first round disappeared on repetition.
DNA hashes, scores, deviation and move counts matched in all paired cases.

These runs used one thread, one restart and `max_relax=0.01`. GFP seed 1 reached zero
reported violations; seeds 2 and 3 each retained one. Random 1,600-aa runs retained
16–21 violations and 6,400-aa runs retained 229–300. These timings measure termination,
not validated time to success. The sample is too small to estimate a success rate;
the synthetic cases may contain impossible constraints.

On alternating-AT DNA corresponding to 100/200/400 aa, hairpin reports contained
4,380/17,760/71,520 records. Repeat detail strings alone occupied
204,400/858,400/3,631,428 bytes. Those observations motivate removing expanded reports
from the search state before pursuing small allocation optimizations elsewhere.

Environment, raw results, source probe and experimental patch are in the
[evidence bundle](../../benchmarks/2026-10-17-engine-review/README.md). The patch is
an experiment, not an implementation shipped by this design PR.

## 3. Validity contract and compatibility

### 3.1 Hard requirements

A successful result must preserve translation, use allowed synonymous codons, satisfy
the requested absolute deviation limit and whole-gene GC band, and pass every enabled
sequence constraint. CAI is a ranking objective after feasibility, not a substitute
for it. Check this contract on initialization, on committed count-changing moves, and
again with a fresh independent validation pass before reporting success.

Treat a multi-codon move atomically: intermediate assignments may be explored privately,
but no intermediate state is emitted as a solution. Transaction rollback must restore
sequence, counts, indices, scores, pool membership, and scheduler versions.

### 3.2 Deviation and early lower bounds

Keep the current rounded-target metric explicit. For residue count n_a and normalized
usable frequency f_c, define r_c = round_half_even(f_c * n_a), R = sum_c r_c, and:

```text
L(counts) = sum_c abs(count[c] - r_c)
D(counts) = L(counts) / (2 * R), when R > 0
```

Precompute r_c and R once because synonymous edits preserve n_a. Update L from only the
affected codon counts. Convert a positive limit to an admissible integer numerator
bound once, with explicit boundary tests for exact equality and adjacent representable
values; do not introduce an undocumented slack epsilon. Reject negative/non-finite
limits. The existing R=0 behavior needs a named compatibility test: preserve its
reported D=0 convention, but expose that the rounded metric is uninformative in this
degenerate case rather than dividing by zero or claiming an informative lower bound.

Ignoring GC and arrangement, the minimum numerator is:

```text
L_min = sum_a abs(n_a - sum_{c encoding a} r_c)
```

Each group's counts must sum to n_a, so the triangle inequality gives this bound;
adding or removing counts from the rounded targets achieves it when usable synonyms
exist. If R > 0 and L_min exceeds the permitted numerator, report proven infeasibility
for the profile limit before searching. A bound below the limit does not prove joint
feasibility. Failure of one GC-adjustment heuristic or one composition allocation also
does not prove infeasibility.

Derive attainable whole-gene GC bounds from the allowed codons and residue counts,
including integer GC-count limits. Bound composition-adjustment work and detect repeated
count states; do not allow a greedy GC loop to oscillate indefinitely. Recheck global GC
after substitutions, since local 20/50/100-base bands do not enforce the global band.

### 3.3 CLI migration and result status

For positive `--max-relax`, enforce an absolute final limit. Preserve the documented
`--max-relax 0` freeze behavior during a compatibility period, introduce an explicit
`--freeze-composition`, and warn that the zero sentinel is not an absolute zero target.
A future major interface can remove the sentinel. In explicit freeze mode, freeze the
selected multiset; if a separate positive absolute limit is supplied, it still applies.
The library should represent these as separate policy fields, not infer freeze mode
from a nonpositive `budget_left`. Zero remaining headroom must not prohibit an atomic
change that leaves the deviation unchanged or reduces it.

Introduce structured termination states, with names finalized in implementation:

- `Feasible`: all requested invariants independently validated.
- `ProvenInfeasible`: accompanied by a checkable bound or complete proof for the full
  declared problem, not just a failed neighborhood.
- `BudgetExhausted`: best candidate and residual counts are available; feasibility is unknown.
- `InvalidInput`: unsupported/invalid configuration or input.

Failed candidates may still be written for diagnosis, matching the current CLI convention,
but must be marked incomplete and yield a nonzero process exit. A bounded neighborhood
can be infeasible with its outside positions frozen without proving the whole gene
infeasible. A failing independent validator is an internal error, never success.

### 3.4 Freeze the constraint semantics

Use a versioned rule set and record it in benchmark output. Preserve the Rust baseline's
full 10-base stems at offsets 0,10,20,…, targets at every complete 10-mer start,
directional G-to-C/T wobble rule, overlap/self-match behavior, and exclusion of partial
tail stems. Removing the 64-variant truncation is a correctness fix within that rule.
Changes to overlap, stem spacing, partial tails, symmetric wobble or thermodynamic
folding require a separate model change and separate comparison.

Likewise record the chosen start-site rule, all GC windows, homopolymer threshold,
repeat-gap rule, restriction motifs and splice settings. Python and Rust defaults
differ; speed comparisons must use equivalent rules or explicitly identify differences.

## 4. Algorithm A: exact incremental global search

### 4.1 Data layout and edit transactions

Retain byte-encoded codons/bases and fixed 64-codon counts. Replace heap-allocated
one/two-position patch buffers with fixed-capacity buffers whose capacities are derived
from the move arity. Deduplicate overlapping affected ranges and k-mer starts. Extend
the patch abstraction to bounded multi-codon moves before using cycles or neighborhood
repair; a capacity derived for two codons is not safe for larger moves.

Keep `delta` side-effect free. An evaluated move may carry a bounded edit summary tied
to the state's generation, allowing `apply` to reuse work without accepting a stale
evaluation. Check generation before committing. Use u64 totals and checked signed
delta arithmetic for pair-counting constraints; define a supported sequence-length
limit that prevents overflow in all index and score types.

### 4.2 Repeat scoring and accepted updates

For a packed 9-mer bucket with sorted positions P:

```text
contribution(P) = len(P), if max(P) - min(P) >= 9
                  0, otherwise
```

A codon edit changes at most 11 starts; a two-codon move changes at most 22. For each
affected key, compute the new length from departures/arrivals, find the surviving
extrema by skipping departing endpoints, then include arriving positions. The number
of skipped endpoints is bounded by the touched starts, regardless of bucket length.
No bucket-sized copy or sort is needed for proposal scoring.

The experimental patch demonstrates this change. It does not fix accepted vector
insert/remove costs. Measure those separately. If they dominate, migrate large buckets
to ordered sets or use a global ordered `(key, position)` index with per-key metadata.
Either provides O(log N) accepted updates and O(N) storage, but may lose cache locality
on normal short buckets. Select the representation using evidence, not a blanket claim
that trees are faster. Guard sequences shorter than nine bases before generating starts.

### 4.3 Exact grouped hairpin counts

Let S[s] count aligned source stems with packed pattern s, T[t] count all target 10-mers,
and A(s,t) be the frozen compatibility predicate. Then:

```text
H = S^T A T
delta_H = delta_S^T A T + S^T A delta_T + delta_S^T A delta_T
```

Implement S and T as sparse count maps, plus source/target position pools for selecting
witnesses. Never materialize A or all matching occurrence pairs. Enumerate compatible
packed keys on demand. A ten-base pattern has at most 2^10 = 1,024 compatible keys,
independent of N; typical patterns have fewer. The inverse query is not the same
directional mapping: target T admits source A or G at the reversed position, while
target C admits G only. Test forward/inverse enumeration separately.

For each proposal, merge all changed source/target keys, evaluate the two linear terms
against the old maps, and add the cross term. Omitting the cross term double-counts or
misses interactions when both endpoints change. Apply the map shifts atomically and
reuse the evaluated total delta. At fixed stem length, scoring depends on the bounded
number of changed keys and compatible patterns, not their occurrence multiplicities.

Do not retain a compatibility-edge list per occurrence. Optional pattern caches need
explicit capacity limits and eviction; remove empty map entries and recycle position
storage. This keeps persistent storage O(N), instead of O(N²) pairs or an unbounded
cache accumulated over all historical states.

### 4.4 Local constraints and homopolymers

Keep direct rolling counts for the fixed 20/50/100-base GC windows. Their edit work is
bounded by configured widths, not O(1) independently of arbitrary user-selected widths.
Use fixed-pattern masks/automata for restriction, splice and start-site constraints.

For homopolymers, preserve the existing one-violation-per-maximal-run score using
ordered run boundaries. A bounded base edit splits/merges only a bounded number of
runs and can update the score in O(log N), with O(N) memory. Counting every forbidden
(limit+1)-base window would also be local but changes the penalty landscape; do not
silently substitute that score in a semantics-preserving optimization.

### 4.5 Scheduler separate from reporting

Replace `report() -> Vec<Violation>` as a scheduling interface with compact records:

```text
Witness { family, anchor_or_pattern, dependency_version, multiplicity }
```

Store at most O(N) active records and bounded pending edits. Format diagnostic strings
only when serializing a final report, with an explicit example/byte cap and exact
aggregate counts. A full expanded pair list is optional output whose output-sensitive
cost must not be charged to the normal search-memory guarantee.

Invalidate local records using the actual disjoint edited spans. For repeats and
hairpins, maintain generations for affected pattern counts and revalidate a selected
witness lazily. Remote changes can create or clear a conflict without touching its
anchor, so position overlap alone is insufficient. Include constraint family in keys.

Discover new witnesses from touched local anchors and changed patterns. A changed target
may affect many source occurrences: schedule the source class, not every occurrence.
Use bounded periodic complete witness scans as a safety net; charge them to a scan
budget. Do not trigger a full refresh after every accepted move. If weighted selection
is desired, represent multiplicity as a weight rather than duplicate queue entries;
unweighted class selection is a deliberate search-policy change and must be benchmarked.

A failed capped repair records a temporary failure at a dependency generation. Retry
after relevant changes or a bounded diversification interval. A cap does not establish
that every possible move was examined, and must not permanently classify the span as
unfixable.

### 4.6 Useful partners, budgets, and best state

Maintain pools indexed by amino acid and codon, with inverse slots for O(1) membership
updates. Select an alternative codon class and then a position, avoiding identical-codon
swaps. Rotate capped exhaustive candidates or use a reproducible permutation so retries
cover different positions.

Charge every evaluated proposal and fallback node to explicit budgets, including
rejections and nested paired repairs. Accepted-move counts alone do not bound work.
Preserve seeded deterministic ordering; never depend on HashMap iteration order.

Maintain best codons and a dirty-position list/bitmap. At improvement, copy only positions
changed since the previous best, then clear the list. This adds O(N) space and amortizes
best-buffer writes over edits rather than copying the entire gene on every improvement.
At termination rebuild and independently validate that best state. Compare hard
feasibility first, followed by documented score, deviation, CAI and deterministic ties.

## 5. Algorithm B: bounded coordinated repair

Trigger repair when local search stalls at a dependency-stable conflict. Escalate from
three-codon cycles and coordinated disjoint swaps to a small neighborhood containing
conflict positions and selected remote same-amino-acid donor positions. Initial tuning
candidates are 12–24 variable positions; these values are not benchmarked optima.

Use depth-first search with an undo trail, variable ordering and constraint propagation.
For count-preserving repair, the selected positions must use exactly their original
codon multiset. For budgeted repair, include fixed outside counts in the global
deviation and GC checks. When quotas are exhausted, prune; also prune when the remaining
positions cannot supply a required quota or attain GC bounds.

Use partial assignments only for sound pruning: a fully assigned forbidden motif can
be checked immediately; a motif containing unassigned bases cannot be rejected merely
because its old sequence was forbidden. Include boundary windows and nonlocal
repeat/hairpin interactions with the fixed outside sequence. Once a full assignment is
available, evaluate its atomic global delta and validate before committing.

Bound neighborhood size, nodes, time and retained scratch data. Roll back completely on
failure or timeout. With at most six codons per position, worst-case time is exponential
in neighborhood size B; depth-first traversal and bounded trails avoid retaining 6^B
states. Account for global scoring work per visited node. Do not store an unrestricted
beam or all learned conflicts if claiming bounded scratch memory.

If one multiset repeatedly stalls, explore a bounded number of alternative allocations
within the absolute deviation and GC constraints. These can include equal-deviation
allocations even when no extra headroom remains. Preserve an explicit freeze request.
Failure to arrange one allocation is not a certificate that all allocations fail.

CP-SAT is an optional alternative for the bounded neighborhood, not a required new
dependency or a proven faster whole-gene replacement. A solver's learned clauses need
memory limits, and its timeout/unknown result must remain distinct from infeasibility.
The initial implementation should be native Rust to retain predictable integration
with the same incremental indices and undo machinery.

## 6. Time and space targets

Let N be nucleotide length, M evaluated fixed-size proposals, k the fixed motif length,
F complete witness/validation scans, W active workers and B the repair-neighborhood cap.
The number of sequence records and configured restarts are separate quantities.

After the proposed representations and bounds are implemented, the target is:

```text
time = O(N * C_init(k) + M * (C_local(k) + log N)
         + F * N * C_scan(k) + sum(repair_node_costs))
working memory = O(W * N) + bounded per-worker repair/cache workspace
```

Hash-map lookup bounds are expected, not adversarial worst-case guarantees. The log N
term allows ordered repeat/run indices and scheduling structures. C_local includes
bounded compatibility enumeration; exponential dependence on configurable k must not
be hidden if stem lengths become user-configurable. Keep F bounded or amortized, and
include initialization, validation, reporting, and failed attempts in measurements.

M proportional to N gives a near-linear or N log N work allowance under fixed rule sizes
and bounded scans/repair budgets. It does not guarantee a solution within that allowance.
The repeat-only experiment still has O(N) accepted updates and does not meet this full
target. Current expanded reports can also exceed the proposed memory bound.

Share immutable host/rule data. Limit active searches to W and reduce completed restarts
immediately using a deterministic ordering instead of retaining all outcomes. Stream
batch records with bounded in-flight and reorder buffers. A fixed proposal budget gives
reproducibility across thread schedules; wall-clock cutoffs may change how far a run
progresses and must be reported as such. Avoid claiming wall-clock and work-budget
reproducibility are identical.

## 7. Validation and benchmark acceptance

### 7.1 Correctness gates

The final reference validator must independently implement the frozen predicates rather
than call the optimized matcher or its variant generator. A full rebuild of the same
incorrect matcher is not independent evidence. Use a brute-force oracle on short inputs;
for long benchmarks, use an independently implemented exact index/scan validator and
cross-check it against the brute-force oracle on short cases.

Required cases:

- GFP below its rounded-profile lower bound; limits exactly at attainable values;
  non-finite/negative limits; R=0; zero-sentinel migration and explicit freeze behavior.
- Initial GC adjustment and final substitutions against global GC, including empty
  integer GC bands and a bounded failure to find a feasible composition.
- Seven through ten Gs per stem, mixed C/T targets, forward/inverse wobble direction,
  repeated stems, self/overlapping targets, partial tails and sequences shorter than k.
- Overlapping edits, same-position/no-op swaps, all bucket entries removed, creation
  and deletion of extrema, and high-multiplicity repeat buckets.
- Score/delta consistency after long random edit streams; u64 pair totals and supported
  length limits; complete rollback after a capped multi-position search.
- Remote pattern invalidation, distinct constraint families at one span, repeated
  unsuccessful capped sweeps, and discovering newly created conflicts.
- Best-state retention after a worsening move; deterministic reduction across worker
  counts with a fixed work budget; bounded queue/report/cache sizes.

### 7.2 Fair performance comparison

Maintain two comparisons: a semantics-preserving kernel comparison against the pinned
baseline, and whole-engine comparison under the corrected, independently validated
contract. Fixing missed hairpins may initially cost time and lower apparent success;
do not count previously invalid outputs as successful reference solutions.

Use real protein inputs, weighted and uniform synthetic inputs, repetitive/G-rich
stress cases, and known-infeasible cases. Include a known-feasible corpus with saved
witness DNA meeting the fixed host/rules/limit, so failure to solve it is distinguishable
from impossible requirements. Compare multiple lengths and host tables, and positive
limits at feasible margins above their lower bounds; test freeze mode separately.

Record compiler, hardware, commit, rule version, sequence/host hashes, seed, threads,
restarts, proposals/nodes, scans, elapsed time, peak RSS, validity, residual counts,
deviation and termination reason. Use interleaved baseline/candidate repetitions after
warmup, at least 30 search seeds for success-rate evaluation, and independent validation
of every alleged success. Report uncertainty, deadline failures and censored time to
success, rather than averaging only successful runs or comparing failed-run durations.

Acceptance gates for a production implementation:

1. All correctness regressions pass and every reported success passes the independent validator.
2. Stress-index and scheduling memory scales approximately linearly at fixed workers;
   no default diagnostic expansion proportional to matching occurrence pairs.
3. Tight-limit success at equal budgets is non-inferior on the agreed feasible corpus,
   with a demonstrated improvement on the targeted stalled cases before claiming it.
4. An end-to-end speed claim requires repeated time-to-valid-solution improvement,
   beyond run-to-run uncertainty. A faster microbenchmark alone does not meet this gate.

## 8. Delivery sequence

| Phase | Modules | Deliverable and exit criterion |
|---|---|---|
| P0: validity | `composition`, `search`, `main`, reference validator | Absolute positive limit, global GC gates, typed outcomes, independent adversarial tests; exact full wobble enumeration as the initial correctness fix |
| P1: compact state | constraints interface, scheduler, reporting | No formatted strings in proposals/queue refresh; compact witnesses and bounded diagnostics; memory stress tests |
| P2: exact indices | `repeat`, `hairpin`, `homopolymer`, `state` | Extrema delta, grouped hairpin cross term, ordered runs, bounded edit buffers; independent delta/oracle checks |
| P3: search quality | `search`, position pools | Useful partner sampling, dependency retries, charged proposal budgets, best-state retention |
| P4: repair | bounded patch/undo API, neighborhood solver, composition alternatives | Atomic quota-preserving cycles and bounded DFS with complete rollback and globally validated commits |
| P5: integration | restart reduction, batch I/O, benchmark tooling | Bounded W-worker memory, deterministic work-budget runs, interleaved independent quality/performance evaluation |

Land each phase separately with its regression tests and measurements. P0 can increase
cost by detecting more matches; compare P1–P5 against that corrected reference as well
as recording historical results. Keep the Rust 1.73 toolchain and pinned dependencies
unless a separate compatibility change is justified. The experimental repeat patch
can be considered during P2; it must not silently become the baseline for its own
microbenchmark.

## 9. Risks and alternatives

- Fixed composition can be unarrangeable even when another permissible composition is
  feasible. Bounded alternatives improve coverage but do not prove completeness.
- Grouping and weighting witnesses changes search trajectories. Preserve scores exactly,
  version the selection policy, and measure success as well as throughput.
- Enumerating 1,024 partners is exact and bounded at k=10, but may cost more on ordinary
  inputs than the old truncated path. Use profiles to decide on bounded caching or tries.
- Ordered sets improve worst-case updates but have allocation/cache overhead. Keep the
  compact-vector option until accepted-update measurements justify migration.
- Depth-first repair limits memory, not combinatorial difficulty. Report exhausted
  budgets honestly; a whole-gene CP model remains an optional separate investigation.
- Biological validation of a designed gene is beyond these sequence predicates. This
  proposal improves compliance with declared computational requirements, not a guarantee
  of expression or manufacturability.
