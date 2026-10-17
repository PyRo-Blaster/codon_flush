# Rust application: five ranked candidates with one CAI retry

Status: proposed; no engine changes made. Inspected revision: `fb14f64`.

## What already exists

| Requirement | Current Rust implementation |
| --- | --- |
| Multi-entry protein FASTA | `src/fasta.rs` parses records; `src/main.rs` processes them in parallel and restores input order. |
| Python-compatible codon usage table | `--host-table` accepts nested amino-acid/codon JSON from `python_codon_tables`, as used by the Python implementation, and flat codon JSON. DNA and RNA codon spellings are supported. |
| Optimized DNA FASTA | `--output` writes one sequence per protein, preserves headers, and checks translation. |
| Five candidates ranked by CAI | Missing. The engine currently returns one selected outcome; restarts can share a composition and produce identical DNA or CAI. |
| Retry and low-CAI fallback | Missing. CAI is reported, but does not drive a retry policy. |

Extend the existing command-line application. Input remains protein FASTA; a codon usage table supplies frequencies, not an alternative genetic code. DNA-input mode and a web UI are outside this increment.

## Agreed behavior

For **each input FASTA record independently**:

1. Run optimization to obtain five distinct, valid DNA candidates, exploring multiple codon compositions and arrangements.
2. Rank candidates by independently recomputed, full-precision CAI, highest first.
3. If at least one candidate has **CAI > 0.80**, return the five ranked candidates without retrying. The other four need not exceed 0.80.
4. If none exceeds 0.80, change the allowed parameters and run optimization once more. There are at most two runs per protein.
5. Pool valid, distinct candidates from both runs and return the best five by CAI. If none exceeds 0.80 after the second run, still return the five candidates and emit a clear low-CAI warning.

The threshold is a retry trigger and quality warning, not a hard filter on output. Equality at 0.80 does **not** pass. This replaces the earlier proposal for a two-sided CAI target and tolerance; do not implement `--target-cai` or `--cai-tolerance` for this workflow.

“Valid” means exact translation, allowed codons, applicable deviation/whole-gene GC limits, and no unresolved sequence constraints. The low-CAI fallback relaxes only the CAI expectation. It must not hide translation errors or unresolved synthesis/expression constraints.

Five distinct sequences cannot always exist, for example for a Met/Trp-only protein. Never duplicate a sequence to pad output. Write the valid candidates available, report the actual count and `insufficient_candidates`, and exit nonzero. Likewise, a bounded search that finds fewer than five candidates must report the shortage rather than claim there are no other possible encodings. This exception is separate from an ordinary five-candidate low-CAI warning.

## Proposed command and retry settings

New workflow example (not yet implemented):

```sh
cargo run --release -- \
  --input proteins.fasta --output genes.fasta \
  --host-table tests/data/e_coli_316407.json \
  --candidates 5 --cai-threshold 0.80 \
  --max-relax 0.10 --retry-max-relax 0.20 \
  --restarts 4 --retry-restarts 8 \
  --report-json report.json
```

CLI defaults for the new workflow: five candidates and threshold 0.80. Keep the existing single-outcome library API as a compatibility wrapper; document the CLI's change from one output sequence to five per input protein.

| Parameter | First run | Retry |
| --- | --- | --- |
| Candidate output count | 5 | Best 5 from pooled results |
| CAI threshold | 0.80 | 0.80, unchanged |
| Maximum profile deviation | `--max-relax`, default 0.10 | `--retry-max-relax`, default `min(1.0, initial + 0.10)` for a positive initial limit; normally 0.20 |
| Arrangement restarts per composition | `--restarts`, default 4 | `--retry-restarts`, default twice the initial value; normally 8 |
| Composition search budget | Configured initial budget | Twice the initial budget |
| Random stream | Derived from seed, record ordinal, candidate/composition, restart | Separate stream including run index |
| Host table/threshold, genetic code, GC bounds, sequence constraints | User settings | Unchanged |

These retry values are engineering defaults, not biologically validated thresholds or guarantees that CAI will improve. Show the effective retry policy before optimization and record actual settings in the report. A user can set an equal retry deviation limit to spend more search effort without increasing allowed deviation. Require an explicit retry limit to be at least the initial limit. Validate finite numeric settings, threshold in `[0, 1)`, deviation limits in `[0, 1]`, positive candidate/restart/search counts, and checked budget multiplication.

Preserve the special meaning of `--max-relax 0`: freeze the initial apportioned composition. Its default retry deviation limit also remains zero; retry only increases arrangement effort. An explicitly positive `--retry-max-relax` authorizes unfreezing on retry. Report the small integer-rounding deviation that can exist even when composition is frozen. Positive limits are hard bounds; do not enlarge them implicitly beyond the declared retry policy.

## 1. Define shared metrics, requirements, and result types

Files: `src/metrics.rs`, new `src/requirements.rs` and `src/batch.rs` (or equivalent), `src/lib.rs`.

- Define `BatchOptions`, `RetryPolicy`, `Candidate`, `RunReport`, and `RecordReport`. Separate hard validity requirements from the CAI quality threshold.
- Retain the Rust CAI convention: `w_c = f_c / max(f_synonyms)` after normalization/thresholding; `CAI = exp(sum(n_c * ln(w_c)) / N)`. Exclude Met and Trp; the current code includes stops. If no codons count, CAI is 1. Record this convention. Table compatibility with Python does not imply identical legacy Biopython CAI values, whose denominator differs.
- Compute CAI from canonical codon counts, reject nonfinite metrics, and compare the unrounded value strictly with the threshold. Use small documented numerical tolerances for hard GC/deviation comparisons; never round CAI to decide whether retry is required.
- Centralize validity checks for allowed codons, GC, and deviation. Positive `max_relax` must hold even for the initial integer composition; do not inherit an oversized budget from an over-limit initial allocation. Preserve zero's freeze behavior explicitly.
- Give each candidate provenance: input ordinal, run, composition, restart, seed, effective deviation limit, and final metrics. A retry candidate can legitimately exceed the first run's deviation limit; always report which bound applied.

## 2. Explore different codon compositions

Files: `src/composition.rs`, new `src/composition/candidates.rs` or equivalent.

Different arrangements of the same multiset have the same CAI. Simply returning five restarts of the current builder will not meaningfully optimize or rank CAI. Add a bounded composition-candidate generator while preserving `composition::build` for existing callers.

- Use nonnegative integer counts, with synonymous-group totals fixed by the protein. Select only `HostTable::usable_synonyms`.
- Seed exploration with the existing host-profile allocation and allocations biased toward higher relative-adaptiveness weights. Search CAI and whole-gene GC jointly, because changing one can undo progress on the other.
- Implement deterministic bounded beam search over synonymous one-count transfers, with coordinated two-transfer proposals for barriers. Deduplicate count vectors and cap beam width, retained states, and expansion count. Rank hard-constraint distance first, then higher CAI and lower host deviation, with stable count-vector tie-breaking. Intermediate composition-search states may violate limits; only valid compositions enter arrangement.
- Retain multiple feasible compositions, including the valid host-profile baseline when available, rather than just the highest-CAI multiset. Arrangement can fail for one composition and succeed for another.
- In frozen mode, use only the fixed composition, while exploring distinct arrangements.
- Expose/document initial composition-search limits; double the expansion budget on retry and report consumption. Measure coverage/runtime before fixing defaults.
- Use inexpensive independent CAI/GC bounds only when they actually prove a limit. A failed bounded search is `search_exhausted`, not proof that the requested CAI or five encodings are impossible. Even if a CAI upper bound is below the threshold, perform the specified second run to seek better fallback candidates.

This generator is a bounded heuristic; do not claim it finds the global top five encodings.

## 3. Collect valid arrangements and rank candidates

Files: `src/arrange/search.rs`; `src/arrange/state.rs` only if necessary.

- Add an API that collects several outcomes across compositions/restarts instead of immediately reducing to one. Generate a bounded surplus of arrangements so deduplication still has a chance to leave five distinct candidates.
- Preserve hard requirements during substitutions, including whole-gene GC and the run's positive deviation limit. Swaps preserve composition. Apply the same checks to paired substitutions, initially requiring intermediate and final counts to remain valid; document this coverage restriction.
- CAI is a ranking objective, not a hard boundary during repair. Permit a CAI decrease that clears a violation, then score the actual finished sequence. Do not discard a valid low-CAI candidate needed for fallback.
- Recompute translation and final metrics and rebuild sequence constraints before marking a candidate valid. Coordinate budget/validity corrections with the existing engine correctness design without assuming its proposals are implemented.
- Deduplicate by complete DNA sequence per input record, across compositions, restarts, and runs. Candidates with equal CAI remain useful if their DNA differs; no minimum sequence distance is required.
- Sort valid candidates by descending CAI, ascending deviation, then lexicographic DNA for deterministic ties. Keep provenance deterministically when duplicate DNA is encountered. Pool and rerank both runs so a retry cannot replace better first-run results with worse ones.

## 4. Implement per-record retry orchestration

Files: new `src/batch.rs`, `src/main.rs`.

- Run the initial candidate search, validate/deduplicate/rank it, and test whether any valid candidate strictly exceeds the CAI threshold.
- Trigger exactly one retry for records without a qualifying candidate, including a first run yielding no valid candidates. Reuse the input/table and apply the declared retry settings.
- If a first run has a qualifying candidate but fewer than five valid distinct sequences after its search budget, report a shortage; the automatic second run is triggered by low CAI as specified, not by shortage alone.
- Merge valid candidates and retain the best five. If the final best CAI is at most 0.80, emit a warning such as: `protein_x: neither run produced CAI > 0.80; returning 5 ranked candidates (best CAI 0.793421).` Display enough precision to explain boundary decisions and put full precision in JSON.
- Mark a complete set as `ok` or `low_cai`; mark shortages as `insufficient_candidates` and records with none as `no_valid_candidates`. Include separate quality-warning and search-exhaustion fields so multiple conditions remain visible.
- A complete five-candidate set with low CAI is a successful fallback: warn on stderr and in the report, and exit zero unless another record has an actual error/shortage. Input, optimization-validity, shortage, and output errors produce a nonzero exit. Successful peers still produce output.
- Derive random streams from stable numeric indices, not randomized hashes or thread scheduling. Parallelize per-record work while retaining deterministic output ordering.

## 5. Integrate FASTA, reports, and input validation

Files: `src/main.rs`, `src/fasta.rs`, `src/host.rs`.

- Output up to five DNA FASTA records for each input, grouped by input order and ranked within each group. Use collision-free IDs such as `record_000001_candidate_01`; preserve the exact original ID/description in JSON and annotate the FASTA description with source ID, rank, CAI, run, and low-CAI warning when applicable. Input ordinal distinguishes duplicate source IDs.
- Add optional `--report-json FILE`: run-level version/convention/seed/policy, and one entry per input protein with original header, threshold, actual candidate count, attempt settings/outcomes, retry reason, status/warnings, and every returned candidate's metrics/provenance/output ID. Always show low-CAI and shortage warnings even without JSON or with `--quiet`.
- Reuse existing multi-record FASTA parsing and JSON table loading. Reject blank headers and internal stops, accept optional terminal `*`, and report malformed input with record/line context. Treat malformed FASTA/table as an input error before optimization.
- Add nested amino-acid-key/codon consistency and finite normalization-total checks. Retain RNA-to-DNA normalization and normalized duplicate detection. Document the current complete-table requirement, including usable stop codons, and provide a full working Python export fixture.
- Reject input/output/report path collisions before writing. Write files through temporary siblings and rename each after completion; do not claim multiple file renames are an atomic transaction.

## 6. Verify and document

Files: tests in changed modules, `tests/engine.rs`, new `tests/cli.rs`, `tests/data/`, `scripts/validate.py`, `README.md`.

Acceptance coverage:

- Multi-entry protein FASTA plus the existing Python-format table produces five distinct valid DNA sequences per ordinary test protein, correctly grouped/ranked, with exact translation and independently recomputed CAI.
- First-run success skips retry even when only one of five candidates exceeds 0.80; first-run low CAI triggers exactly one retry; retry success returns the pooled top five; two low-CAI runs still output five with a warning and zero exit status.
- CAI exactly 0.80 triggers retry/warning; values on either side are judged at full precision. Nonfinite metrics cannot qualify. Test retry branching with controlled search results as well as real CLI fixtures to avoid stochastic assertions.
- Changed deviation limits, doubled budgets/restarts, distinct random streams, and recorded provenance match the declared policy. Frozen zero remains frozen by default; an explicit positive retry limit permits composition changes.
- Candidate pools deduplicate across both runs, preserve better first-run candidates, sort ties deterministically, and retain distinct DNA with equal CAI. Fewer than five possible/found candidates produce an honest count and nonzero exit, never duplicate padding.
- Hard GC, allowed-codon, deviation, sequence-constraint, and translation requirements remain valid even for low-CAI fallback. Cover single/paired substitutions and rollback.
- Exhaustive tiny-protein oracles verify metrics/validity and measure search coverage. Distinguish bounds proving impossibility from bounded search exhaustion.
- Test Met/Trp-only inputs, short proteins, terminal/internal stops, duplicate/blank headers, malformed tables, RNA/DNA table equivalence, negative/nonfinite values, unusable groups, normalized duplicates, path collisions, and write failures.
- Mixed batches retry only affected records, preserve input/rank order, output valid peers, and aggregate exit status correctly. Same seed gives identical output across thread counts.
- Extend the Python auditor to map multiple candidate output IDs to each source record, independently check translation/CAI/GC/deviation/ranking/uniqueness and retry policy, and retain separately labeled legacy Biopython CAI results.

Run `cargo test`, `cargo clippy --all-targets`, and independent audit of real output fixtures. Measure runtime, success rates, candidate counts, and CAI improvements across lengths and retry settings. Preserve Rust 1.73 and existing dependency pins.

Document five-per-input output, CAI convention, strict threshold, retry defaults, effective deviation limits, warning/exit behavior, search limits, and a complete Python table export example. Completion means the example command executes the specified two-run policy and every emitted candidate passes hard validity checks, with low CAI clearly reported when necessary.

## Delivery order

1. Shared metrics, validity, retry-policy, and report types.
2. Composition exploration and multi-outcome arrangement collection.
3. Candidate validation, deduplication, deterministic ranking, and one-retry orchestration.
4. CLI, multi-candidate FASTA/report output, and end-to-end batch tests.
5. Independent audit, performance measurements, and documentation.
