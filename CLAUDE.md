# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Codon Harmony** reverse-translates a protein and optimises the gene for a host's codon
usage, removing features that make genes hard to synthesise or express. The CLI returns
several ranked candidates per protein, not one gene; see *Candidates and the retry*.

The optimiser is the **Rust engine in `src/`**. The Python package in `codon_harmony/`
is the original implementation, kept as the reference and as an independent auditor of
the Rust output; it is no longer the engine. Both are MIT licensed, derived from the
original Codon Harmony.

## Commands

### Rust engine (the thing you are usually working on)

```bash
cargo build --release
cargo test                        # 126 tests, about 5 min in debug
cargo clippy --all-targets        # expected to be clean
./target/release/codon_harmony \
    --input proteins.fasta \
    --output genes.fasta \
    --host-table tests/data/e_coli_316407.json \
    --report-json report.json
```

Useful flags: `--candidates` and `--cai-threshold` size the candidate set and set the
retry trigger, `--retry-max-relax` / `--retry-restarts` / `--composition-budget` tune
the two runs, `--report-json` writes the machine-readable report,
`--max-relax 0` freezes the codon profile (no substitutions at all),
`--start-site-rule purine-run` restores the original start-site definition, and
`--seed` / `--restarts` / `--threads` control the search. `--cycles` and
`--inner-cycles` are accepted and ignored, with a warning: they described the old
sampling loop.

`.github/workflows/ci.yml` runs exactly these three on every push to `main` and every
pull request, on 1.73 with `--locked` and `-D warnings`. Clippy is clean on 1.73; a
newer clippy warns about `map_or` calls it would write as `is_none_or`, which is why
the toolchain is pinned rather than floating.

Host tables are JSON, either nested `{"A": {"GCA": 0.21, ...}}` as produced by
`python_codon_tables`, or a flat codon map. Frequencies are renormalised per amino acid,
so a table whose single-codon amino acids do not sum to 1 is repaired rather than
trusted. There is no table download: pass a file.

### Auditing the output with the Python

```bash
python scripts/validate.py proteins.fasta genes.fasta tests/data/e_coli_316407.json
python scripts/validate.py proteins.fasta genes.fasta table.json report.json
```

Re-derives every number from the sequences and checks the engine's claims *against*
them, rather than printing both and trusting the claim: a falsified CAI, a wrong
ranking or an out-of-band sequence must fail the audit. Candidates map back to their
source protein through the `source=` field, not through the Rust tool's grouping.
Biopython's CAI divides by `len - 1` where the engine uses `len`, so it is a separate
column and never what a claim is checked against; the engine-convention CAI is
recomputed in plain Python from the host table and agrees with the Rust value to the
last bit. The metric and audit helpers import nothing from Biopython, so
`pytest tests/test_validate_audit.py` can exercise them on falsified input anywhere.

### Python package

```bash
pytest tests/        # 37 pass, 3 fail
flake8
```

The three failures are pre-existing and unrelated to this work: two need `misc/`
fixtures that are not in this copy of the repository, and one expects a network 404 for
a bogus taxonomy ID.

## Environment constraints

- **Rust 1.73** is the toolchain here, which is older than several crates' MSRV.
  `rayon` is pinned to `=1.8.1` with `rayon-core` at 1.12.1, `half` to 2.2.1 and
  `criterion` to `=0.5.1`. Bumping them breaks the build. `Option::is_none_or` and other
  recent std additions are unavailable.
- The Python side needs **Biopython < 1.78**, because it uses `Bio.Alphabet`.

## Architecture

The problem splits in two, and that split is the whole design:

- **Composition** (`src/composition.rs`) decides *which* codons the gene uses and how
  many of each. That alone fixes the deviation from the host profile, the CAI and the
  whole-gene GC. Solved exactly, once, by largest-remainder apportionment.
- **Arrangement** (`src/arrange/`) decides *where* each codon goes. That alone fixes
  every sequence constraint. Searched by min-conflicts local search.

The default move swaps the codons at two positions holding the same amino acid, which
cannot change the multiset and so cannot change deviation or CAI. Substitutions, which
do change the multiset, are charged against the `--max-relax` budget and are only made
when no swap can clear a violation.

```
src/code.rs                        base and codon encoding, genetic code
src/host.rs                        table parsing, thresholding, CAI weights
src/composition.rs                 apportionment, GC band, budget, composition search
src/batch.rs                       candidate sets, ranking, the one-retry policy
src/metrics.rs                     deviation, CAI, GC of a finished gene
src/arrange/state.rs               SeqState, Move, and the View/Patch machinery
src/arrange/search.rs              the search driver
src/arrange/constraints/           one file per constraint, plus the Constraint trait
src/fasta.rs  src/rng.rs  src/small.rs
```

## Invariants to preserve when editing

- **`Constraint::delta` must not mutate anything.** It is handed a `before` and an
  `after` view of the same pre-move state and must return the exact change in the
  violation total. `apply` then folds the same move in.
- **Evaluate affected dependencies incrementally and state their actual cost.** Fixed
  windows have bounded local work, but an index does not by itself make a nonlocal
  constraint O(1): occurrence-list traversal and long runs still grow with multiplicity.
  Reports are linear in the violations as of commit 3402d01. See the corrective design.
- **Touched positions are kept sorted** (`sorted_positions` in `state.rs`). Constraints
  merge adjacent regions as they walk that list, so an unsorted pair silently drops a
  region. This caused every incremental test to fail once already.
- **Fixed-capacity buffers** (`src/small.rs`) have capacities derived from how far a move
  can reach. A full buffer means an assumption changed, not that the capacity should grow.
- Every constraint needs **both** differential tests from
  `src/arrange/constraints/mod.rs`: `assert_matches_oracle` against a brute-force
  reference, and `assert_incremental_matches_rebuild` over tens of thousands of moves.
  These are the backbone of the suite; do not add a constraint without them.
- Restarts are reduced **deterministically** by index, so output never depends on thread
  scheduling. A seed must reproduce a run exactly.

## Candidates and the retry

Arrangement cannot change CAI, because every arrangement of one multiset has the same
one. `composition::explore` therefore searches multisets too: a bounded beam over
synonymous one-count transfers from the apportioned baseline and from the
all-preferred-codon allocation, with a paired GC-repair proposal so a one-transfer
barrier does not stop it. Intermediate states may break the limits; only compositions
inside the deviation bound and the GC band are handed to arrangement, and the baseline
is always among them, since arrangement can fail for one multiset and succeed for
another. A spent budget is `search_exhausted`, never a claim that nothing else exists.

`batch::run_record` arranges every offered composition (`optimise_all` keeps all the
restarts rather than reducing to one), re-derives every metric from the finished
sequence, admits only genes that translate exactly and leave no violation, deduplicates
by complete DNA, and ranks by descending CAI, then ascending deviation, then the
sequence itself. If nothing strictly exceeds the threshold it runs the whole thing once
more under the retry policy and pools both runs, so a retry can never replace a better
first-run candidate with a worse one. Two runs per protein, never three.

Invariants here: the threshold is a trigger and a warning, never a filter, and equality
does not pass; a complete low-CAI set exits zero while a short set or an empty one
exits non-zero; a sequence is never duplicated to pad a set; random streams come from
`stream_seed`, which mixes the base seed with the record, run, composition and round
**indices** only, so output never depends on thread scheduling.

**The deviation limit and the GC band are hard, everywhere.** A positive `--max-relax`
bounds every composition, every substitution (`best_substitution` checks both) and
final acceptance. The baseline is not exempt: when its rounding floor exceeds a
positive limit the answer is `CompositionError::DeviationFloor`, not a wider limit, and
nothing is written. The band is hard the same way, with `GcOutsideBudget` when it
cannot be reached within the limit. Only the CAI expectation is ever relaxed, and
`--max-relax 0` still freezes the multiset and reports its rounding deviation. Reviews
have caught both of these being quietly widened; do not reintroduce
`max_relax.max(baseline.deviation)` or a GC allowance inherited from the baseline.

A run that is short of candidates keeps arranging on fresh streams up to
`ARRANGEMENT_SURPLUS * wanted` attempts, because one frozen multiset with four restarts
cannot give five distinct genes. That is not a retry; the one-retry CAI rule is
separate.

## Constraint definitions

GC uses true sliding 20/50/100 nt windows; whole-gene GC belongs to composition, where
no rearrangement can disturb it. A homopolymer is a run of identical bases, not the
commonest base's count in a 6-nt window. An internal start site is, by default, four
contiguous bases of the `AGGAGG` consensus placed 4-14 bases ahead of an ATG, GTG or
TTG; the original purine-run rule is still selectable. Repeats, hairpins, restriction
sites and splice motifs keep their original meanings.

Constraints the engine cannot satisfy are reported with position and cause, and the
process exits non-zero. A capped search failure is not proof of infeasibility: other
failed cases may still have a feasible arrangement. Some are provably impossible, such
as three consecutive methionines, which are ATG ATG ATG in every encoding, and four
lysines with three AAA in the multiset, which give six consecutive A's however they are
arranged.

## Design documents

- `docs/superpowers/specs/2026-10-17-engine-correctness-performance-design.md` — proposed
  corrections to validity and complexity claims, bounded-memory scoring/scheduling,
  coordinated repair and phased validation. This is a proposal, not implemented behavior.
- `docs/benchmarks/2026-10-17-engine-review/README.md` — baseline/prototype evidence and
  reproduction instructions, including known positive-limit and hairpin failures.
- `docs/superpowers/specs/2026-10-17-codon-harmony-rust-redesign-design.md` — the
  measured problems with the original, the approach, and the results.
- `docs/superpowers/plans/2026-10-17-codon-harmony-rust-redesign.md` — the task
  breakdown the rewrite followed.
