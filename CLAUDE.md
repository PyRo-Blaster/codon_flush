# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Codon Harmony** reverse-translates a protein and optimises the gene for a host's codon
usage, removing features that make genes hard to synthesise or express.

The optimiser is the **Rust engine in `src/`**. The Python package in `codon_harmony/`
is the original implementation, kept as the reference and as an independent auditor of
the Rust output; it is no longer the engine. Both are MIT licensed, derived from the
original Codon Harmony.

## Commands

### Rust engine (the thing you are usually working on)

```bash
cargo build --release
cargo test                        # 80 tests, about 2 s
cargo clippy --all-targets        # expected to be clean
./target/release/codon_harmony \
    --input proteins.fasta \
    --output genes.fasta \
    --host-table tests/data/e_coli_316407.json
```

Useful flags: `--max-relax 0` freezes the codon profile (no substitutions at all),
`--start-site-rule purine-run` restores the original start-site definition, and
`--seed` / `--restarts` / `--threads` control the search. `--cycles` and
`--inner-cycles` are accepted and ignored, with a warning: they described the old
sampling loop.

Host tables are JSON, either nested `{"A": {"GCA": 0.21, ...}}` as produced by
`python_codon_tables`, or a flat codon map. Frequencies are renormalised per amino acid,
so a table whose single-codon amino acids do not sum to 1 is repaired rather than
trusted. There is no table download: pass a file.

### Auditing the output with the Python

```bash
python scripts/validate.py proteins.fasta genes.fasta tests/data/e_coli_316407.json
```

Re-derives translation, deviation and CAI from the Rust output using Biopython and the
original package, so results are never self-certified.

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
src/composition.rs                 apportionment, GC band, budget
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
  constraint O(1). Occurrence-list traversal, long runs and expanded reports still
  grow with sequence length or multiplicity; see the corrective design below.
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
