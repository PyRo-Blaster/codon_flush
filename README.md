# Codon Harmony

Reverse-translates a protein and optimises the gene for a host's codon usage, removing
the sequence features that make genes hard to synthesise or express.

This repository holds a Rust engine (`src/`) that replaces the original Python
implementation's optimisation algorithm. The Python package (`codon_harmony/`) is kept
as the reference implementation and as an independent auditor of the Rust output.

## Why it was rewritten

The original algorithm resampled codons, harmonised them against the host profile once,
then ran six cleanup stages in a `cycles x inner_cycles` loop. Measured on GFP with an
*E. coli* table:

- **Time was quadratic**, exponent 2.0: 100 residues took 0.063 s and 1600 residues
  17.3 s at reduced settings, extrapolating to roughly seven minutes per sequence at the
  defaults. Most of it went to rescanning the whole gene per window.
- **The codon profile was never repaired.** Harmonisation ran once per sample, *before*
  six cleanup stages that mutate codons at random. Over 12 seeds, five produced nothing
  at all and all seven that succeeded ended at 12-16% deviation against a stated 10%
  limit. Twelve constraint violations survived in the finished gene, unreported.

## The idea

Deviation from the host profile and the codon adaptation index depend only on *which*
codons a gene uses. Every sequence constraint - GC windows, homopolymers, repeats,
hairpins, restriction sites, start sites, splice motifs - depends only on *where* they
are placed. The two halves separate:

- **Composition** is solved exactly, once, by largest-remainder apportionment against
  the host profile. This is the L1-optimal integer allocation.
- **Arrangement** is then searched using moves that swap two positions holding the same
  amino acid. Such a swap cannot change the codon multiset, so it cannot change the
  deviation or the CAI, and the search repairs constraints without ever disturbing the
  profile. Budget is spent on codon substitutions only when no swap can clear a
  violation.

Constraints are evaluated incrementally. Fixed-window checks do a bounded amount of
local work, but an index does not by itself make a non-local constraint constant time:
a 9-mer that occurs many times still has an occurrence list that grows with it. The
[correctness and performance design](docs/superpowers/specs/2026-10-17-engine-correctness-performance-design.md)
records which of those limits have been fixed and proposes bounded-memory replacements
for the rest.

## Results

Historical measurements on the same machine and host table. `scripts/validate.py`
independently recomputes translation, deviation and CAI; it does not validate every
sequence constraint or enforce the requested deviation limit.

| GFP, 238 residues | Python, 10x10 | Rust, one gene | Rust, 5 candidates |
|---|---|---|---|
| deviation from host profile | 0.1375 | **0.0083** | 0.1000 |
| codon adaptation index | 0.701 | 0.738 | **0.843** |
| produced a gene | 6 of 8 seeds | always | always |
| unresolved constraints | 12 | **0** | **0** |
| wall time | ~7 s | **31 ms** | 0.93 s |

The last column is the current default: the composition search spends the deviation
limit it was given to buy CAI, arranges eight multisets instead of one, and returns the
five best distinct genes. It costs a second per gene rather than 31 ms, and it is the
deviation limit, not the apportionment optimum, that the result sits at. `--max-relax 0`
still produces the middle column's behaviour.

Scaling, single thread, random proteins of realistic composition: 1600 aa 1.05 s,
3200 aa 2.04 s, 6400 aa 4.55 s, 12800 aa 9.38 s - a fitted exponent of **1.05** against
the Python's 2.0.

Anything the engine cannot fix is reported with its position and cause, and the process
exits non-zero. A search that stops is not a proof of infeasibility, and other failed
cases may still have a feasible arrangement. Some are provably impossible: three
consecutive methionines are ATG ATG ATG in every encoding, and four lysines with three
AAA in the multiset give six consecutive A's however they are arranged.

## Five candidates per protein, with one retry

Every arrangement of one codon multiset has the same CAI, so a search that only
rearranges cannot optimise it. The engine therefore explores several multisets as well:
a bounded beam search over synonymous one-count transfers, starting from the host
apportionment and from the all-preferred-codon allocation, keeping the compositions
that stay inside the deviation limit and the GC band. Each is then arranged, and the
finished genes are ranked.

For each protein the CLI returns `--candidates` distinct valid genes, best CAI first.
If none of them exceeds `--cai-threshold`, the whole search runs **once** more with
more deviation allowed and twice the effort, and the best candidates of both runs are
returned. There are at most two runs per protein.

```
./target/release/codon_harmony \
    --input proteins.fasta --output genes.fasta \
    --host-table tests/data/e_coli_316407.json \
    --candidates 5 --cai-threshold 0.80 \
    --max-relax 0.10 --retry-max-relax 0.20 \
    --restarts 4 --retry-restarts 8 \
    --report-json report.json
```

| | first run | retry |
|---|---|---|
| deviation limit | `--max-relax`, default 0.10 | `--retry-max-relax`, default `min(1, limit + 0.10)` |
| arrangements per composition | `--restarts`, default 4 | `--retry-restarts`, default twice that |
| composition search budget | `--composition-budget`, default 2000 | twice that |

These are engineering defaults, not biologically validated thresholds, and nothing
guarantees the retry improves anything. The effective policy is printed before the run
and recorded in the report.

The threshold is a retry trigger and a quality warning, never a filter. Equality at the
threshold does **not** pass; the comparison is made on the unrounded value. A complete
set whose best CAI is low is a successful fallback: it is written, a warning naming the
best CAI goes to stderr, and the exit status stays zero. What does exit non-zero is a
record with **fewer** candidates than asked for (`insufficient_candidates`) or none at
all (`no_valid_candidates`) - a sequence is never duplicated to pad a set. Five distinct
genes do not always exist; `MWMWM` has exactly one. A bounded search that finds fewer
is reported as a shortage, not as proof that no other encoding exists.

Valid means the gene translates exactly, uses only allowed codons, stays inside the
deviation bound and the GC band, and has no unresolved sequence constraint. The
low-CAI fallback relaxes the CAI expectation and nothing else.

Output IDs are `record_000001_candidate_01`: the input ordinal, not the name, since two
input records may share an ID. Each description carries `source=`, the rank, the CAI,
the deviation, the GC and which run produced it; `--report-json` writes the same with
full precision, plus the policy, both runs' settings and outcomes, and every
candidate's provenance.

A short protein is the one case where the deviation bound is not `--max-relax`: the
largest-remainder baseline is L1-optimal, so when its integer rounding alone costs more
than the limit it is still admitted, and the bound that actually applied is reported as
`dev_limit` rather than the limit being silently enlarged.

## Usage

```
cargo build --release
./target/release/codon_harmony \
    --input proteins.fasta \
    --output genes.fasta \
    --host-table tests/data/e_coli_316407.json
```

Useful options:

- `--max-relax 0` freezes the profile: no substitution is ever made, so deviation stays
  at the apportionment optimum and the CAI is identical across every seed and restart.
  It stays frozen on the retry, which then buys arrangement effort only; an explicit
  positive `--retry-max-relax` authorises unfreezing.
- `--start-site-rule` chooses how an internal start site is recognised. The default,
  `shine-dalgarno`, asks whether the upstream sequence could actually pair with the 16S
  rRNA. `purine-run` restores the original rule.
- `--seed`, `--restarts`, `--threads` control the search. A seed reproduces a run
  exactly, whatever the thread count.

Host tables are JSON, in either the nested `{"A": {"GCA": 0.21, ...}}` shape that
`python_codon_tables` produces or a flat codon map. RNA spellings are accepted and
normalised to DNA, and a codon filed under the wrong amino acid is rejected rather
than renormalised into the wrong group. The table must be complete: every amino acid
and the stop codons need at least one codon the threshold leaves usable.
`tests/data/e_coli_316407.json` is a full working export and doubles as the example.

## Verifying the output

`scripts/validate.py` reads a protein FASTA and the genes the Rust tool produced, maps
each candidate back to its source protein through the `source=` field rather than
trusting the tool's grouping, and independently re-checks translation, deviation, GC
and CAI using Biopython and the original Python package, so nothing is self-certified.
Given the JSON report as well, it also checks that each set is distinct, ranked, no
longer than was asked for, within the deviation bound that applied, and that the retry
happened exactly when the declared policy says it should. It needs Biopython < 1.78.

```
python scripts/validate.py proteins.fasta genes.fasta tests/data/e_coli_316407.json
python scripts/validate.py proteins.fasta genes.fasta table.json report.json
```

The CAI it prints is Biopython's, which divides the log sum by `len - 1` where the Rust
engine uses `len`. The two are labelled separately and are expected to differ in the
last digits rather than to match.

## Tests

`cargo test` runs 114 tests. The backbone is differential: every constraint detector is
checked against a brute-force oracle on random sequences, and its incremental
bookkeeping is checked against a full rebuild after every move across 20,000 random
moves. Integration tests assert that the protein always survives, that the deviation
budget is never exceeded, that swaps leave deviation and CAI bit-identical, that runs
are reproducible from a seed, and that the search only stops where no single move
helps. End-to-end tests in `tests/cli.rs` run the binary and check the candidate FASTA
against the JSON report: five distinct genes per protein, ranked, each translating back;
exactly one retry when nothing clears the threshold and none when something does; an
honest shortage and a non-zero exit for a protein with one encoding; a failing record
that does not stop its peers; and identical output at one thread and at four.

## Design documents

- [Correctness and performance proposal](docs/superpowers/specs/2026-10-17-engine-correctness-performance-design.md) -
  corrections to the validity contract, bounded-memory indices, coordinated repair,
  delivery phases and validation gates; [reproducible evidence](docs/benchmarks/2026-10-17-engine-review/README.md).
- [Design](docs/superpowers/specs/2026-10-17-codon-harmony-rust-redesign-design.md) -
  the measured problems, the approach, the constraint definitions and the results.
- [Implementation plan](docs/superpowers/plans/2026-10-17-codon-harmony-rust-redesign.md) -
  the task breakdown the rewrite followed.

## Credit

Derived from the original Codon Harmony by Brian Weitzner and contributors, MIT
licensed; see `AUTHORS.rst` and `LICENSE`. The original README is kept as `README.rst`.
