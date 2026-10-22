# Codon Harmony

Reverse-translates a protein and optimises the gene for a host's codon usage, removing
the sequence features that make genes hard to synthesise or express.

This repository holds a Rust engine (`src/`) that replaces the original Python
implementation. That implementation has now been removed; what survives of it is
`scripts/validate.py`, a standalone auditor that re-derives the engine's numbers from
its output and checks its claims against them.

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

Both limits are hard, for every candidate including the low-CAI fallback. A positive
`--max-relax` bounds every composition and every substitution; the largest-remainder
baseline is L1-optimal, so when even its integer rounding costs more than the limit,
the run fails saying what the floor is rather than raising the bound to fit. The GC
band is hard throughout composition search, arrangement and acceptance; a protein that
cannot reach it within the deviation limit is an error, never a widened band. Only the
CAI expectation is ever relaxed. `--max-relax 0` keeps its own meaning: the multiset is
frozen and the integer-rounding deviation it already has is reported as `dev_limit`.

The candidate count is a target the run works at, not a count of restarts. One frozen
multiset with four restarts cannot yield five distinct genes, so a run that is short
keeps arranging on fresh streams, up to a bounded surplus, before it reports a
shortage. That is not a retry: the one-retry CAI rule is untouched.

## Usage

```
cargo build --release
./target/release/codon_harmony \
    --input proteins.fasta \
    --output genes.fasta \
    --host-table tests/data/e_coli_316407.json
```

Useful options:

- `--splice-rule` chooses what counts as a splice problem. The default, `no-motif`,
  rejects every donor and every acceptor. `no-pair` rejects only a donor and an
  acceptor *together*. See *Splice sites that no encoding can avoid* below.
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

## Splice sites that no encoding can avoid

Thirty tripeptides force a splice donor in **every** synonymous encoding:

```
{M, W} x {C, W, Y} x {A, D, E, G, V}
```

M is `ATG` and W is `TGG`, so both end in G with no alternative; C, W and Y all begin
`T` with A or G second; A, D, E, G and V all begin with G. The nine-base window is
therefore forced to `??G|T[AG]?|G??`, which matches the donor consensus `..GT[AG].G.`
whatever synonyms are chosen. The GT that defines a donor straddles a codon boundary,
where no synonymous substitution reaches. This is a property of the genetic code, not
of the host table: the same 30 come out with every codon available. Nothing forces an
*acceptor* — across all 160,000 four-residue peptides and every encoding of each,
there is always an escape.

Antibody constant regions contain these tripeptides, so under the default rule the
engine correctly reports that it cannot build them. `--splice-rule no-pair` is the way
through. Splicing needs a pair: a 5' donor and a downstream 3' acceptor define an
intron, so a lone half site cannot be spliced out of the coding sequence. Under
`no-pair` the engine clears every acceptor — always possible — and keeps the donors it
cannot avoid. `grep -oE "[MW][CWY][ADEGV]"` on a protein finds them.

**The caveat, which is why this is opt-in.** A lone donor is inert only relative to the
coding sequence. A transgene is expressed from a vector, and a donor inside the insert
can still pair with an acceptor downstream in the vector's expression cassette — the
usual mammalian arrangement, and exactly the context `--remove-splice-sites` is for. So
every retained site is reported with its position and kind, in the FASTA description
(`retained_splice_sites=donor@27,donor@117`) and per candidate in the JSON, so the gene
can be checked against the vector it is going into.

A second cost: `no-pair` does not *minimise* half sites, it only forbids pairs. Once a
gene has no acceptors the objective is flat, so the search stops removing donors it
could have avoided — GFP is spotless under the default and carries lone donors under
`no-pair`. Among candidates of equal CAI and deviation the one with fewer is ranked
first, but the rule is for proteins that need it, not a better default.

Counting is per site, not per pattern. The 20 donor patterns overlap heavily, so one
eight-base window can match several; it is one site, and is counted and reported once.

## Verifying the output

`scripts/validate.py` reads a protein FASTA and the genes the Rust tool produced, and
maps each candidate back to its source protein through the `source=` field rather than
trusting the tool's grouping, so nothing is self-certified. Every number it prints is
re-derived from the sequence — translation, GC, CAI and
deviation — and the claims in the FASTA descriptions and in the report are then checked
*against* those, so a wrong CAI, a wrong ranking or an out-of-band sequence fails the
audit instead of passing it. Given the report it also checks set size, distinctness,
ranking by the recomputed CAI, the requested GC band, the deviation bound that applied,
and whether the retry fired exactly when the declared policy says it should.

It needs **only the Python standard library**, so it runs anywhere and CI runs it. Each
metric is a second implementation written from its definition rather than shared with
the engine, which is what makes the check worth running; on real output the two agree
to the last bit. Deviation follows the original tool's `compare_profiles` at
`relax=1.0`.

`tests/test_validate_audit.py` feeds the auditor deliberately falsified reports and
requires each to be caught.

```
python scripts/validate.py proteins.fasta genes.fasta tests/data/e_coli_316407.json
python scripts/validate.py proteins.fasta genes.fasta table.json report.json
```

A protein the report accounts for with `no_valid_candidates` or `error` is reported as
such rather than as missing output; any other absence is a failure. An `error` record
never searched, so it is required only to claim nothing. Every input protein must have
a record, so a report cannot drop one to escape being checked.

## Tests

`cargo test` runs 138 tests, and `pytest tests/test_validate_audit.py` runs 23 more
against the auditor itself. The backbone is differential: every constraint detector is
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
