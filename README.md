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

Constraints are evaluated incrementally. Fixed-window checks have bounded local work,
but repeat/hairpin indices and diagnostic queues still have costs that grow with
occurrence multiplicity. The [correctness and performance design](docs/superpowers/specs/2026-10-17-engine-correctness-performance-design.md)
documents those limits and proposes bounded-memory replacements.

## Results

Historical measurements on the same machine and host table. `scripts/validate.py`
independently recomputes translation, deviation and CAI; it does not validate every
sequence constraint or enforce the requested deviation limit.

| GFP, 238 residues | Python, 10x10 | Rust |
|---|---|---|
| deviation from host profile | 0.1375 | **0.0083** |
| codon adaptation index | 0.701 | 0.738 |
| produced a gene | 6 of 8 seeds | always |
| unresolved constraints | 12 | **0** |
| wall time | ~7 s | **31 ms** |

Scaling, single thread, random proteins of realistic composition: 1600 aa 1.05 s,
3200 aa 2.04 s, 6400 aa 4.55 s, 12800 aa 9.38 s - a fitted exponent of **1.05** against
the Python's 2.0.

The current implementation has known correctness limitations: an initial composition
can exceed a positive `--max-relax` limit, and stems with more than six G bases can
have mixed wobble matches missed by the hairpin detector. The linked design includes
reproducers, correction requirements and benchmark evidence. Its proposed fixes are
not yet implemented; the historical results above are not a guarantee of validity for
every input or configuration.

Detected unresolved constraints are reported with their positions, and the process
exits non-zero. A failed search does not prove infeasibility. Some requirements are
unavoidable: three consecutive methionines force three identical ATG codons, and four
lysines with three AAA in the multiset give six consecutive A's however they are arranged.

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
- `--start-site-rule` chooses how an internal start site is recognised. The default,
  `shine-dalgarno`, asks whether the upstream sequence could actually pair with the 16S
  rRNA. `purine-run` restores the original rule.
- `--seed`, `--restarts`, `--threads` control the search. A seed reproduces a run
  exactly, whatever the thread count.

Host tables are JSON, in either the nested `{"A": {"GCA": 0.21, ...}}` shape that
`python_codon_tables` produces or a flat codon map.

## Verifying the output

`scripts/validate.py` reads a protein FASTA and the genes the Rust tool produced, then
independently re-checks translation, deviation and CAI using Biopython and the original
Python package, so nothing is self-certified. It needs Biopython < 1.78.

```
python scripts/validate.py proteins.fasta genes.fasta tests/data/e_coli_316407.json
```

## Tests

`cargo test` runs 80 tests. The backbone is differential: every constraint detector is
checked against a brute-force oracle on random sequences, and its incremental
bookkeeping is checked against a full rebuild after every move across 20,000 random
moves. Integration tests assert that the protein always survives, that the deviation
budget is never exceeded, that swaps leave deviation and CAI bit-identical, that runs
are reproducible from a seed, and that the search only stops where no single move helps.

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
