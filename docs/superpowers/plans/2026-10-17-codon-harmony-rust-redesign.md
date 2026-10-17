# Codon Harmony Rust Redesign Implementation Plan

> **Execution note:** this plan is executed inline, in the session that wrote it and the
> spec (`docs/superpowers/specs/2026-10-17-codon-harmony-rust-redesign-design.md`), so
> tasks carry interfaces, algorithms, test cases and exact commands rather than a full
> transcription of every line of Rust. Steps use checkbox (`- [ ]`) syntax for tracking.
> The repository is **not** a git repository, so each task ends in a verification
> checkpoint (tests green, `cargo clippy` clean) instead of a commit.

**Goal:** Replace the `src/` Rust port with an engine that optimises codon choice in
`O(n)` expected time and lands within ~1% of the host codon profile instead of 12-16%.

**Architecture:** Split the problem in two. *Composition* (which codons, how many)
fixes deviation and CAI exactly, once, by largest-remainder apportionment. *Arrangement*
(where each codon goes) is then searched by min-conflicts local search using only
same-amino-acid swaps, which provably cannot change the multiset — so deviation and CAI
stay frozen at their optimum while every sequence constraint is repaired. All constraint
checks are incremental, O(1) per move.

**Tech Stack:** Rust 1.73, clap 4.4 (derive), rayon, serde/serde_json, criterion.
No regex: splice motifs use bitap over character classes.

---

## File structure

| file | responsibility |
|---|---|
| `src/lib.rs` | public API, re-exports |
| `src/host.rs` | codon table parse (both JSON shapes), threshold, renormalise, CAI weights |
| `src/composition.rs` | apportionment, whole-gene GC adjustment, deviation budget |
| `src/metrics.rs` | deviation, CAI, GC of a finished sequence |
| `src/arrange/state.rs` | `SeqState`, `Move`, swap/substitute primitives |
| `src/arrange/search.rs` | min-conflicts driver, escalation, restarts |
| `src/arrange/constraints/mod.rs` | `Constraint` trait, registry, `Violation` |
| `src/arrange/constraints/gc.rs` | sliding 20/50/100 nt GC windows |
| `src/arrange/constraints/homopolymer.rs` | runs of identical bases |
| `src/arrange/constraints/repeat.rs` | identical codon triples, duplicated 9-mers |
| `src/arrange/constraints/hairpin.rs` | 10-nt stems vs reverse-complement partners |
| `src/arrange/constraints/restriction.rs` | enzyme recognition sites |
| `src/arrange/constraints/start_site.rs` | purine runs upstream of xTG |
| `src/arrange/constraints/splice.rs` | 30 donor/acceptor motifs via bitap |
| `src/arrange/reference.rs` | brute-force O(n^2) detectors, `#[cfg(test)]` only |
| `src/main.rs` | CLI, FASTA I/O, rayon over records, reporting |
| `benches/optimize.rs` | criterion benchmarks, n = 100..10000 |

Delete: `src/optimize.rs`, `src/codon.rs` (replaced), `src/data/ecoli_codons.json`
(corrupt: TGG = 0.48 for single-codon Trp).

---

## Shared encodings (used by every task)

```rust
// base: A=0, C=1, G=2, T=3.  codon index = b0*16 + b1*4 + b2, so 0..64.
pub type Base = u8;
pub type CodonIdx = u8;
pub type AaIdx = u8;      // 0..20, 20 = stop
```

A 9-mer packs into 18 bits and a 10-mer into 20 bits; both are built by
`(prev << 2 | base) & mask`, so rolling updates need no hashing.

---

### Task 1: Crate scaffold and host table

**Files:** Create `src/host.rs`; rewrite `src/lib.rs`, `Cargo.toml`; delete
`src/optimize.rs`, `src/codon.rs`, `src/data/ecoli_codons.json`.

- [ ] **Step 1: Cargo.toml** — deps `clap 4.4.18` (derive), `rayon 1.8`, `serde 1` (derive),
  `serde_json 1`, `anyhow 1`; dev-dep `criterion 0.5`; keep `[profile.release] opt-level=3, lto=true`.
- [ ] **Step 2: Write failing tests** in `src/host.rs`:
  - `parses_nested_table`: `{"A": {"GCA": 0.21, ...}, ...}` loads and `freq("GCA") == 0.21`.
  - `parses_flat_table`: `{"GCA": 0.21, ...}` loads to the same thing.
  - `renormalises_per_amino_acid`: a table where Trp `TGG` is 0.48 yields 1.0, because
    each amino acid's frequencies are divided by their own sum.
  - `threshold_zeroes_and_renormalises`: threshold 0.10 zeroes E. coli `AGA`(0.04),
    `AGG`(0.02), `ATA`, `CGA`, `CTA`; the surviving Arg codons sum to 1.0.
  - `cai_weights_are_relative_to_group_max`: weight of the most-used codon in each
    group is 1.0.
- [ ] **Step 3: Run** `cargo test --lib host` — expect failure (no `host` module).
- [ ] **Step 4: Implement** `HostTable { freq: [f64; 64], weight: [f64; 64], by_aa: [Vec<CodonIdx>; 21] }`
  with `load(path, threshold)`, `freq`, `weight`, `synonyms(aa)`. Parse into
  `serde_json::Value`, accept either shape, group by amino acid via the standard code,
  renormalise, apply threshold, renormalise again, then compute weights.
- [ ] **Step 5: Run** `cargo test --lib host` — expect PASS.
- [ ] **Step 6: Checkpoint** — `cargo clippy --all-targets -- -D warnings` clean.

**Cross-check:** `tests/data/e_coli_316407.json` at threshold 0.10 must produce exactly
the profile the Python produces; the known-good values are in
`viz-demo/pydata.json` (`hostProfile`, `relativeAdaptiveness`) in the shared notes
folder. Assert a handful of them directly.

---

### Task 2: Composition

**Files:** Create `src/composition.rs`.

- [ ] **Step 1: Write failing tests**
  - `apportionment_sums_to_group_size`: for random profiles and n in 1..200, the counts
    for each amino acid sum to exactly its position count.
  - `apportionment_beats_rounding_on_l1`: largest-remainder L1 distance to `f*n` is
    always <= that of per-codon `round()`.
  - `deviation_matches_python`: for GFP's amino acid composition under the E. coli
    table, deviation is <= 0.01 (the Python's own metric; prototype measured 0.0083).
  - `gc_adjust_moves_into_band`: a synthetic profile whose apportionment lands at GC
    0.72 is pulled under 0.65, and the deviation it costs is reported.
  - `gc_adjust_reports_infeasible`: an all-Trp protein (GC fixed at 0.667 by TGG) is
    reported infeasible for a band of 0.30..0.40 rather than looping.
- [ ] **Step 2: Run** `cargo test --lib composition` — expect failure.
- [ ] **Step 3: Implement**
  ```rust
  pub struct Composition { pub counts: [u32; 64], pub deviation: f64, pub budget_left: f64 }
  pub fn build(protein: &[AaIdx], host: &HostTable, max_relax: f64, gc_band: (f64, f64))
      -> Result<Composition, CompositionError>;
  ```
  Largest-remainder per amino acid; then while whole-gene GC is outside the band, apply
  the unit move (one count from codon `x` to synonym `y`) maximising
  `|delta_gc| / delta_deviation`, stopping when in band, when the budget is gone, or
  when no move improves GC (infeasible).
- [ ] **Step 4: Run** `cargo test --lib composition` — expect PASS.
- [ ] **Step 5: Checkpoint** — clippy clean.

---

### Task 3: Sequence state and moves

**Files:** Create `src/arrange/state.rs`, `src/arrange/mod.rs`.

- [ ] **Step 1: Write failing tests**
  - `swap_preserves_translation_and_multiset`: after 10,000 random swaps the translation
    equals the input protein and the codon count vector is unchanged.
  - `substitute_changes_one_count_pair`: `Subst(p, c)` decrements the old codon's count
    and increments `c`, and nothing else.
  - `nt_stays_in_sync`: `state.nt` always equals the codons expanded to bases.
- [ ] **Step 2: Run** `cargo test --lib state` — expect failure.
- [ ] **Step 3: Implement**
  ```rust
  pub enum Move { Swap { p: u32, q: u32 }, Subst { p: u32, codon: CodonIdx } }
  pub struct SeqState { pub codons: Vec<CodonIdx>, pub nt: Vec<Base>,
                        pub counts: [u32; 64], pub pos_by_aa: Vec<Vec<u32>>, pub aa: Vec<AaIdx> }
  impl SeqState { fn apply(&mut self, mv: Move); fn touched(&self, mv: Move) -> [Option<u32>; 2]; }
  ```
  `apply` writes codons, the 3 bases each, and the counts. `touched` returns the codon
  positions a move changes, which is what every constraint's `delta` is scoped to.
- [ ] **Step 4: Run** `cargo test --lib state` — expect PASS.
- [ ] **Step 5: Checkpoint** — clippy clean.

---

### Task 4: Constraint trait and brute-force reference

**Files:** Create `src/arrange/constraints/mod.rs`, `src/arrange/reference.rs`.

- [ ] **Step 1: Define the trait** exactly as in spec section 4, plus
  ```rust
  pub struct Violation { pub kind: Kind, pub at: u32, pub len: u32, pub detail: String }
  pub enum Kind { Gc, Homopolymer, Repeat, Hairpin, Restriction, StartSite, Splice }
  ```
- [ ] **Step 2: Implement `reference.rs`** — naive O(n^2) detectors returning
  `Vec<Violation>` for all seven kinds, written for obvious correctness, not speed.
  These are the oracle for every later task.
- [ ] **Step 3: Write the shared differential harness**
  ```rust
  #[cfg(test)]
  pub fn assert_matches_reference<C: Constraint>(c: &mut C, kind: Kind, seqs: usize);
  #[cfg(test)]
  pub fn assert_incremental_equals_rebuild<C: Constraint>(c: &mut C, moves: usize);
  ```
  The second applies random moves, calling `apply` each time, and asserts the
  incremental counters equal a fresh `rebuild` after every move.
- [ ] **Step 4: Run** `cargo test --lib constraints` — the harness compiles, no
  constraints registered yet.
- [ ] **Step 5: Checkpoint** — clippy clean.

---

### Tasks 5-11: the seven constraints

Each follows the same five steps, so they are listed once and the per-constraint detail
follows:

- [ ] **Step 1:** write the constraint's unit tests (specific hand-built sequences,
  listed below) plus `assert_matches_reference(.., 200 random sequences)` and
  `assert_incremental_equals_rebuild(.., 100_000 moves)`.
- [ ] **Step 2:** run `cargo test --lib <name>` — expect failure.
- [ ] **Step 3:** implement `rebuild`, `violations`, `delta`, `apply`, `report`.
- [ ] **Step 4:** run `cargo test --lib <name>` — expect PASS.
- [ ] **Step 5:** checkpoint — clippy clean.

**Task 5 — `gc.rs`.** Window sets 20, 50, 100 nt with bands (0.15,0.90), (0.15,0.80),
(0.28,0.68). State: one `Vec<u16>` of per-window GC counts per set, plus a violation
counter per set. A base change at `i` shifts windows starting in `[i-len+1, i]`; update
that range and recount those windows only. Unit test: a window at exactly the band edge
is not a violation (`<` and `>`, not `<=`/`>=`), matching the Python.

**Task 6 — `homopolymer.rs`.** New semantics: a run of more than `N` consecutive
identical bases. State: `run[i]` = length of the run ending at `i`, plus a count of
maximal runs over `N`. A change re-walks only the affected maximal run. Unit tests:
`AAAAA` with N=4 is one violation; `GCGCGC` is **zero** (the old rule flagged it);
`AAAAAAAA` with N=4 is still one violation, not four.

**Task 7 — `repeat.rs`.** Two rules: three identical consecutive codons, and any 9-mer
occurring at two positions at least 9 apart. State: `[u32; 262144]` count table indexed
by the packed 9-mer, a `Vec` of first/last position per key, and a violation counter. A
codon change touches the 9 k-mers starting in `[p-8, p]`. Unit test: `ACCGAGTCTAAAACCGAGTCTGGC`
has the duplicated 9-mer `ACCGAGTCT`; `TTGTTGTTG` trips the identical-triple rule.

**Task 8 — `restriction.rs`.** Sites packed as `u16` (6 bases x 2 bits) with a length,
so matching is a masked compare. A change rechecks the windows ending in `[p, p+2]`.
Unit test: `CATATG` is found for NdeI and not for XhoI.

**Task 9 — `start_site.rs`.** `purine_run[i]` (A/G); for each TTG/CTG/ATG at `i >= 18`,
violation if a run of 5 purines lies fully inside `[i-18, i-3)`. A change rechecks
`[p-2, p+21]`. Unit test: `CATATGGGGGGCTCGAGA` + a TTG start reproduces the Python's hit.

**Task 10 — `splice.rs`.** The 20 donor and 10 acceptor motifs are fixed-length patterns
over character classes. Precompute, per base, a mask with one bit per (pattern,
position); matching is shift-or over `u64` lanes. A change rechecks match starts in
`[p-9, p+2]`. Unit test: the known donors/acceptors from the Python test
(`CGGTGAAT`, `CTGCTACAGC`, ...) are each detected, and the bitap result equals
`reference.rs` on 200 random sequences.

**Task 11 — `hairpin.rs`.** Stems at offsets 0, 10, 20, ...; pattern =
reverse complement with C allowed to match C or T. State: `HashMap<u32,u32>` of 10-mer
counts, per-stem variant list and match count, and a reverse variant->stem index.
Variants capped at 64 with a direct-scan fallback. Unit test: `stem + "TTGTTGTTGA" +
revcomp(stem)` is detected; a random sequence matches `reference.rs`.

---

### Task 12: Search driver

**Files:** Create `src/arrange/search.rs`.

- [ ] **Step 1: Write failing tests**
  - `greedy_seed_beats_random`: greedy left-to-right seeding yields fewer initial
    violations than a random arrangement on GFP (prototype saw 42-58 random).
  - `clears_gfp_to_zero`: GFP at defaults reaches 0 violations, deviation <= 0.01, and
    the translation is unchanged.
  - `swaps_do_not_move_deviation`: with substitutions disabled, deviation and CAI are
    bit-identical before and after the search, for 5 seeds.
  - `budget_is_respected`: with `max_relax = 0.0`, no `Subst` is ever applied.
  - `deterministic_for_a_seed`: same seed, same thread count, identical output; and
    identical output at 1 vs 8 threads.
  - `reports_what_it_cannot_fix`: a poly-Trp protein (no synonyms at all) returns its
    residual violations rather than spinning.
- [ ] **Step 2: Run** `cargo test --lib search` — expect failure.
- [ ] **Step 3: Implement** the loop from spec section 7.3 with the stated defaults
  (`p_sideways` 0.3, `p_noise` 0.05, `stall_limit` 2000, `max_moves` 200n, restarts =
  threads). Candidate peers: `min(n_aa, 16)` sampled. Restarts via rayon, results
  collected by index and reduced deterministically (fewest violations, then lowest
  deviation, then lowest restart index).
- [ ] **Step 4: Run** `cargo test --lib search` — expect PASS.
- [ ] **Step 5: Checkpoint** — clippy clean.

---

### Task 13: CLI, FASTA, reporting

**Files:** Rewrite `src/main.rs`; create `src/metrics.rs`.

- [ ] **Step 1: Write failing tests** — `metrics::deviation` and `metrics::cai` match the
  values asserted in Task 1 for a known sequence; FASTA parsing handles multi-line
  records, `*` stops and lowercase; an unknown amino acid letter is a clear error, not a
  panic.
- [ ] **Step 2: Run** `cargo test` — expect failure.
- [ ] **Step 3: Implement** clap derive with the flags from spec section 8;
  `--cycles`/`--inner-cycles` accepted and warned about; rayon over FASTA records; per
  sequence report deviation, CAI, GC, budget spent, moves, time, and any residual
  violation with its kind and position; exit status 1 if any sequence ends with
  violations.
- [ ] **Step 4: Run** `cargo test && cargo run --release -- --input <gfp.fasta> --output /tmp/out.fasta --local-host-profile tests/data/e_coli_316407.json`
  — expect a clean report and a translatable gene.
- [ ] **Step 5: Checkpoint** — clippy clean.

---

### Task 14: Benchmarks and validation against the baseline

**Files:** Create `benches/optimize.rs`, `scripts/validate.py`.

- [ ] **Step 1:** criterion benchmark over n = 100, 200, 400, 800, 1600, 3200, 6400,
  10000; assert in a test that the fitted exponent between n=200 and n=6400 is <= 1.2.
- [ ] **Step 2:** `scripts/validate.py` reads the Rust output FASTA and, using the
  existing Python package, independently re-checks: translation equals the input
  protein, deviation (Python's own `compare_profiles`), CAI (Biopython), and the count
  of violations. This keeps the Python as an outside auditor of the Rust.
- [ ] **Step 3: Run** both; record the numbers against spec section 11 success criteria.
- [ ] **Step 4:** update the spec's section 11 with measured results.
- [ ] **Step 5: Checkpoint** — full `cargo test`, clippy clean, benchmarks recorded.

---

## Self-review

**Spec coverage:** section 4 layout -> Tasks 1-13 (every file assigned); section 5
composition -> Task 2; section 6 semantics -> Tasks 5-11 (homopolymer and GC changes
are pinned by named unit tests); section 7.1 state -> Task 3; 7.2 indices -> Tasks 5-11;
7.3 search -> Task 12; 7.4 complexity -> Task 14; section 8 CLI -> Task 13; section 9
testing -> Task 4 harness, used by every constraint task, plus Tasks 12-14; section 10
risks -> Task 12 (`reports_what_it_cannot_fix`, `budget_is_respected`) and Task 2
(`gc_adjust_reports_infeasible`).

**Placeholders:** none; every task names exact files, exact test names and exact commands.

**Type consistency:** `SeqState`, `Move`, `CodonIdx`, `AaIdx`, `Violation`, `Kind`,
`Constraint`, `HostTable`, `Composition` are defined once (Tasks 1-4) and used with the
same names and shapes in Tasks 5-14. `touched(&self, mv) -> [Option<u32>; 2]` is the
single scoping primitive every `delta` uses.
