# Codon Harmony: Rust redesign

Date: 2026-10-17
Status: historical design underlying the current engine

The [correctness and performance proposal](2026-10-17-engine-correctness-performance-design.md)
corrects this document's unconditional complexity and validity claims and records
differences between the original plan and the implementation. Its changes remain
proposed; this historical document is retained for context.

## 1. Context

`codon_harmony` reverse-translates a protein and then iteratively cleans the DNA:
resample codons from the host profile, harmonise codon counts to that profile, then
loop six cleanup stages (GC windows, alternate start sites, repeats, homopolymers,
hairpins, restriction sites) `cycles x inner_cycles` times, rank the samples by codon
adaptation index (CAI), and finally remove splice sites from the winner.

### 1.1 Measured baseline

Host: *E. coli* 316407, threshold 0.10, local table. Protein: GFP, 238 aa.
Python 3.8, Biopython 1.76, M-series Mac.

Time, at `cycles=2 inner_cycles=2` (default is `10x10`, i.e. 25x this work):

| protein length | wall time | vs half length |
|---|---|---|
| 100 aa | 0.063 s | |
| 200 aa | 0.154 s | x2.45 |
| 400 aa | 0.436 s | x2.83 |
| 800 aa | 5.59 s | x12.8 |
| 1600 aa | 17.3 s | x3.09 |

100 aa to 1600 aa is 16x the length for 274x the time: an empirical exponent of 2.0.
At default settings the 1600 aa case extrapolates to roughly 7 minutes for one sequence.

Where the time goes (cProfile, `cycles=3 inner_cycles=3`): `remove_local_homopolymers`
32%, `gc_scan` 31%, `remove_repeating_sequences` 16%, `harmonize` 9%, `remove_hairpins` 4%.

Quality, default `10x10`, 12 seeds:

- 5 of 12 seeds produced no gene at all.
- All 7 that succeeded ended at **12.1%-16.3% deviation** from the host profile,
  against a stated `max_relax` of 10%. Not one run met its own bar.
- 12 constraint violations survived in the final sequence (median, scored
  independently).
- CAI median 0.701.

Net deviation caused by each stage over one default run:

| stage | net change in deviation |
|---|---|
| `harmonize_codon_use_with_host` | **-0.886** |
| `remove_start_sites` | +0.389 |
| `remove_local_homopolymers` | +0.187 |
| `remove_hairpins` | +0.149 |
| `remove_restriction_sites` | +0.063 |
| `remove_repeating_sequences` | -0.034 |
| `gc_scan` | -0.001 |

### 1.2 Root causes

1. **Profile drift is never repaired.** `harmonize` is the only restoring force and it
   runs once per sample, *before* the cleanup stages. Every cleanup mutation is
   profile-blind. The `max_relax` gate is also checked before splice removal, so the
   limit is violated in 100% of successful runs.
2. **Random walk instead of directed repair.** `mutate_codon` draws a random synonymous
   codon without checking whether it resolves the violation. `gc_scan` picks a random
   position in the window and keeps the result only when GC happens to move the right
   way. The 30-pass caps and the `10x10` cycle nest exist to compensate.
3. **Stages undo each other.** The pipeline is not a fixpoint; 12 violations remain at
   the end of a default run.
4. **Whole-sequence rescans per window.** Repeats run `re.findall(window, whole_seq)`
   for each of n windows; hairpins regex the whole gene per stem; the `twist_long` GC
   window is `x3`, i.e. three times the gene length, so it spans the entire gene. This
   is the quadratic term.
5. **Exponential splice repair.** `product(*synonymous_codons)` enumerates up to 6^k
   combinations per site, each re-scanning the gene with 30 regexes.
6. **`relax` does not relax the gate.** `compare_profiles` derives the reported
   deviation from `ideal_usage_abs`, which ignores `relax`. So `relax` only skews the
   harmonisation target, making later samples worse, while the gate still compares
   against the unrelaxed deviation.

The existing Rust port in `src/` is a loose translation with its own defects (its repeat
detector searches only within the window rather than the whole gene; its
`compare_profiles` differs from the Python). It is replaced wholesale.

## 2. Goals

- Reduce time complexity from `O(cycles x inner x n^2)` to `O(n)` expected.
- Produce sequences much closer to the host codon profile, with the budget enforced
  against the final sequence.
- Finish with zero constraint violations, or report precisely what could not be fixed.
- Keep space at `O(n)` plus small fixed tables.

Non-goals: changing the host-table format, the FASTA interface, or the biological
meaning of the constraints beyond the two corrections in section 6.

## 3. Key insight

Deviation and CAI are functions of the **multiset** of codons alone. Every sequence
constraint (GC windows, homopolymers, repeats, hairpins, restriction / start / splice
motifs) is a function of the **arrangement** alone.

Therefore a move that exchanges the codons at two positions holding the same amino acid
leaves the multiset unchanged, and so leaves deviation and CAI *exactly* unchanged,
while changing the arrangement freely. The two halves of the problem separate:

- **Composition**: which codons, how many of each. Solved exactly, once.
- **Arrangement**: where each one goes. Searched, with composition frozen.

### 3.1 Evidence

A deliberately naive prototype (random position choice, ten random partners, full rescan
per candidate, 4 s budget) on GFP, rearrangement moves only:

| | current, 10x10 | prototype |
|---|---|---|
| deviation | 13.75% | **0.83%** |
| CAI | 0.701 | **0.738** |
| violations left | 12 | 2 |
| produced a gene | 6/8 seeds | always |

Deviation and CAI came out bit-identical across all three seeds, confirming the
invariance. The two residual violations are an artifact of the prototype's unguided
search, not of the approach.

## 4. Architecture

```
src/main.rs                        CLI, FASTA I/O, rayon across records
src/lib.rs                         public API
src/host.rs                        codon table load, threshold, profile, CAI weights
src/composition.rs                 apportionment + GC and budget adjustment
src/metrics.rs                     deviation, CAI, GC reporting
src/arrange/state.rs               sequence state, Swap / Subst primitives
src/arrange/search.rs              min-conflicts driver, escalation, restarts
src/arrange/index.rs               shared incremental structures
src/arrange/constraints/mod.rs     Constraint trait, registry, violation accounting
src/arrange/constraints/gc.rs
src/arrange/constraints/homopolymer.rs
src/arrange/constraints/repeat.rs
src/arrange/constraints/hairpin.rs
src/arrange/constraints/restriction.rs
src/arrange/constraints/start_site.rs
src/arrange/constraints/splice.rs
```

Every constraint implements one trait, so each can be tested against a brute-force
reference in isolation:

```rust
pub trait Constraint {
    fn rebuild(&mut self, s: &SeqState);
    fn violations(&self) -> usize;
    fn delta(&self, s: &SeqState, mv: Move) -> i32;   // local evaluation, no mutation
    fn apply(&mut self, s: &SeqState, mv: Move);      // commit incremental update
    fn report(&self, s: &SeqState) -> Vec<Violation>; // for the final summary
}
```

`delta` must never scan more than a constant number of positions around the move.

## 5. Composition

Input: per-amino-acid position counts from the protein, host profile `f` (thresholded
and renormalised), budget `B` = `max_relax`, whole-gene GC band.

1. **Apportionment.** For each amino acid `a` with `n_a` positions and synonymous codons
   `c` with frequency `f_c`: take `k_c = floor(f_c * n_a)`, then hand the remaining
   `n_a - sum(k_c)` units to the largest fractional parts (largest-remainder / Hamilton).
   This minimises `sum |k_c - f_c * n_a|` subject to `sum k_c = n_a`, so it is the
   L1-optimal integer allocation.
2. **Deviation metric.** Keep the existing definition for comparability. Both sums
   run over all 64 codons, where `n_a` is the position count of the amino acid that
   codon belongs to:
   `D = sum_c |k_c - round(f_c * n_a)| / (2 * sum_c round(f_c * n_a))`.
3. **Whole-gene GC.** Compute `g` from the multiset. Its band is the one `twist_long`
   carries today, 0.30 to 0.65. While `g` is outside that band,
   apply the unit move (shift one count from codon `x` to synonymous codon `y`) with the
   best `|delta g| / delta D` ratio. If the band cannot be reached, report it rather
   than silently continuing.
4. **Reserve** `B_left = B - D` for the arrangement phase.

The budget is spent on feasibility, not on inflating CAI. There is no `relax` ramp: it
was a way of trading profile match for CAI across samples, and it made later samples
strictly worse.

## 6. Constraints

| constraint | definition | change from today |
|---|---|---|
| GC windows | sliding 20 / 50 / 100 nt windows within their bands | `twist_long` (`x3`, the whole gene) moves into composition, where it is exact |
| Homopolymer | run of more than N consecutive identical bases | **was** the count of the commonest base in a 6-nt window, so `GCGCGC` tripped a G/C limit of 4 |
| Repeats | 3 identical consecutive codons; any 9-mer occurring twice | unchanged |
| Hairpins | 10-nt stem at offsets 0, 10, 20, ... with a reverse-complement partner anywhere, G-T wobble allowed | unchanged |
| Restriction | the chosen enzymes' recognition sites | unchanged |
| Start sites | four contiguous bases of the `AGGAGG` consensus, 4-14 nt ahead of an ATG, GTG or TTG | **changed**, see section 11b. The old rule (5 purines before TTG/CTG/ATG) stays available as `--start-site-rule purine-run` |
| Splice | the 20 donor and 10 acceptor motifs | unchanged |

## 7. Arrangement engine

### 7.1 State

- `codons: Vec<u8>` codon indices, length n.
- `nt: Vec<u8>` 2-bit bases, length 3n, kept in sync.
- `pos_by_aa: [Vec<u32>; 20]` positions grouped by amino acid.
- Moves:
  - `Swap(p, q)` with `aa(p) == aa(q)` and different codons. Multiset-preserving, free.
  - `Subst(p, c)` with `c` synonymous and `f_c > 0`. Changes counts, charged against
    `B_left`.

### 7.2 Incremental indices

All are O(1) per move because a move touches 6 nt and every window length is a constant.

- **GC**: per window set, an array of per-window GC counts. A base change shifts a
  contiguous range of at most `cw` entries (34 for the 100-nt window); update and
  recount those directly. No log factor, `O(n)` space per window set.
- **Homopolymer**: `run[i]` = length of the run of identical bases ending at `i`.
  A change re-walks only the affected maximal run.
- **Repeats**: a 9-mer packs into 18 bits, so a direct `[u32; 262144]` count table (1 MB)
  gives hash-free O(1) lookup, plus a counter of violating 9-mers. A 9-mer
  violates when it occurs at two positions at least 9 apart, which is what the Python's
  non-overlapping `re.findall` count means. A codon change touches the 9 k-mers
  starting in `[p-8, p]`.
- **Hairpins**: `HashMap<u32, u32>` of 10-mer counts (n entries), plus, per stem, its
  wobble variant set and match count, and a reverse index from variant to stem. Variants
  are capped at 64 with a direct-scan fallback for G-rich stems.
- **Restriction**: sites packed as u16 (6 nt x 2 bits); check the 6 windows ending in
  the changed bases.
- **Start sites**: `purine_run[i]`, with xTG positions rechecked in `[p-2, p+21]`.
- **Splice**: bitap (shift-or) over character-class masks for all 30 motifs; no regex.
  A change rechecks match positions in `[p-9, p+2]`.

### 7.3 Search

```
composition -> multiset
seed arrangement greedily: left to right, take from the amino acid's remaining pool
   the codon creating the fewest new local violations      O(n * |synonyms|)
rebuild indices, collect violated windows into a queue

while violations > 0 and moves < max_moves:
    w = a violated window from the queue
    candidates = { Swap(p, q) : p in w, q in K sampled peers with aa(q) == aa(p) }
                 K = min(positions of that amino acid, 16)
    score each by the summed local delta over all constraints
    if best delta < 0            -> apply
    elif best delta == 0         -> apply with probability p_sideways
    else                         -> with probability p_noise apply a random swap,
                                    otherwise try Subst candidates if B_left allows,
                                    cheapest in deviation first
restart with a new seed when stuck, meaning no improvement in stall_limit moves;
restarts run in parallel under rayon
```

Defaults: `p_sideways` 0.3, `p_noise` 0.05, `stall_limit` 2000 moves, `max_moves`
200n, `restarts` equal to the thread count. All are CLI-overridable.

Return the best result: fewest violations, then lowest deviation. Restart results are
collected by index so the output does not depend on thread scheduling.

### 7.4 Complexity

| | now | redesign |
|---|---|---|
| per run | `O(cycles * inner * n^2)` = 100n^2 | **O(n)** expected |
| per move | `O(n)` rescan | O(1) amortised |
| splice repair | `O(sites * 6^k * n)` | O(1) per candidate |
| space | O(n) | O(n) + 1 MB k-mer table per worker |

The prototype needed about 1.5 accepted moves per residue, supporting the `O(n)` move
count.

## 8. CLI

Kept: `--input --output --host --host-threshold --local-host-profile --max-relax
--restriction-enzymes --remove-splice-sites/--no-remove-splice-sites
--remove-start-sites/--no-remove-start-sites --local-homopolymer-threshold --verbose`.
`--start-site-rule` picks between `shine-dalgarno` (the default) and `purine-run`.
`--local-homopolymer-threshold` keeps its name but now means the longest run of
identical bases allowed, per section 6.

Added: `--seed`, `--restarts`, `--threads`, `--max-moves`.

Dropped: `--cycles`, `--inner-cycles`. They describe the old sampling loop and have no
meaning here. Passing either prints a warning and is otherwise ignored.

Per sequence the tool reports deviation, CAI, GC, budget spent, moves, time, and any
violation it could not clear, by type and position. The current tool reports none of
this, which is how 12 surviving violations went unnoticed.

## 9. Testing

1. **Differential tests per constraint.** Each incremental detector must agree with a
   brute-force O(n^2) reference over random sequences. This is the method that took the
   JavaScript port to zero mismatches against Python across 30 sequences.
2. **Incremental equals rebuild.** After each move in a random move sequence, every
   constraint's counters must equal a full `rebuild`. This is where incremental code
   breaks, so it gets the most coverage: 10^5 random moves.
3. **Invariants.** Translation is preserved by every move; `Swap` leaves the codon
   multiset, deviation and CAI bit-identical; `Subst` changes deviation by exactly the
   predicted amount.
4. **Composition optimality.** Apportionment beats or equals `round()`-based allocation
   on L1 distance for random profiles, and always sums to `n_a`.
5. **End to end.** GFP and a stress protein (tandem repeats, poly-Lys, GC-rich runs)
   reach zero violations within budget.
6. **Benchmarks.** Criterion, n = 100 to 10,000, asserting the fitted exponent stays
   near 1.

## 10. Risks

- **Infeasible rearrangements.** Some proteins cannot be cleaned by rearrangement alone
  (an odd number of lysines makes some 6-nt windows unavoidable under the old rule;
  tandem protein repeats can force duplicate 9-mers). Handled by budget escalation and
  then an explicit report. No silent failure.
- **Budget exhaustion.** If violations remain once `B_left` is spent, the sequence is
  still returned, with the remaining violations listed and a non-zero exit status.
- **k-mer table memory.** 1 MB per worker thread; at 8 threads, 8 MB. Acceptable. A
  `HashMap` fallback is available if this ever matters.

## 11. Results as built

Measured on the same machine and host table as section 1, with `scripts/validate.py`
re-deriving every figure from the Rust output using the original Python package.

### Speed

Single thread, one restart, random proteins with a realistic amino acid composition:

| protein | Rust | vs half length |
|---|---|---|
| 1600 aa | 1.05 s | |
| 3200 aa | 2.04 s | x1.94 |
| 6400 aa | 4.55 s | x2.23 |
| 12800 aa | 9.38 s | x2.06 |

A fitted exponent of **1.05** between 1600 and 12800 residues, against the Python's 2.0.
GFP (238 aa) takes 31 ms end to end including process start and table load, against
roughly 7 s for the Python at its default settings: about 225x.

### Quality, GFP, audited by the Python

| | Python, 10x10 | Rust, defaults | Rust, max_relax 0 |
|---|---|---|---|
| deviation from host profile | 0.1375 | **0.0083** | **0.0083** |
| codon adaptation index | 0.701 | 0.738 | 0.738 |
| produced a gene | 6 of 8 seeds | always | always |
| unresolved constraints | 12 | **0** | 2 |

With the Shine-Dalgarno rule in place GFP needs no budget at all, so the default and the
profile-frozen mode land on the same deviation.

`--max-relax 0` is pure rearrangement: no substitution is ever made, so deviation sits
at the apportionment optimum and CAI is identical across every seed and restart. On the
1600-residue protein that mode reaches a deviation of 0.0013.

### What is left unresolved, and why

The engine stops only where no single move helps, which an integration test checks by
trying every swap and every substitution across each remaining violation's window.
Three sources of genuinely unsatisfiable constraints turned up:

1. **Start sites, under the old purine rule.** Glu is GAA/GAG and Lys is AAA/AAG, so
   `EK` is six purines in every encoding. Where that lands 18 bases ahead of an in-frame
   ATG, which Met cannot change, no codon assignment removes it. GFP contains exactly
   this at `N E K R D H M`. Section 11b replaced that rule, and under the default GFP
   now finishes with no violations at all; the case remains as a test of the old rule.
2. **Tandem repeats.** A protein made of repeated copies forces repeated 9-mers: each
   copy would need its own encoding of every three-residue window, and windows holding
   Met or Trp have only one.
3. **Lysine runs.** Four lysines with three AAA in the multiset give six consecutive A's
   however they are arranged. Two coordinated substitutions fix it, which the paired
   lookahead finds when budget allows.

Every one of these is now reported with its position and cause, and the process exits
non-zero, rather than being written out silently.

## 11b. The start-site rule

Section 11 noted that the inherited start-site rule is unsatisfiable surprisingly often.
Measured on an uncleaned 6400-residue gene, it flags 230 sites, and:

- only **27%** of the motifs contain `AGGA`, `GGAG` or `GAGG`, so the other three
  quarters cannot pair with the 16S rRNA and are not ribosome binding sites at all;
- **37%** of the flagged "start codons" are CTG, which is under 1% of real E. coli
  initiations. It is in the list only because the rule borrows NCBI translation table
  1's `start_codons`, which describes the genetic code rather than what E. coli
  initiates at;
- **43 GTG sites** with the same purine context are missed completely, and GTG is about
  14% of real E. coli starts.

So the rule is loose where biology is strict and strict where biology is loose. It also
generates the unsatisfiable cases: `E K` is six purines in every encoding, while a
Shine-Dalgarno core needs adjacent G's, which synonymous choice can almost always break.

`Rule::ShineDalgarno` (CLI `--start-site-rule shine-dalgarno`) asks the biological
question instead: four contiguous bases of the `AGGAGG` consensus, placed 4 to 14 bases
ahead of an ATG, GTG or TTG. The consensus and the start-codon set are standard; the
spacer window is a deliberately wide superset of the 5-11 bases the literature reports,
since this is a hazard filter rather than an expression predictor.

Effect, same proteins and seeds:

| | purine run | Shine-Dalgarno |
|---|---|---|
| GFP, unresolved | 1 (unfixable) | **0** |
| GFP, deviation | 0.0167 | **0.0083** |
| 1600 aa, unresolved start sites | 7 | **1** |
| 1600 aa, deviation | 0.0221 | **0.0122** |
| 6400 aa, unresolved start sites | 49 | **6** |
| 6400 aa, deviation | 0.0581 | **0.0433** |

Deviation improves as a side effect: budget that was being spent chasing motifs that
were never ribosome binding sites stays unspent.

`shine-dalgarno` is the default. `--start-site-rule purine-run` restores the old
behaviour exactly, and both rules are covered by the same differential tests.

## 12. Notes

The repository is not a git repository, so this document is not committed.
The measurement scripts behind section 1 live in the session scratchpad
(`prof.py`, `scaling.py`, `quality.py`, `probe.js`, `compare.js`).
