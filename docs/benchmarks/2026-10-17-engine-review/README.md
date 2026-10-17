# Evidence for the correctness and performance design

These are exploratory measurements, not a production performance guarantee. The
[technical design](../../superpowers/specs/2026-10-17-engine-correctness-performance-design.md)
explains their interpretation, limitations and planned acceptance criteria.

## Environment and inputs

- Baseline: `06061f40d4822b35795743900b5438235abbf161`.
- macOS 26.6.2, arm64; rustc 1.73.0; repository release profile (opt-level 3, LTO).
- Engine: one Rayon thread, one restart, seeds 1/2/3, `max_relax=0.01`.
- Host: `tests/data/e_coli_316407.json`, threshold 0.10; remaining Rust options default.
- GFP is the repository's 238-residue fixture. Synthetic proteins sample the 20 amino
  acids uniformly with seed 2026, independently restarted at each length. These are
  not the README's weighted-composition scaling inputs.
- Kernel: 10,000 evaluations of the same synonymous substitution in `GCT` repeats,
  with `black_box`, excluding initialization and accepted-move application.

`comparison_micro.txt` contains the baseline/prototype kernel comparison, a 20,000-move
differential check, hairpin oracle discrepancy, report-size measurements and the positive
limit initialization discrepancy. The probe's hairpin oracle scans all pairs independently
of the optimized detector, using full 10-base stems and its directional compatibility rule.

The four CSVs contain two sequential rounds, each baseline then patched. First-round
apparent improvements disappeared on repetition. These were not randomized, statistically
powered benchmarks. Most engine cases ended with violations, so these are termination
times rather than validated time to success. The DNA hash is a reproducibility check,
not a proof of validity or a collision-free equality check. Peak RSS was not measured.

Both baseline and patched suites passed 72 unit tests and 8 integration tests. The patch
only changes repeat proposal scoring; it preserves baseline correctness limitations.
All probe/patch code is covered by the repository's MIT license.

## Reproduce in a disposable checkout

The probe's `baseline` implementation comes from the crate being built. Build and save
that executable **before** applying the patch. Applying the patch first would compare
the optimized repeat scorer against itself.

From this design branch, run the following using a new temporary directory. This leaves
the current checkout untouched. If dependencies are not cached, omit `--offline`.

```sh
review_root=$(pwd)
review_tmp=$(mktemp -d)
git clone --no-hardlinks "$review_root" "$review_tmp/codon_flush"
cd "$review_tmp/codon_flush"
git checkout --detach 06061f40d4822b35795743900b5438235abbf161
mkdir -p examples/probe_support
cp "$review_root/docs/benchmarks/2026-10-17-engine-review/comparison_probe.rs" examples/
cp "$review_root/docs/benchmarks/2026-10-17-engine-review/probe_support/fast_repeat.rs" examples/probe_support/
cargo test --release --offline
cargo build --release --offline --example comparison_probe
cp target/release/examples/comparison_probe "$review_tmp/baseline"
"$review_tmp/baseline" > "$review_tmp/micro.txt"
"$review_tmp/baseline" engine > "$review_tmp/baseline_engine.csv"
git apply "$review_root/docs/benchmarks/2026-10-17-engine-review/repeat_delta.patch"
cargo test --release --offline
cargo build --release --offline --example comparison_probe
cp target/release/examples/comparison_probe "$review_tmp/patched"
"$review_tmp/patched" engine > "$review_tmp/patched_engine.csv"
"$review_tmp/baseline" engine > "$review_tmp/baseline_engine_round2.csv"
"$review_tmp/patched" engine > "$review_tmp/patched_engine_round2.csv"
```

The probe needs the repository root as its current directory to find the host table.
The reproduction creates result files under `review_tmp`; their timings will differ
across machines and runs. The fast scorer source is a copy of the baseline repeat
module with the documented delta change and its test module omitted. It is kept in
this evidence directory for side-by-side compilation, not linked into the engine.

For the CLI positive-limit reproduction, write `src/testdata.rs`'s GFP sequence to a
protein FASTA and run the baseline (or repeat-only patched) engine with:

```sh
cargo build --release --offline
./target/release/codon_harmony \
  --input gfp.fasta --output gfp_dna.fasta \
  --host-table tests/data/e_coli_316407.json \
  --threads 1 --restarts 1 --seed 0 --max-relax 0.001 --quiet
```

The measured repeat-only patched run exited 0 and printed:

```text
GFP    deviation 0.0083    CAI 0.738    GC 0.462    violations 0
```

The positive threshold is intentionally below the achievable rounded-profile minimum;
the expected corrected behavior is a profile-limit infeasibility result, not an attempt
to make that impossible threshold achievable. `--max-relax 0` has a separate documented
freeze-composition meaning and is not the failing positive-limit case.
