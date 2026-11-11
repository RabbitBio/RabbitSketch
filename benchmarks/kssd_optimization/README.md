# KSSD construction optimization, 2026-11-11

This experiment compares the preserved implementation with the optimized KSSD
implementation. It is separate from the stopped 233 thread-scaling experiment;
its completed 47 points remain intact.

## Full-data comparison

- Input: the original 233 list, 208,753 files, 945,020,243,458 input bytes.
- Both programs use native compilation with identical compiler options.
- 64 threads, CPU affinity 0–63, fixed OpenMP placement.
- Effective k=20, half_subk=6, drlevel=3, distance cutoff 0.05, inverted index.
- Complete CLI wall time includes reading, sketching, index construction,
  comparison and TSV output. Output validation is excluded from timing.
- One completed run per version, serial order after/before. The user stopped
  the second pair; its incomplete attempt is excluded. Ordinary filesystem
  caching; no cache eviction or concurrent benchmark jobs.
- Compare row counts and order-independent 128-bit sum/XOR fingerprints of both
  (genome pair, distance) and all output columns. Also compare pair/distance
  fingerprints with the earlier validated full-data result. These are integrity
  fingerprints, not a mathematical proof of equality.
- Retain raw logs, command lines, `/usr/bin/time -v`, samples and digests. Delete
  only this experiment's large TSV files, after successful validation.

The user's subsequent request switched directly to full-data 64-thread runs;
the optional 1/8/64-thread subset CLI grid in `run.py cli` was not run.

## Construction microbenchmarks

`bench.cpp` loads identical FASTX records before timing. Each case has one warmup
and three measured repetitions. Sorted per-genome hash dumps must match by SHA256
across versions, lookup methods and deduplication methods. The unified builder
case compares its complete pairwise query matrix because its state is type-erased.

Cases include 128 real genomes, the first genome split into 1000-base records,
and eight independently built file sketches. Parameter initialization is recorded
separately except in the unified builder test, where its per-builder cost is
part of the problem under study. Microbenchmark peak RSS includes preloaded input;
it must not be presented as standalone sketch memory.

Default lookup uses an exact admission bitmap and compact rank map. Direct
permutation lookup is retained for dense drlevel=0. The map/direct/bitmap paths
are tested against the same permutation. Buffered append/sort was also measured;
the CLI defaults to the faster reusable hash set on the measured real panel.
The optional bounded buffer spills to a set on large or repetitive inputs.

## API compatibility

`Kssd::update(sequence, length)` accumulates across records and materializes sorted
hashes at getters, queries, `toLite()`, projection, or `finalize()`. The historical
single-argument `update(sequence)` continues to materialize public hash arrays
before returning. The builder and Python binding use the length-aware overload.
Generated parameters are shared by configuration; custom shuffle tables retain
their original behavior. The permutation, reduced tuples and estimators are unchanged.

## Reproduction

The before snapshot and raw artifacts are under
`/home/bigssd/zt/rabbitsketch2_benchmark/runs/kssd_optimization_20261111`.
`before_manifest.json` records the original source hashes and Git revision;
`before.patch` preserves the pre-existing workspace modifications. Versioned
source snapshots distinguish the initial buffered candidate from the final CLI.

```bash
cmake -S benchmarks/kssd_optimization -B <build-after>
cmake --build <build-after> --target rabbitsketch_main kssd_bench -j 8
# Configure another build with -DRABBITSKETCH_SOURCE=<run>/before
# and -DRABBITSKETCH_BUILD_TESTS=OFF for the preserved reference.
python3 benchmarks/kssd_optimization/validate.py
python3 benchmarks/kssd_optimization/run.py micro
python3 benchmarks/kssd_optimization/run.py full
```

The scripts use explicit local data/build paths. Review those paths before
reproducing the experiment on another machine. Do not overwrite an existing run
directory when preserving a previous experiment. `run.py full` now runs one
paired comparison, as requested by the user. The original two-pair plan is
archived with the stop record.
