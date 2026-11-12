# Full-233 inverted-index construction comparison

This benchmark compares local-map insertion + 64-shard merging with the
1024-bucket scatter/sort builder, on identical full-233 FastKMV and weighted
ProbMinHash keys. All jobs run sequentially at 64 threads, with K/M=1024,
k=21, seed=42 and the CLI's strict distance cutoff 0.05.

Each algorithm first calls the actual CLI's FASTX/sketch constructor once and
saves a versioned key snapshot. Each index method then runs in its own process,
loads the same snapshot, constructs the index, and runs the actual native
similarity/output kernel. FastKMV restores and finalizes its native sketches
before index timing; ProbMinHash uses its native positional winner counts.

The index timer includes **all** local-map insertion, shard merging and CSR
flattening for the old method. For the new method it includes counting,
scattering, bucket sorting, singleton removal and CSR construction. Snapshot
loading and native-object restoration are outside this timer. Query/output
has a separate timer. Consequently these results isolate the index change;
they are not measurements of two complete FASTX-to-output CLI executions.

The process peak up to index completion includes the loaded keys and, for
FastKMV, native sketch objects. Full-process peak and wall time are recorded
separately. Inputs, source hashes, binary hashes, commands, ordinary-cache
policy and one-run-per-method order are archived under
`/home/bigssd/zt/rabbitsketch2_benchmark/runs/index_build_233_20261112`.

Before full runs, a 128-input real panel must match the current production
CLI. The parallel full-output fingerprint must match the previous serial
checker on real outputs and artificial line-boundary cases. Full outputs
must match between index methods and the preserved full-233 reference by row
count, byte count, and normalized XXH3-128 sum/XOR fingerprints. Fingerprints
are integrity checks, not mathematical proofs of equality. Fully verified
large TSVs are removed after samples and fingerprints are saved; key snapshots
are retained. Failed comparisons retain the full outputs for investigation.

Build and run:

```bash
cmake -S benchmarks/index_build_233 \
  -B /home/bigssd/zt/rabbitsketch2_benchmark/build_index_build_233_20261112 \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_C_COMPILER=/usr/bin/gcc
cmake --build /home/bigssd/zt/rabbitsketch2_benchmark/build_index_build_233_20261112 \
  --parallel 8 --target index_build_233 output_digest
python3 -u benchmarks/index_build_233/run.py
```

The controller holds the same benchmark lock used by the paused thread-scaling
experiment. It refuses to overwrite an already registered measurement plan.
