# RabbitSketch2 metric/workload validation

The completed Chinese report and CSV/PDF/PNG outputs are at
`/home/bigssd/zt/rabbitsketch2_benchmark/analysis/regimes_20261001/report.md`.
This is the follow-up to the frozen ordinary-Jaccard screening experiment.
MinHash means FastKMV. KSSD is included; OrderMinHash is excluded by request.

Read [PROTOCOL.md](PROTOCOL.md) and [METRIC_AUDIT.md](METRIC_AUDIT.md) before
interpreting the columns. SetSketch requires the isolated experimental patch
from `../screening/patches/setsketch_partitioned_joint.patch`.

The batch adapter is compiled into the existing RabbitSketch `examples/main.cpp`
entry point via `../screening/CMakeLists.txt`. `--screening-manifest` detects
`#regimes-v1` manifests; legacy screening manifests retain their old path.
Each benchmark is one process with 64 OpenMP threads. The shared flock prevents
multiple benchmark processes from running concurrently. Inputs are read first,
each sketch is constructed once, mutable query caches are materialized, and
then immutable sketches are shared across distinct pair queries. Main matrices
have 256/1024/4096/10000 different assembly accessions, three complete passes,
and no diagonal. The compact supplement has 1024/10000 objects.

## Frozen runs used by the final analysis

All paths below are under `/home/bigssd/zt/rabbitsketch2_benchmark/runs/`:

| Directory | Final contribution |
| --- | --- |
| `regimes_20261001` | ANI, cardinality, workload, weighted, FastANI and 74-pair NUCmer audit |
| `regimes_supplement_20261001` | KSSD dr=4 ANI/workload; its first containment measurements are exploratory |
| `containment_holdout_20261001` | All final containment: donor and reference splits are disjoint |
| `compressed_validation_20261001` | Compact BinDash ANI; small MinHash/BinDash actual all-pairs |

The initial containment and CRLF preflight are retained for audit and excluded
from final summaries. Each run has its own `suite_manifest.json`, input SHA256,
binary and source snapshots, commands, raw JSONL and batch completion hashes.
The runner is intentionally strict about identity changes. Do not use a later
working script to overwrite a frozen run; the snapshot records the exact code
used there. The working scripts are specific to this workspace/data inventory.

## Regenerate analysis from completed runs

Use the existing benchmark venv (system pandas is not the pinned environment):

```bash
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/rabbitsketch_regimes_mpl
/home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/regimes/analyze_suite.py
/home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/regimes/analyze_alignment.py
/home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/regimes/paired_intervals.py
/home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/regimes/write_report.py
/home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/regimes/finalize_delivery.py
```

`analyze_suite.py` scores all final containment from the corrected run only.
It refuses incomplete phases. Figures have standalone PDF and PNG versions.
FastANI and NUCmer references remain separate. Exact set-to-ANI conversion is
a model baseline, never an independent alignment truth.
The finalizer verifies all batch-output hashes and freezes the report, plots,
tables, source snapshots and run manifests in the analysis delivery directory.

## Advisor prototype

```bash
/home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/regimes/advise.py \
  --family substitution --ani-bin '99–100' --pair-budget 8192 \
  --target-error 0.001 --queries-per-sketch 10000
```

The budget is **pair-total query payload**. Error is an ANI fraction (0.001 =
0.1 percentage point). Reuse is incident comparisons per sketch, so amortized
cost uses half that number; `T=N*B+P*t` and `R=2P/N`. The calibration choice is
checked on a disjoint holdout. Failure remains a failure, not a silent fallback
chosen using the holdout. Cost is effective 64-thread mixed-bacterial matrix
throughput, excluding parsing/output, and is not single-pair latency. This
prototype handles the tested ANI regimes only; it is not production API wiring.

## Checks

`test_regimes.py` covers cached multi-object traversal, FASTA/CRLF boundaries,
identity invariants, BinDash+HLL composition, weighted scale invariance, and an
independent quadratic probability-Jaccard calculation. Together with the
screening/native-input regressions, 62 tests pass. The existing partitioned
SetSketch model test remains available in the validated CMake build.
