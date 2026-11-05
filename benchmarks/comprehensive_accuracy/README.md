# Comprehensive accuracy and scenario selection

This experiment expands the two-algorithm 2026-10-12 follow-up into a common
seven-algorithm grid and a fuller real/task validation. Runtime comparisons,
scalar/SIMD experiments and performance rankings are excluded. The native source
is the frozen validated checkout with the isolated SetSketch correction; the
production estimator is unchanged.

The registered [protocol](PROTOCOL.md) and [config](config.json) were written on
2026-10-30. Native preflight passed before the user paused the task; the full
suite resumed on 2026-10-31. The workspace keeps the registration date in its
name: `/home/bigssd/zt/rabbitsketch2_benchmark/runs/comprehensive_accuracy_20261030`.

| Panel | Coverage | Independence and limits |
| --- | --- | --- |
| Controlled canonical sets | 4 k values × 3 larger cardinalities × 7 size ratios × 10 small-set containments × 128 inputs; 4 budgets; all 7 algorithms | 840 nominal regimes, 832 distinct realized count conditions after rounding; 107,520 independent generated inputs with 8 balanced hash seeds |
| Natural genomes | Same 900 genomes / 180 species / 60 genera; all within-genus pairs plus frozen cross-genus controls; 6,426 pairs; k=20/24/32 | Added pairs/configurations, not new organisms; genus clusters retained; retrospective reuse of the existing inventory |
| Substitution/indel | Same 320 controlled pairs; full configuration sweep at k=24 | Reference clusters; mutation-model and sampling errors remain distinct |
| Directed containment | Corrected 700 pairs at k=24; all declared parameters | Original reference/donor split preserved; unavailable samples remain failures |
| Weighted probability Jaccard | Existing synthetic abundances and official CAMI profiles; m=32/128/512/2048/8192; 8 hash seeds | Separate estimand; no claim of a best weighted method from one eligible backend |
| Downstream | Species nearest neighbor, exact-J neighbor recovery, single-linkage ARI | 30 calibration / 30 validation genera; clustering threshold chosen only on calibration; fractional credit for ties |

The synthetic plan contains 3,118,080 native configuration/pair measurements;
the real/task plan has 505 native/exact jobs. These are repeated measurements,
not millions of independent biological samples. `status.py` reports completed
archives separately from partial active output.

```bash
# Dependencies: ../screening/requirements.txt
python benchmarks/comprehensive_accuracy/build.py
python benchmarks/comprehensive_accuracy/run.py --stage plan
python benchmarks/comprehensive_accuracy/run.py --stage preflight
python benchmarks/comprehensive_accuracy/run.py --stage synthetic
python benchmarks/comprehensive_accuracy/run.py --stage real

OPENBLAS_NUM_THREADS=1 python benchmarks/comprehensive_accuracy/analyze_synthetic.py
python benchmarks/comprehensive_accuracy/constraint_profiles.py
OPENBLAS_NUM_THREADS=1 python benchmarks/comprehensive_accuracy/analyze_real.py
OPENBLAS_NUM_THREADS=1 python benchmarks/comprehensive_accuracy/downstream.py
python benchmarks/comprehensive_accuracy/figures.py
python -m pytest -q benchmarks/comprehensive_accuracy/test_analysis.py
python benchmarks/comprehensive_accuracy/report.py
python benchmarks/comprehensive_accuracy/verify_publish.py
python benchmarks/comprehensive_accuracy/status.py
```

`analyze_when_ready.py` is an optional incremental analysis driver while native
jobs run. There is one native process at a time, using the existing main
`--screening-manifest ... -t 64` entry point and the shared workspace lock.
Building once against the pinned source, checking identity on resume, and
keeping snapshots protect the native experiment from analysis-code changes.
The new adapter generates each controlled input once and passes it through
every backend's native DNA update interface. It does not replace the estimators.
The independent preflight enumerates exact sets for 24 controlled inputs and
checks all seven algorithms against the earlier adapter on matching input.

Selection uses calibration only. Fixed-size configurations are declared in
advance; FracMinHash has specified expected budget headroom; KSSD resolution
uses the densest tested payload fitting all calibration inputs. Real variable
parameters are fixed across the workload after calibration. Validation budget
violations and missing outputs are reported, never repaired using test data.

The synthetic primary error families are Jaccard absolute error, bounded
directed-containment absolute error and mean per-object cardinality relative
error. The latter averages the two objects within one independent input pair;
per-object cardinality tables separately retain true object size and bytes.
KSSD cardinality/containment are diagnostic-only and excluded as public metric
candidates; unweighted ProbMinHash lacks those outputs. BinDash cardinality
saturation is unavailable, not a zero error. Direct native containment is also
archived so it is not confused with the unified derived estimate.

Error surfaces, candidate sets and validity masks serve different purposes.
A supported unique winner must beat every eligible competitor with a paired
simultaneous interval excluding zero and both >=10% and >=0.0001 MAE improvement.
Intervals share resampling across budgets within each regime/metric; they do
not adjust across the entire 840-regime grid. Real and downstream cluster
intervals are explicitly pointwise. Descriptive minima and calibration nominees
are not automatically statistically supported winners.

The [analysis addendum](ANALYSIS_ADDENDUM.md) explicitly records a diagnostic
added after early results: compare MAE winners with constant-zero/one baselines
and publish illustrative absolute/relative precision profiles. These profiles
are empirical operating conditions, not registered confirmatory endpoints or
error guarantees. The original grid and inferential results remain unchanged.

Full native outputs, per-input analysis and frozen genome inventory remain in
the external workspace. The scripts require that inventory and do not download
it. Large analyses write to `analysis/comprehensive_accuracy_20261030`; compact
publication tables/figures are produced only after completion and verification.

The publication audit checks all 517 native output archives, reconciles the
registered observation counts, recomputes synthetic MAE/P95 from per-input
rows, reruns analysis invariants, and fingerprints the compact evidence.
Figures cover synthetic candidate sets and precision requirements, natural
genomes, controlled containment, mixed-mutation model error, weighted profiles
and downstream tasks. The three-dimensional precision/memory/runtime study
remains outside this accuracy-only run; the descriptive frontier here contains
only error and query payload.
