# Observable scenario selection and independent validation

This is a precision-oriented research extension, not a runtime benchmark or a
production automatic-selection API. The [protocol](PROTOCOL.md),
[pre-validation amendments](ANALYSIS_ADDENDUM.md), [native grid](config.json) and
[four unseen validation conditions](validation_extra_conditions.json) define the
experiment. The Oct30 evidence remains unchanged.

Development uses the previous 900 genomes / 60 genera. Validation uses 450 new
genomes / 90 species / 30 genera disjoint from that development panel, excluding
all prior accession aliases, BioSamples and strand/order-invariant sequence
fingerprints. All 30 genera occurred in older broad panels: see the explicit
[historical-scope correction](REPORTING_NOTES.md). There are ten development and
14 validation conditions, k=20/24/32, four pair-payload budgets and seven backends.
Independent validation means **30 biological genus clusters**, not millions of
independent observations. Read/fragment/incomplete-assembly inputs are controlled
derivatives of these real genomes, not new clinical or real sequencing cohorts.

The questions are separate:

* How accurately can the configuration estimate Jaccard, query containment or
  query cardinality under a final sketch-payload limit?
* Can it recover exact-Jaccard neighbors, distinguish source species, or cluster
  nine query genomes/derivatives using a six-reference known-genus candidate bank?
* Does a frozen observable rule improve on the strongest development-fixed
  configuration per task/input kind/k/budget, including its failures?
* Does the choice transfer to previously unseen coverage, length, error and
  completeness conditions, and when should its prediction be treated as unsupported?

Use FastKMV, FracMinHash, HLL, **experimentally corrected SetSketch**, KSSD,
unweighted ProbMinHash and BinDash. Weighted probability Jaccard and order
similarity are separate estimands. Legacy MinHash is the traditional implementation
baseline. None of these experiments ranks speed, SIMD, construction memory or RSS.

## Reproduction order

The scripts use `/home/bigssd/zt/rabbitsketch2_benchmark` for large data/raw files.
The native executable is the Oct30 frozen accuracy adapter, with its existing
isolated SetSketch correction. `../comprehensive_accuracy/build.py` documents its
build identity. Local NCBI assemblies and taxonomy files are needed for regeneration;
the published accession/provenance manifests identify the exact inputs.

```bash
python benchmarks/scenario_selection/prepare.py
python benchmarks/scenario_selection/extend_validation.py
python benchmarks/scenario_selection/audit_history.py
python benchmarks/scenario_selection/run.py --stage development
python benchmarks/scenario_selection/measure.py --split development --stage features
python benchmarks/scenario_selection/measure.py --split development --stage scores
python benchmarks/scenario_selection/advisor.py --stage fit
python benchmarks/scenario_selection/run.py --stage validation-pilot
python benchmarks/scenario_selection/measure.py --split validation --stage features
python benchmarks/scenario_selection/advisor.py --stage freeze-decisions
python benchmarks/scenario_selection/run.py --stage validation-outcomes
python benchmarks/scenario_selection/measure.py --split validation --stage scores
OPENBLAS_NUM_THREADS=1 python benchmarks/scenario_selection/evaluate.py
python benchmarks/scenario_selection/figures.py
python -m pytest -q benchmarks/scenario_selection/test_selection.py
python benchmarks/scenario_selection/report.py
python benchmarks/scenario_selection/verify_publish.py
```

`continue_pipeline.py` can continue after calibration, preserving the same
dependency order. Native processes share the existing benchmark lock and use
64 OpenMP threads. Resumption checks completed raw hashes and registered manifest/
binary identities. Validation outcomes cannot start without frozen model artifacts
and precomputed choices. Scoring can overlap a native runner using `--wait`; only
completed archives are read. Errors stop the pipeline instead of producing partial
success reports.

## Advisor inputs and costs

The rule chooses one compatible backend/parameter for a workload before the final
sketches/index are built. It cannot switch algorithms inside an existing index or
compare sketches from different backends.

Allowed inputs are input kind, k, base/record-length statistics, and a K=32 FastKMV
pilot using seed 919. The pilot provides estimated cardinalities, size ratios,
similarities, zero rates and top-score margins. True Jaccard/containment/counts,
species/genus identifiers, actual mutation rates, contamination and generation
coverage are not model features. Training group identifiers only enforce distinct
genus neighbors. The rule uses eight nearest distinct development genera, fixed
median/IQR scaling and no hyperparameter search.

Pilot payload is at most 256 bytes per object, or 3,840 bytes for the 9+6 object
workload. It is released before building the final sketches. This is additional
temporary query-state payload, not RSS; retaining both pilot and final sketches
requires adding their payloads. The pilot's input scan and building time are not
claimed to be free or beneficial. Large-bank/global species search is not validated.

The metadata-only variant is an ablation. The fixed baseline and seven individual
backend baselines are calibrated on development data only. A variable-payload
guard uses 1.5 times pilot cardinality (or input-base upper bounds for metadata).
This guard is heuristic; actual memory-budget violations remain failures.

The portable model can score a JSON object containing the explicitly named
observable features in `measure.PILOT_FEATURES`:

```bash
python benchmarks/scenario_selection/advisor.py --stage recommend \
  --model benchmarks/results/scenario_selection_20261105/model \
  --features-json benchmarks/results/scenario_selection_20261105/example_features.json \
  --k 24 --kind assembly --metric jaccard --budget 8192
```

Returned neighboring MAE/P95 values are empirical predictions, not confidence
bounds or accuracy guarantees. `out_of_support_candidate` is explicitly unvalidated.
Only the registered k/budget grid is supported. Capability mismatch is rejected.
The library's production API is not changed by this experiment.

## Reading the evidence

Unavailable/over-budget outputs have loss 1. Cardinality loss is clipped relative
error for selection; reported conditional MAE/P95 remain unclipped. Neighbor and
species losses are one minus fractional tie credit; clustering loss is
`(1 - ARI) / 2`. These objectives must not be averaged into a single score.

Precision attainment requires both error limits, at least 99% valid in-budget
outputs and at most 1% budget violations. Report recommendation coverage and
false-feasible rate together. Constant-zero/one and positive-truth relative errors
diagnose apparently good absolute errors with little information. Exact k-mer
task scores are an empirical reference, not aligned ANI or a biological upper bound.

The offline oracle minimizes the same failure-penalized loss over all measured
configurations. A second feasible-oracle column requires complete payload fit and
99% availability. Neither is used to make choices. Held-out candidate maps are
descriptive results after validation; they are **not** the pre-frozen selection rule.

Intervals average seeds within workload/genus before genus resampling. Primary
advisor-vs-fixed intervals are simultaneous over 72 task/k/budget contrasts for
each advisor. Scenario intervals are pointwise; algorithm candidate intervals are
simultaneous only within a scenario cell. They are not global guarantees over all
scenario cells. Out-of-support inputs, unsuccessful transfer and negative downstream
results must remain visible.
