# Accuracy and task conditions, 2026-10-12

Registered before running the new native measurements. The question is which
algorithm/configuration to consider for a specified metric and input regime.
Construction time, comparison time, SIMD/scalar ablations and performance
rankings are outside this experiment. The existing adapter emits incidental
timings; none enter the analysis or recommendations.

## New experiment

The earlier screening and independent-genome results motivate a **targeted**
comparison of BinDash (16-bit signatures) and SetSketch with the isolated
partitioned-joint estimator correction. This is not a new seven-way competition
and cannot eliminate the other candidates. Do not change the production
SetSketch implementation or the frozen source/build used in earlier runs.

The complete grid is in `config.json`: k=24; larger exact set cardinality 100,000
or 1,000,000; larger/smaller cardinality ratios 1, 2 and 4; 15 requested Jaccard
values between 0 and 1. Reject impossible cells with J > 1/ratio. There are 60
feasible regimes, 30 independently generated data replicates per regime, six
configurations, and 10,800 native measurements. Hash seed is fixed at 42;
independent data seed 20261012 differs from the original screen. These are
exact canonical k-mer sets encoded as DNA carriers, not simulated whole genomes.

Compare **A+B query payload** budgets 2, 8 and 32 KiB. Matched configurations
are BinDash bins 512/2048/8192 and SetSketch precision 10/12/14. This payload is
not total resident or construction memory. Pair every comparison on the same
input digests, truth and replicate. Any missing, failed, duplicated, nonfinite
or over-budget row fails the analysis, rather than silently changing a cohort.
J=0 and J=1 are controls; a zero error in both algorithms is not a win.

Primary outcome is absolute Jaccard error. Report MAE, RMSE, bias, P95 and zero
estimate rate; report relative MAE only for J>0. For every regime and budget,
estimate d = MAE(SetSketch-fixed) - MAE(BinDash). Use 10,000 paired bootstrap
resamples of the 30 independent inputs. Reuse resampling indices across budgets
within a regime. Alongside descriptive pointwise 95% intervals, use a centered
bootstrap maximum standardized deviation to construct approximate simultaneous
95% intervals across all nondegenerate regime/budget contrasts. A supported
preference requires this simultaneous interval to exclude zero **and** at least
10% lower observed MAE. Zero-variance contrasts remain unresolved unless both
algorithms have zero error, which is reported separately. Resampling is
conditional on these implementations, input generator and fixed hash seed.

Do not fit a universal J threshold to these outcomes. Report supported cells
and unresolved cells; any apparent transition may depend on cardinality, ratio
and budget. Exact J is an explanatory variable here, not a free feature an
online selector can observe. A deployable selection policy would require
observable features, calibration and independent validation.

## Existing evidence, reinterpreted without speed

Keep existing cohorts distinct from the new experiment:

* Seven-algorithm synthetic screen: merge original native rows with the frozen
  SetSketch correction and FracMinHash budget extension; retain KSSD. Use exact
  truth, the largest tested payload fitting the budget (without consulting
  error), all expected replicates, and the same maximum-byte eligibility rule.
* Independent real genomes: use the existing 900-genome/1,967-pair validation
  and its genus-cluster intervals. No new independent biological validation is
  implied by reanalysis. Missing parameter coverage is not algorithm failure.
* Directed containment: use the corrected held-out cohort, preserve unavailable
  outputs, and distinguish native containment, derived containment and composite
  BinDash+HLL. KSSD derived containment is diagnostic, not a supported API feature.
* Cardinality: compare the same register payload and input size. Reuse of a
  HLL/SetSketch-only panel cannot establish an all-algorithm optimum.
* Weighted task: probability Jaccard and min/max weighted Jaccard are different
  estimands. The existing weighted experiment validates ProbMinHash's own error;
  one weighted candidate cannot establish superiority over alternatives.

The final guide must distinguish an experimentally supported preference, a
metric/API eligibility decision, and an unproven hypothesis. Do not invent a
scenario for every backend or import old speed-based recommendations. ANI
conversion and clustering/near-species decisions require their own reference
and downstream validation; Jaccard gains alone do not prove them.

## Reproduction

Use the frozen validated binary and its source manifest, with one native
process at a time and 64 OpenMP threads. Acquire the existing workspace-wide
benchmark lock. Keep new raw outputs under a new `runs/selection_accuracy_20261012`
directory; preserve the config, protocol, runner/adapter snapshots and hashes
before launch. The runner registers a plan and fingerprints the binary, source,
inputs and execution environment. Source scripts and accuracy tables/plots are
published separately from the large raw workspace.

The legacy runner requires several performance-analysis config fields. They
are schema placeholders in this config and are not used for selection here.
