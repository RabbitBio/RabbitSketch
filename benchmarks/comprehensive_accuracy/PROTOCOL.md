# Comprehensive accuracy and selection experiment — 2026-10-30

Written before launching new measurements. This replaces neither the frozen
2026-10-01 runs nor the 2026-10-12 targeted experiment. No runtime comparison,
ISA ablation, performance ranking or production estimator modification is part
of this study. The task is to identify conditional candidate sets, explain
failures, and test whether accuracy differences matter to downstream decisions.

## Registered synthetic grid

Seven candidates: FastKMV (`minhash`), FracMinHash, HLL, SetSketch-fixed,
KSSD, unweighted ProbMinHash, BinDash (16-bit). SetSketch uses the already
validated isolated partitioned-joint correction, not the production estimator.
Legacy MinHash and OrderMinHash remain outside the comparative scope.

The complete grid in `config.json` crosses four k values (20/24/28/32), three
larger cardinalities (10^4/10^5/10^6), seven larger/smaller ratios (1 to 500),
ten target small-set containments C (0 to 1), and 128 independent input
replicates. Use the nearest feasible integer intersection C*|small| and record
realized C and J=I/(|small|+|large|-I). Small-set integer rounding is explicit;
different nominal C cells that round to the same truth are not independent
overlap conditions. This parameterization avoids impossible J > 1/ratio cells.
There are 840 nominal regimes and 107,520 input pairs. Data seeds are derived
from the registered grid and replicate. Eight hash seeds cycle evenly across
the independent replicates, with identical input/seed pairing across backends.
KSSD's hash seed does not change its reduction table; never count artificial
seed duplicates as independent evidence for it.

At A+B query payload budgets 512/2048/8192/32768 bytes, fixed-size configurations
are K=m=B/16 (FastKMV/ProbMinHash), bins=B/4 (BinDash), and p=log2(B/2)
(HLL/SetSketch). FracMinHash scaled is ceil(8*(n_a+n_b)/(0.8*B)), determined
without looking at errors, providing explicit expected headroom; observed bytes
and violations remain outcomes. All five KSSD dr/half_subk pairs are measured;
budget choice uses the densest configuration fitting every calibration input.
KSSD's coarse resolution ladder and integer width are reported, not hidden.

Each independent input is generated once, then consumed by every native backend
and configuration through its DNA update interface. No substitution of a
generic estimator, reduced hash-width accounting, or exact-cardinality oracle
in estimates. Report Jaccard, native cardinalities, native/derived containment,
payload, empty samples and input digests. KSSD cardinality/containment are
diagnostics, not supported public metrics; ProbMinHash-unweighted has no
cardinality output. Containment uses each backend's cardinality estimate and
the same common-intersection cap as the public unified query, with direct
FracMinHash containment also recorded. Cardinality uses per-object payload.

The new adapter reuses the existing canonical k-mer generator and calls the
frozen corrected native library through the existing main `-t 64` entry point.
It omits repeated timing loops. A small native preflight independently enumerates
the generated truth and compares estimates/payloads with the frozen adapter.
Do not interpret this as a new optimized implementation or performance result.

## Broader real and task panels

Reuse the frozen 900-genome, 180-species, 60-genus inventory, but measure all
seven candidates and complete the parameter coverage at k=20/24/32. Extend
the pair set to every within-genus pair plus the previously frozen cross-genus
controls. This adds comparisons and configurations, not independent organisms.
Use seeds 42/101/202, one KSSD reduction. FracMinHash and KSSD sweeps are fully
declared in the config; no parameters are chosen based on test error. Preserve
all size/availability failures and independently computed exact canonical truth.

Also measure the existing 320 substitution/indel pairs and corrected 700-pair
containment panel at k=24 with the complete configuration sweep. These panels
test sequence structure and input imbalance; their mutation truths, donor split
and earlier provenance remain unchanged. Weighted ProbMinHash is evaluated
separately on the existing synthetic abundances and official CAMI profiles at
m=32/128/512/2048/8192 and all eight hash seeds. Its truth is probability Jaccard,
not min/max weighted Jaccard; a single eligible weighted method cannot be named
the best weighted algorithm. This is not a reads/error/coverage benchmark.

## Calibration, uncertainty and clear decisions

Synthetic replicates 0–63 are calibration; 64–127 are validation, balanced by
hash seed. Candidate choice and parameter/budget eligibility use calibration
only. Validation reports budget violations, unavailable/empty fractions, MAE,
RMSE, signed bias, P95, zero estimates and relative errors for positive truth.
Never turn an unavailable estimate into a perfect zero or drop a failed row.
The KSSD configuration is fixed per regime/budget by calibration bytes. Other
configurations are fixed in advance. A held-out failure remains a failure.

For each condition compare all eligible candidates, not only the prior two
favorites. Report the calibration nominee, validation candidate set, and a
separate supported-winner flag. Paired bootstrap over independent inputs, with
shared draws across algorithms/budgets within a regime. Control all candidate
contrasts within each regime with a maximum standardized bootstrap deviation;
these intervals are simultaneous within a regime, NOT across all 840 regimes.
Require both an interval excluding zero and a practical improvement of at least
10% and 0.0001 absolute MAE to claim a unique winner. Grey means unresolved,
not equal. Also provide interval widths and descriptive MAE heatmaps, so the
figure conveys more information than a grid of significance labels.

Predeclared pooled summaries stratify by k, payload, size ratio, and containment
zone: zero; low (0<C<=0.1); middle (0.1<C<0.9); high (0.9<=C<1); full (C=1).
Give equal weight to registered cells, retain individual cells, and show seed,
cardinality and k sensitivity. Pooled averages are distribution-specific and
cannot establish a universal switching boundary.

Real accuracy uses genus-cluster uncertainty (reference clusters for controlled
mutations/containment), with seed errors averaged inside biological pairs.
Fixed-size configurations are exact matched bytes. Variable configurations are
chosen on calibration genera/references, fixed for validation; additionally
publish per-configuration tables so poor fixed-grid utilization is visible.
True J and C may explain regimes but are not free online features.

## Downstream checks

Deterministically hash-split the 60 genera into 30 calibration and 30 validation
genera. Within each genus's 15 genomes, evaluate leave-one-out nearest-neighbor
species assignment, exact-J nearest-neighbor recovery, and single-linkage
clustering against species labels. Choose a global J threshold on calibration
genera only from 0.01:0.01:0.99 to maximize genus-mean adjusted Rand index;
freeze it for validation. Include exact-J and fixed-algorithm baselines and
average across tested seeds. Handle nearest-neighbor ties by fractional credit,
not object-order tie-breaking. Cluster resampling unit is genus. Labels are
species taxids, not an alignment-derived definition of species.

Compare descriptive best-in-hindsight errors only as a labeled oracle bound;
do not deploy an oracle selector. Test any task-level calibration choice on
validation genera without reselecting. Better Jaccard MAE need not improve
species retrieval or clustering; record null/reversed downstream results.

## Integrity and delivery

One native process at a time under the existing global benchmark lock. New
workspace `runs/comprehensive_accuracy_20261030`; independent new build and
code/config/source/binary fingerprints. Preserve frozen experiments. Validate
completion, exact input fingerprints, real exact truth and matched cohorts.
Publish compact CSV/gzip observations and summaries, coverage matrices,
candidate maps, error surfaces, downstream results, and a Chinese report with
explicit recommendations, unresolved conditions and scope limits.
