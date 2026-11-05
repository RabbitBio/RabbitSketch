# Diagnostic addition, 2026-10-31

Added during the full run after inspecting the first completed synthetic
batches. This is a transparent **post-registration diagnostic**, not a change
to the native grid, original endpoints, paired inference, or winner criterion.
The registered comparisons remain published unchanged.

Some FastKMV containment MAE wins in small-overlap / high-imbalance cells arise
when nearly all estimates are zero. Thus a significant absolute-error advantage
need not imply useful relative precision or detection of weak overlap. Report
the MAE of constant-zero and constant-one predictions for Jaccard/containment,
zero-estimate rates, and whether observed MAE improves on those baselines. These
constant predictions are diagnostics, not selectable sketch algorithms. Do not
advertise a low-MAE winner as a capable overlap estimator solely from this test.

Provide descriptive constraint profiles with calibration-only nomination and
untouched validation results, preserving all failures:

* Jaccard: absolute MAE <=0.001/0.005/0.01/0.05; or relative MAE <=10%/20%/50%.
* Directed containment: absolute MAE <=0.01/0.05/0.1; or relative MAE <=20%/50%.
* Cardinality: mean per-object relative MAE <=1%/2%/5%.

Each profile additionally requires P95 <= twice its MAE limit, output
availability >=99%, and budget violations <=1%. Relative profiles are undefined
at zero truth; label them inapplicable, never divide by epsilon. A nominee is
chosen from calibration-feasible candidates by lowest calibration MAE; ties
use algorithm name. Report all validation-feasible candidates, the fixed
nominee's success/failure, and no-candidate outcomes. These are empirical
constraints on the measured distribution, not statistical error guarantees.
The additional thresholds are illustrative operating requirements, not a
retrospectively optimized selector or universally preferred error tolerances.

This diagnostic helps distinguish a supported MAE winner, a configuration that
meets a stated precision requirement, and a downstream-useful algorithm. None
of the three automatically implies the other two.

Also report the smallest calibration-feasible tested budget for each profile,
retaining validation failure instead of switching budgets afterward. A separate
descriptive error/payload frontier uses observed validation MAE and maximum
payload of complete configurations. It is a measured-point visualization, not
a statistical dominance claim or a deployable validation-tuned policy.
