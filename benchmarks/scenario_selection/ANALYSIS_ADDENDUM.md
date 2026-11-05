# Implementation review before validation (2026-11-05)

The registered datasets, native grid, objectives and advisor rule are unchanged.
No validation algorithm/exact outputs existed when these notes were written.

1. An open-set scoring review found that availability was initially checked on
   the original six-reference gallery. Once the query species is removed, the
   availability and budget checks must use only the remaining four references.
   This has been corrected and independently tested. The first four development
   configuration summaries were archived outside the publication and recalculated;
   no model or task threshold was frozen from them. Native measurements are unchanged.
2. Report constant-zero and constant-one MAE and positive-truth relative MAE for
   Jaccard/containment. These are diagnostic columns, not new trained objectives or
   additional candidate algorithms. They expose an already known failure mode of
   absolute error at very small Jaccard: outputting zero can pass an absolute-error
   target without retaining useful ranking information. Task-specific retrieval
   results and zero rates must accompany such cells.
3. The first scoring launch preceded completion of pilot feature extraction and
   stopped with a missing-file error before producing a score. It was restarted
   after feature extraction completed. This was an execution-order issue, not a
   change to data or outcomes.

## Additional held-out workload conditions, before model/validation measurements

The ten original conditions alone test new organisms within the same workload
grid. To test transfer of an unchanged rule, add four **validation-only** conditions
from `validation_extra_conditions.json`: a 5 kb fragment, 30% tile assembly with
10% same-genus other-species contamination, 150 bp reads at 0.5×/2% substitution,
and 250 bp reads at 0.5×/1% substitution. No native or pilot results from these
conditions are available when registering them. They use the same 30 independent
validation genera and do not add biological replicates. Calibration stays at ten
conditions; no hyperparameter or threshold is retuned. Report original-condition
and unseen-condition aggregates separately, alongside all 14 conditions. The
primary 72 objective/k/budget contrasts cover all registered validation conditions.

The unchanged original data-preparation artifacts are archived before extending
the validation manifests. The extension has its own code/config hashes and a
completion marker. Genomes, accession exclusions, native configurations, k values,
budgets and random seeds are unchanged. The long-read-length condition is expected
to challenge the observable support check; it is still a simulated short-read
profile, not a real sequencing dataset or a long-read technology model.

Precision-profile gates use at least 99% valid in-budget outputs and at most 1%
budget violations, in addition to both registered mean and P95 error targets.
Neighboring P95 values are averaged for prediction; they are not a statistical
upper confidence bound. Conditional false-feasible rates must be reported.
