# Reporting clarifications (model/decisions remain frozen)

The advisor minimizes predicted failure-penalized mean loss. Accuracy profiles
screen its recommendation afterwards; this is not a trained MAE/P95-constrained
optimizer. Therefore a selected configuration failing a target cannot establish
that the memory budget is intrinsically insufficient.

Before the final validation analysis, add a diagnostic that checks **every measured
native alternative** against the same registered accuracy/availability/budget gate.
Report separately: (1) the chosen configuration fails although some measured
alternative passes; (2) no measured alternative passes; (3) a prediction of success
is wrong. This uses validation truth only for retrospective scoring and never
changes a recommendation, parameter, feature, threshold, or fitted model. Even
case (2) is limited to the measured configuration grid, not all possible algorithms
or parameter values. Constant-zero diagnostics remain necessary for small-J cells.

The score for an out-of-support candidate can still be measured retrospectively;
the range check is heuristic, not a guarantee against distribution shift. The four
pre-registered unseen conditions are always reported separately irrespective of
whether that heuristic labels their observable features as within range.

## Taxonomic novelty: scope correction

The initial registration/commentary called the 30 validation genera previously
unused. The implementation excluded the 60 genera in the development inventory,
whereas accession/BioSample/sequence exclusions used **all historical inputs**.
A later metadata-only audit of all earlier input manifests found 1,047 mapped
historical genera; all 30 validation genera had appeared in older broad panels.
The 450 selected genomes remain new by accession, BioSample and sequence, and no
validation genus overlaps the current 60-genus development set. The claim is
therefore **new genomes and a genus-disjoint development/validation split**, not
genera never studied anywhere in the project's history. This is an explicit
scope correction to the overbroad registration wording. No inputs, model,
thresholds or decisions are changed in response to this audit.

`independence_audit.json`'s zero `overlap_genera` field refers specifically to the
development inventory. `historical_genus_audit.json`, reproducible with
`audit_history.py`, reports the broader historical overlap. Do not omit either.

## Conservative simultaneous inference

The implementation additionally reports genus bootstrap max-t intervals across
the 72 primary task/k/budget contrasts and all eligible algorithm contrasts within
each scenario cell. A material improvement also requires a 10% relative and
0.0001 absolute loss reduction, consistent with the earlier accuracy study.
This conservative reporting rule does not change any frozen decision or target.
The original registration did not specify this full simultaneous-inference
implementation; distinguish the registered experimental grid from these analysis
details. There is no claim of global multiplicity control across every scenario
cell or across both the primary rule and its metadata ablation together.

## Candidate-map display

The candidate contrast table ranks relative loss among available backends. It
must not recommend a statistically best configuration when every measured
configuration still fails the registered precision target. Display this failure
directly in the map, using the already defined feasible-alternative diagnostic.
Mixed-feasibility cells are marked separately. The underlying contrasts and all
frozen decisions remain unchanged, including the three primary containment
mean-loss improvements and their unchanged practical attainment rates.
