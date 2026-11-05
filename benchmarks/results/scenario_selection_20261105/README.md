# Scenario selection: independent genome validation

Start with the [Chinese selection guide](../../../RabbitSketch2_scenario_selection.md)
or the [complete report](report.md). The [candidate map](figures/scenario_candidates_k24_8k.png)
explicitly marks cells where no measured configuration meets the standard target.

The experiment uses 450 new genomes from 30 genera held out from the 60-genus
development panel, 14 validation conditions, three k values and four query-payload
budgets. These genera occurred in older broad panels; see
[the scope correction](REPORTING_NOTES.md). SetSketch results use the isolated
corrected implementation. No runtime ranking is made.

The frozen rule does not establish a universal selector or reliable downstream
gains. Relative error reduction and practical target attainment are reported
separately. The full tables, model, frozen decisions, truth/pilot data, source and
input identities, and [verification results](publication_verification.json) are
included here. [Workflow and reproduction](../../scenario_selection/README.md).
