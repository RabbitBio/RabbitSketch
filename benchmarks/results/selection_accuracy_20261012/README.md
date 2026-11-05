# Accuracy-only selection, 2026-10-12

Read the [selection guide](../../../RabbitSketch2_accuracy_selection_guide.md)
and [registered protocol](../../selection_accuracy/PROTOCOL.md).

| Evidence | Files | Scope |
| --- | --- | --- |
| New native measurements | `new_measurements.csv.gz`, `new_configuration_accuracy.csv` | All 10,800 observations; input fingerprints, exact counts, estimates, payload and empty flags; no timing fields |
| New paired inference | `new_paired_comparisons.csv`, `overlap_selection_map.png` / `.pdf` | 180 matched contrasts; simultaneous bands, pointwise bands and explicit unresolved cells |
| Frozen seven-way screen, reanalyzed | `frozen_seven_configurations.csv`, `frozen_seven_budget_accuracy.csv`, `frozen_seven_coverage.csv` | 92,960 rows summarized; coverage and empty samples retained; minima are descriptive |
| Existing independent genomes | `real_budget_accuracy.csv`, `real_paired_cluster_intervals.csv`, `real_budget_coverage.csv` | Only Jaccard endpoints; related pairs use the frozen genus-cluster intervals, not an IID pair bootstrap |
| Existing held-out containment | `heldout_containment_budget_accuracy.csv` | Separate API eligibility, availability and pointwise minimum MAE; no new significance claim |
| FracMinHash support audit | `fracminhash_query_support.csv`, `fracminhash_query_support_rows.csv.gz`, `containment_sample_support.png` / `.pdf` | Retained query hashes reconstructed from raw payload; reproduces frozen containment errors and availability |
| Other metrics | `cardinality_accuracy.csv`, `weighted_accuracy.csv` | Original task-specific summaries; probability Jaccard is distinct from min/max weighted Jaccard |
| Provenance | `analysis_manifest.json`, `native_manifest.json`, `publication_manifest.json` | Input/source/binary digests, analysis policy/version, and published artifact digests |

All budgets in pair-comparison tables are A+B query payload, not resident or
construction memory. Cardinality uses per-object register payload. SetSketch
means the experimental corrected estimator. Real-genome and task-specific tables
are reused frozen evidence; only the targeted synthetic measurements are new.

Excluding J=0 and J=1 controls, the 156 new contrasts yield four supported
BinDash preferences, four supported SetSketch-fixed preferences, and 148
unresolved cells. These results do not define a universal switching threshold
or a trained selector. The raw run is stored at
`/home/bigssd/zt/rabbitsketch2_benchmark/runs/selection_accuracy_20261012`.
