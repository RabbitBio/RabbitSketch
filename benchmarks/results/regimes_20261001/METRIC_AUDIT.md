# Metric contracts used in this experiment

Native source: isolated `RabbitSketch-validated`, parent commit
`8db783ea75f597068b72705272ae207c94194e83`, experimental SetSketch patch.
No production API or public enum was renamed as part of this benchmark.

* **ANI:** native `ani()` is `exp(-D)`, while `1-D` is a different Mash-based
  identity approximation. Both are reported. Neither exact-set transformation
  is itself an alignment-derived ANI ground truth. Signed error decomposes
  exactly into sampling plus model/reference components; absolute MAEs do not
  add, and bias cancellation can lower total MAE below the exact-J baseline.
* **BinDash containment:** `BinDash::cardinalityEstimate()` returns infinity
  when every bin is occupied. Its low-level `.containment()` can then return
  zero; `UnifiedQuery::populateCardinalityMetrics` instead warns that derived
  metrics are unavailable. The benchmark treats saturation as unavailable.
  The separately labelled BinDash+HLL composition charges both sketches and
  both construction passes and does not use exact-cardinality side information.
* **SetSketch/HLL containment:** use joint Jaccard and marginal cardinalities,
  then bound the common intersection by both cardinalities, as UnifiedQuery
  does. Low-level native containment is recorded separately. The SetSketch
  results require the experimental partitioned-joint estimator patch.
* **KSSD containment:** the pinned UnifiedQuery supports KSSD Jaccard/ANI and
  their distances, but not directed containment. KSSD containment figures here
  are an explicitly experimental sampled-intersection diagnostic, not a claim
  that the public API already supports that metric. The adapter's generic raw
  `containment_native_*` field also contains this diagnostic for KSSD; it is
  excluded from the published native-containment summaries.
* **Weighted ProbMinHash:** the winner array estimates probability Jaccard,
  `J_P = sum_i 1 / sum_j max(a_j/a_i, b_j/b_i)`, summing over positive common
  coordinates. One aggregated weight is inserted per key. The header/comment
  claims that winner collisions estimate `sum(min)/sum(max)` are incorrect.
  Scale invariance is a direct counterexample: `b=10a` gives `J_P=1`, whereas
  min/max Jaccard is 0.1. This is covered by the native test and full experiment.

The current ProbMinHash native/UnifiedQuery weighted-containment derivation
also assumes min/max Jaccard and cannot generally be derived from `J_P`.
For `a=(1,3)`, `b=(3,1)`, exact `J_P=1/2`; plugging it into that derivation
gives containment `2/3`, whereas the min-mass containment is `1/2`.
Weighted ANI and those derived containment quantities are outside this
benchmark's recommendation. The generic frozen weighted raw records contain
unused `containment_derived_*` numbers; these are not valid weighted estimands
and are intentionally absent from weighted analysis/figures. The scored raw
`jaccard` field in weighted mode is explicitly **probability Jaccard**.

The exact J_P implementation uses sorted abundance ratios and prefix/suffix
sums. It is checked against an independent quadratic formula on 100 random
weight pairs. CAMI validation uses the ten official per-genome abundance
profiles, verified byte-for-byte against the local original setup archive;
it does not claim read-to-k-mer abundance reconstruction or coverage validation.

Primary sources: [ProbMinHash author repository](https://github.com/oertl/probminhash),
[probability Jaccard definition](https://arxiv.org/abs/1809.04052),
[Mash distance definition](https://mash.readthedocs.io/en/latest/distances.html).
