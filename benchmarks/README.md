# RabbitSketch2 benchmark and selection workflows

The benchmark suite compares algorithms by their supported metric, measured
payload, accuracy requirement, construction cost and expected query reuse.
MinHash in benchmark results means FastKMV. KSSD is included in the final
validation; OrderMinHash is outside the benchmark scope but remains a library
backend.

## Workflows

| Workflow | Entry points | Purpose |
| --- | --- | --- |
| Synthetic and real-set screening | [screening](screening/README.md) | Exact Jaccard truth, native resolution sweeps, Pareto diagnostics and paired comparisons |
| Metric and workload validation | [regimes](regimes/README.md) | ANI, cardinality, directed containment, probability Jaccard and actual all-vs-all |
| Final external validation | [robustness](robustness/README.md) | Independent bacterial genomes, known substitution/indel edits, frozen advisor evaluation and alignment audits |
| Constraint-based ANI advisor | [advise.py](regimes/advise.py) | Fixed calibration choice under error, payload and workload constraints; holdout failures remain failures |

The main adapter constructs sketches once, materializes query caches, then
shares immutable sketches across pair queries. Measured jobs run sequentially,
one native main process with 64 OpenMP threads. Exact sets, raw outputs,
commands, input fingerprints and source/binary hashes support reproducibility.
The scripts retain the original experiment workspace paths and frozen data
manifests. Large-scale data preparation is specific to that local inventory;
it is not a zero-download public dataset installer.

## Reports and evidence

- [Detailed algorithm-selection analysis](../RabbitSketch2_algorithm_selection_analysis.md)
- [Independent genomes and mixed-mutation validation](../RabbitSketch2_robustness_validation.md)
- [Published evidence index](results/README.md)
- [Publication validation](PUBLICATION_VALIDATION.md)

The external validation contains 900 genomes, 180 species taxids, 60 genera,
1,967 natural pairs, 320 controlled mutation pairs, 94 native/exact batches and
140 NUCmer reference checks. At k=24 and 8 KiB/pair, BinDash's Jaccard MAE is
0.484 times MinHash's; SetSketch-fixed's high-overlap MAE is 0.816 times
BinDash's. These sampling advantages do not imply transferable ANI error
guarantees: only 6 of 18 frozen natural-ANI advisor choices passed the external
constraints. Performance evidence remains limited to the measured CPU.

The SetSketch results require the isolated
[partitioned-joint correction](screening/patches/setsketch_partitioned_joint.patch).
The suite ships this correction and its independent model/API tests as a
reviewable patch against the recorded source revision. Publishing this suite
does not silently change the production estimator or the frozen experiment.

## Try the advisor with published tables

Install the analysis dependencies in [requirements.txt](screening/requirements.txt),
then run from the repository root:

```bash
python benchmarks/regimes/advise.py \
  --family substitution --ani-bin '99–100' \
  --pair-budget 8192 --target-error 0.001 --queries-per-sketch 1
```

The advisor automatically uses `benchmarks/results/regimes_20261001` when
present. `--data-dir` or `RABBITSKETCH_ADVISOR_DATA` selects another frozen table
directory. Plain or gzip-compressed `ani_pair_config.csv` is supported without
changing the statistical decision rule. Budgets are pair-total bytes; error
is an ANI fraction; query reuse is incident comparisons per object.

This remains a research prototype. Read the external-validation failures and
the reference-definition limits before using its output as a recommendation.
