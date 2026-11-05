# Accuracy-only algorithm selection evidence

This workflow asks which algorithm/configuration to consider for a metric,
set-size ratio, overlap structure and query payload budget. It does not use
runtime, SIMD comparisons, construction memory or query reuse to recommend an
algorithm. Read the [Chinese selection guide](../../RabbitSketch2_accuracy_selection_guide.md)
and the preregistered [protocol](PROTOCOL.md).

The new targeted experiment has 10,800 native measurements: 60 feasible regimes,
30 independent generated inputs, three pair budgets, and two algorithms. It
uses k=24, fixed hash seed 42 and an independent input seed. BinDash uses 16-bit
signatures. SetSketch requires the isolated partitioned-joint correction;
the production source is not modified. The binary matches the frozen external
validation binary, and the native source matches the corrected screening source.

The analysis also reuses 92,960 frozen synthetic rows from seven candidates,
the independent-genome results and task-specific validation. These are separate
cohorts, not newly collected genomes. Budget selection uses maximum observed
payload across replicates, then takes the largest fitting tested payload without
looking at error. An empty sample cannot become a winner via a zero estimate.

From the repository root, using the existing external experiment workspace:

```bash
python benchmarks/selection_accuracy/run.py --plan-only
python benchmarks/selection_accuracy/run.py
OPENBLAS_NUM_THREADS=1 python benchmarks/selection_accuracy/analyze.py
python -m pytest -q benchmarks/selection_accuracy/test_accuracy.py
```

Dependencies are pinned in [screening/requirements.txt](../screening/requirements.txt).
The default workspace is `/home/bigssd/zt/rabbitsketch2_benchmark`. `run.py` accepts
`--binary`, `--source`, `--output`, and `--threads`; `analyze.py` accepts `--root`,
`--run`, and `--output`. The runner acquires the existing workspace-wide benchmark
lock, snapshots source/config/protocol before launch, checks resume identity,
and uses the existing native adapter. Config fields needed only by the legacy
performance-analysis schema are placeholders, not decision inputs.

Published [results](../results/selection_accuracy_20261012/README.md) contain every
new accuracy observation, exact counts and input digests, per-configuration
errors, paired comparisons, task summaries, and PNG/PDF figures. Complete raw
JSONL, plans, batch logs and snapshots remain in the external workspace. Full
reanalysis requires that workspace; this directory is not a dataset downloader.

Simultaneous bands use a centered paired bootstrap maximum standardized
deviation across the complete family, sharing resampling indices between budgets
within a regime. They are approximate, conditional on 30 independent inputs and
one hash seed. Supported preferences also require >=10% observed MAE improvement.
Unresolved cells are not evidence of equivalence. The map is explanatory: exact
Jaccard is unavailable to a deployable selector without additional measurement.
