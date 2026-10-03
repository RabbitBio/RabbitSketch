# Final robustness validation

This is the final, restricted supplement to `../regimes`: independent bacterial
assemblies and mixed substitution/indel perturbations. It does not extend the
general resolution grid, test a second CPU, or include OrderMinHash.

The frozen design is in `PROTOCOL.md`. Absolute data, run and analysis locations
are defined in `common.py`. Data live on `/home/bigssd/zt`; generated reports are
also written to the repository root. The measured executable is the existing
isolated `rabbitsketch_pilot` main adapter with the validated SetSketch patch,
not a new estimator implementation.

Preparation and execution:

```bash
OPENBLAS_NUM_THREADS=1 /home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/robustness/prepare_natural.py
OPENBLAS_NUM_THREADS=1 /home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/robustness/prepare_mixed.py
OPENBLAS_NUM_THREADS=1 /home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/robustness/run_pipeline.py
```

The pipeline waits for preparation completion markers, computes FastANI in both
directions, freezes the natural pair list before any native sketch results,
runs 94 batches, and performs the 140-pair NUCmer/dnadiff audit. Native/reference
processes share the existing global benchmark lock. Each native batch and
FastANI/NUCmer invocation requests 64 threads. The jobs use the current main
entry point `--screening-manifest … -t 64`. KSSD uses one seed, because this
implementation does not change its reduction with the seed argument.

Analysis and report generation after execution:

```bash
OPENBLAS_NUM_THREADS=1 /home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/robustness/audit_provenance.py
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/rabbitsketch_robustness_mpl /home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/robustness/analyze.py
OPENBLAS_NUM_THREADS=1 /home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/robustness/write_report.py
OPENBLAS_NUM_THREADS=1 /home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python benchmarks/robustness/finalize.py
```

If `analyze.py` runs while the alignment audit is pending, it writes native
statistics but explicitly marks `alignment_complete=false`; it is not a final
delivery. Reports and `finalize.py` require the completed alignment artifacts.

Preparation tests check known-edit identities, independent event replay,
nonoverlap, fingerprint normalization, directional reference missingness,
deterministic balanced selection, and alignment coverage filtering:

```bash
OPENBLAS_NUM_THREADS=1 /home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python -m pytest benchmarks/robustness/test_robustness.py -q
```

Existing data are resumable and guarded by completion markers and manifests.
The native runner rejects altered suite identity. Do not change inputs or
configuration lists under the frozen run directory to add experiments. A new
design needs a separate data/run/analysis directory and a new protocol.

Interpretation rules:

- The natural panel excludes old accession roots, BioSamples and identical
  contig-level sequence fingerprints. It is not guaranteed to exclude near
  clones, circular-origin rotations or every previously seen taxonomic group.
- ANI-reference selection is permitted; sketch-error selection is not used.
  The stratified panel is not a random estimate of the RefSeq population.
- Exact canonical k-mer Jaccard is the sampling truth. FastANI, NUCmer,
  known gapped edit identity and known ungapped identity are distinct columns.
- Paired intervals resample genera/references, conditional on tested sketch
  seeds, and are pointwise without multiple-comparison adjustment.
- Budget envelopes choose configurations using measured bytes alone; frozen
  advisor evaluation keeps the old choice and records failures without fallback.
- This run reproduces construction measurements on new input data, but does
  not fit a new all-vs-all crossover or establish cross-CPU performance claims.
