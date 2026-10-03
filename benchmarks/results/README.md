# Published benchmark evidence

These are the frozen summary tables, plots and reports underlying the algorithm
selection guide. Raw FASTA/FASTQ, native binaries, full per-job output archives,
and the largest per-measurement tables remain outside Git. No estimates were
recomputed or configurations refitted for publication.

- [Stage I: systematic characterization](systematic_20261001/report.md)
- [Stage II: metrics, workloads and advisor](regimes_20261001/report.md)
- [Final external validation](robustness_20261001/report.md)

`regimes_20261001/ani_pair_config.csv.gz` is a lossless compression of the
original CSV and supplies the advisor's observed payload constraints. Markdown
links were made repository-relative. CSV, PNG and PDF evidence is unchanged.
The original `completion.json` files describe the complete external archive,
not only the subset shipped here; use `publication_manifest.json` for the
published subset's hashes and explicit omissions.

The recorded source baseline is `8db783ea75f597068b72705272ae207c94194e83`.
SetSketch-fixed requires the supplied experimental patch. This publication
does not assert that the production estimator has been replaced by that patch.

## External archive

The original runs and data were frozen under
`/home/bigssd/zt/rabbitsketch2_benchmark`. Protocols and code snapshots document
the local inventory and exact commands. A clone contains the evidence and
advisor tables, but not that full genome inventory. Links to omitted local
archives point here rather than pretending that those files are downloadable.
