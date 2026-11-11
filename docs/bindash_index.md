# BinDash inverted-index comparisons

`--bindash --index` in list mode builds and queries a CSR inverted index. It
does not fall back to direct all-pairs comparison. Without `--index`, BinDash
retains the existing packed-signature direct path.

```bash
examples/rabbitsketch --bindash --index -l -i genomes.txt -t 64 -d 0.05 -o distances.tsv
```

Each finalized (densified) bin contributes the exact key `(bin position,
b-bit value)`. Posting traversal counts equal values at equal positions;
the count is converted to Jaccard using BinDash's native b-bit random-match
correction, then to Mash distance. This is the same sketch estimator as
`BinDash::distance`, not a count of common full k-mer hashes. Empty signatures
follow the native packed representation as well.

The index uses 32-bit keys when possible (including the default 1024 bins and
16-bit values), 64-bit keys for wider configurations, and complete composite
keys when position plus value exceed 64 bits. All supported b-bit widths
1 through 64 preserve their full values. Composite hash-table collisions are
resolved by full-key equality. The match counters are 32-bit, supporting more
than 65535 bins.

The minimum matching-bin count is found once per search by binary search over
the native distance formula. This avoids inverse-transform rounding at the
cutoff. Final output uses the original strict condition `distance < cutoff`,
where `cutoff` is the value of `-d`.
There is no extra witness sample or estimated-cardinality ratio filter in
this BinDash index path. Inputs must use the same sketch parameters, as they
do when built by one CLI invocation.

`bindash_index_test` compares posting counts against independent per-position
values and compares filtered output against the native packed estimator. It
covers multiple widths, shifted positions, identical and empty signatures,
thread counts (including environment-limited OpenMP teams), 65536 bins, and
adjacent floating-point distance boundaries. The CLI entry point can be checked
with `python tests/test_bindash_index_cli.py examples/rabbitsketch`.
