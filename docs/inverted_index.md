# Inverted-index construction

The inverted-index CLI paths for FastKMV (the current MinHash implementation),
legacy MinHash, ProbMinHash, FracMinHash, BinDash and KSSD, along with the
HLL/SetSketch witness paths, use bucket-based parallel construction by default:

1. Extract each sketch's native keys.
2. Count `(key, sample ID)` pairs per input partition and bucket.
3. Prefix-sum the counts and scatter into disjoint ranges of a flat array.
4. Sort the 1024 buckets independently, group equal keys, and drop singletons.
5. Fill a contiguous CSR posting array and record each retained key's range.

FastKMV, legacy MinHash, ProbMinHash and FracMinHash share `buildCSRIndexFromKeys` in
`src/InvertedIndex.h`. BinDash uses the same builder with an exact positional
key and its bucket functor. The MinHash, FastKMV, ProbMinHash, HLL and SetSketch
standalone tests use this builder for their inverted-index paths as well.
KSSD's standalone test delegates to the same CLI implementation.
`buildCSRIndex` remains available for older callers that already own local maps
and for explicit A/B measurements; the CLI and standalone tests no longer call it.

The integer-key default mixes each key **only to choose a bucket**. Stored and
queried keys remain unchanged. This is needed for FastKMV's 53-bit bottom-k
values and FracMinHash's threshold-limited range: selecting the high bits of
the raw 64-bit value would concentrate work in very few buckets. FracMinHash
therefore no longer needs a separate key permutation and re-sort.

KSSD follows the same count/scatter/sort/CSR workflow in `src/KssdIndex.h`,
partitioning on the reduced code's effective bits. Its ordered key space
allows direct offset arrays when appropriate; larger or sparse spaces use
a hash map. This specialization preserves direct addressing without forcing
all algorithms to use the same key representation or lookup structure.

The final key-to-range map is populated serially; sorting and CSR payload
writes are parallel. Flat construction trades per-key allocations and map
merging for a temporary pair array and bucket sorting. Neither construction
method is universally fastest: key count, frequency distribution, thread
count and allocation costs determine the result.

This unification changes construction, not the choice of index versus direct
comparison mode, native estimators, or cutoff. In the CLI's indexed paths,
FastKMV still verifies candidates with its bottom-k union-prefix estimator;
FracMinHash/KSSD use sampled-set intersections; ProbMinHash uses positional
matches; BinDash uses positional b-bit matches with its native correction.
Singleton removal never changes the original sketch sizes used as denominators.

`inverted_index_test` checks original keys, posting ranges and pair counts
against independent map/set references, including low-range and full-width
keys, empty and identical sketches, and ProbMinHash winner keys. Both ordinary
and environment-limited OpenMP teams are covered. BinDash's composite-key
and KSSD's direct-offset checks remain in their existing tests.

Run `python tests/test_index_build_cli.py examples/rabbitsketch` to compare
FastKMV, ProbMinHash and FracMinHash indexed output with their direct CLI
paths, including partial/empty sketches, gzip, nondefault parameters,
threads exceeding input count and limited OpenMP teams. An optional second
binary argument compares with the previous index implementation as well.
