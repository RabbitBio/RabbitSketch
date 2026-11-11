# Threshold-aware pairwise comparison

FastKMV (the current Bottom-k MinHash implementation) and KSSD use early
rejection automatically in the CLI's **list-mode direct all-pairs** path.
The existing `-d` cutoff is passed into the native comparison; its default is
`0.05`. Output still uses the strict condition `distance < cutoff`.

```bash
examples/rabbitsketch --fastkmv -l -i genomes.txt -o fastkmv.tsv -t 64
# KSSD normally uses its inverted index. Force direct comparison to use this path:
RABBIT_FORCE_ALLPAIRS=1 examples/rabbitsketch --kssd -l -i genomes.txt -o kssd.tsv -t 64
```

This does not change the choice of indexed versus direct search. Pair mode
(two files without `-l`) continues to report the full distance. Legacy MinHash
is a separate implementation and is not changed by this feature.

## Library API

```cpp
double d1 = fastkmv_a.distance(fastkmv_b, 0.05);
double d2 = kssd_a.distance(&kssd_b, 0.05);
// Always filter the returned distance: early rejection is an optimization,
// not a replacement for the final cutoff check.
if (d1 < 0.05) { /* retain pair */ }
```

Both Python classes accept `distance(other, max_distance=...)` and
`jaccard(other, min_jaccard=...)`. Ordinary calls without a threshold compute
the full native estimate, preserving existing C++ and Python behavior.
Jaccard thresholds must be finite and in `[0,1]`; distance thresholds must be
nonnegative, with positive infinity disabling pruning.

Rejected pairs may return Jaccard `0` or distance `+infinity`. Pairs meeting
the requested threshold retain the same complete estimate as an unfiltered
query. A returned zero from a thresholded call is not evidence of zero overlap.
KSSD's historical unfiltered zero-overlap distance is `1`; cutoffs at least
`1` use the full estimate to preserve that behavior.

## Rejection bounds

For Mash distance cutoff `D` and k-mer length `kappa`, the necessary similarity
is `Jmin = exp(-kappa*D) / (2-exp(-kappa*D))`. The inverse is conservatively
rounded so pairs at a floating-point boundary survive; the caller applies
the original distance cutoff after evaluation.

FastKMV uses its union Bottom-k estimator, not Jaccard of the two complete
stored sample sets. A full sketch fixes the estimator's denominator at the
common `K`. When both sketches are partial, `max(size_a,size_b)` is a lower
bound on that denominator and gives a conservative necessary match count.

KSSD first projects differing reduction levels to their common sparse
resolution. For aligned sample sizes `a,b`, its set Jaccard is
`c/(a+b-c)`, so a qualifying pair needs at least
`ceil(Jmin*(a+b)/(1+Jmin))` shared tuples. Integer boundaries are adjusted
using the estimator's actual floating-point division, avoiding an unsafe
round-up at equality.

During sorted intersection, the largest possible final overlap is bounded
by `matches_seen + min(remaining_a,remaining_b)`. If this is below the
necessary count, the pair can be rejected without finishing. SIMD kernels
check after blocks; portable loops and scalar tails check periodically.
The optimization preserves threshold recall relative to the same sketch
estimator; it does not claim agreement with exact sequence similarity or
alignment-derived ANI.

## Correctness checks

`pairwise_threshold_test` compares results with an independent materialized
union/membership oracle. It covers full/partial and differing-K FastKMV
sketches, KSSD 32/64-bit tuples, mixed reductions, empty/identical/disjoint
sets, SIMD block boundaries, asymmetric sizes, exact thresholds and adjacent
floating-point cutoffs. It also verifies strict-cutoff output membership
and unchanged unfiltered calls after a thresholded query.
