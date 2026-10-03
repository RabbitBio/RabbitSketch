# RabbitSketch2 Stage I protocol

This protocol registers a screening experiment, not a predetermined list of
winners. MinHash means the current FastKMV implementation (`Sketch::FastKMV` in
the pinned upstream source). Legacy MinHash is **not** an additional candidate.
No implementation is deleted by this experiment.

## Scope and estimands

| Family | Estimand | Candidates / follow-up |
|---|---|---|
| A: general sets | Exact canonical-k-mer set Jaccard | MinHash, FracMinHash, HLL, SetSketch, unweighted ProbMinHash, BinDash |
| B: unequal sizes | Directed containment, separately from Jaccard | Native supported estimators; a later dedicated experiment |
| C: abundance | Probability Jaccard | Weighted ProbMinHash; never enters A's frontier |
| D: cardinality | Relative error of distinct-element count | HLL / SetSketch; native estimates are recorded in this pilot |

The implemented runner and pruning decision concern **A only**, including its
unequal-size set regimes. It does not yet benchmark a directed containment
query, abundance, ANI, FASTA throughput, or all-vs-all search.
Family D measurements are descriptive columns, not another pruning decision.

Scope amendment (2026-09-30, user decision): OrderMinHash and the dedicated
sequence-order/rearrangement experiments are excluded from this benchmark
suite and its selection guide. This is a scope choice, not an experimental
pruning result. The completed historical seven-candidate pilot already excludes OrderMinHash.

Scope amendment (2026-10-01, user decision): KSSD is also excluded from subsequent
experiments and the recommendation space. Historical raw data and frozen
protocols remain intact. The six-candidate view must recompute all frontiers;
deleting KSSD rows from an old frontier is incorrect. This is **not** a finding
that KSSD is dominated. `configs/pilot_six.json` specifies the new view.

The 2026-10-01 follow-up additionally validates native Jaccard on 200 local
genome pairs at k=20/24 and three hash seeds. It uses a separately corrected
SetSketch estimator. See [VALIDATION.md](VALIDATION.md) for the amendment,
measurement populations and limits on algorithm-level pruning.

Probability Jaccard, min/max weighted Jaccard, and ordinary set Jaccard must have
separate truth formulas. Winner-label collision probabilities for ProbMinHash
must be checked against the actual algorithm, not inferred from an API method
name or an old comment. See the [ProbMinHash authors' repository](https://github.com/oertl/probminhash).

## Exact synthetic sets

Each pair is constructed from a collision-checked pool of random, canonical
2-bit DNA words. The pool contains exactly `n_a + n_b - intersection` distinct
elements; shared elements are inserted into both sets, and the remaining
elements into just one. Each set is independently permuted. Elements are encoded
as `k-mer + N` carrier strings and passed to each native sequence-update API.
There are no overlapping genomic windows or mutation-model approximations.
This retains KSSD's native input semantics without replacing it with a generic
hash-threshold sampler. It is **not** a shared prehashed-integer microbenchmark.

For larger-set cardinality `n_b` and ratio `r`, `n_a = n_b/r`. The runner requires
integer sizes. Feasibility is checked **before** rounding:

```
J <= min(n_a,n_b)/max(n_a,n_b) = 1/r
intersection = nearest_integer(J * (n_a + n_b) / (1 + J))
J_truth = intersection / (n_a + n_b - intersection)
```

Infeasible cells are written to `plan.json`, not clamped to another regime.
Errors use realized `J_truth`, not the target J. No ANI truth is assigned to
these sets. The generator rejects near-saturated k-mer universes and duplicate
realized regimes. Native verification independently enumerates both carrier
strings, including reverse complements and boundaries.

Every backend and native parameter sees the same pair in a replicate. Replicate
data seeds are independent and derived from the master seed and regime; the
sketch hash seed is fixed. Replicates measure variation over random sets under
that fixed hash implementation, not independently reseeded sketches of a
single genome. Input digests are checked across configurations during analysis.

## Registered grids

`configs/pilot.json` contains k=20/24, larger cardinalities 10^5/10^6/10^7,
ratios 1/10/100, the eight proposed J values, and 20 replicates. After removing
infeasible combinations, it has **84 regimes, 54 configurations, 90,720 jobs**.
The registered grid visits 537,541,920,000 elements
including repeated construction across configurations. It regenerates inputs
per pair, so this is a substantial run, not a quick unit test.

The KSSD primary grid fixes `half_subk=6` and uses `drlevel=0..3`.
The inspected API requires `3 <= half_subk <= 7`, `half_subk <= k/2`, and
`half_subk - drlevel >= 3`. Thus levels 4/5 cannot simply be added at fixed 6;
level 4 needs a different, much larger dictionary, and level 5 is invalid for
this API. Such a secondary experiment would have to be registered separately.
The smoke grid uses half_subk=4 and is expressly ineligible for pruning.

BinDash fixes bbits=16 and varies bins, converting to native 64-bin groups.
SetSketch fixes base=2, a=5 and disables witnesses. ProbMinHash uses unweighted
input, unlimited updates (`max_L=0`), and its native register-value comparison.
The other grids are the native resolutions specified in the JSON.

## Time and memory

Execution amendment (2026-09-30, user requested existing main and all threads):
`examples/main.cpp` accepts a synthetic sample manifest through a small adapter.
One main process handles one algorithm/native configuration/k batch; OpenMP
parallelizes the list using `-t 64` on this host. There are 108 batches of 840
pairs. Only one batch process runs at a time. Data, native parameters, seeds,
truth and pruning thresholds are unchanged. Batch order follows first appearance
in the original randomized plan; jobs keep their randomized relative order.
The earlier isolated measurements are archived separately, never pooled with
this profile. Within-batch timings include concurrent input-generation load.

KSSD's dictionary construction uses global RNG state, so it is initialized once
before the OpenMP loop, following the existing main's shared-parameter policy.
Output lines are synchronized and each row carries a stable job ID. The runner
records the source/CLI/adapter and executable digests, actual threads, compiler/
runtime SIMD paths, CPU model/affinity, configuration, and seeds. Resume rejects
changes to the plan, binaries, sources, runner, execution profile, or host.

| Field | Definition |
|---|---|
| `payload_bytes_a/b` | Measured finalized query-array element count × native element width |
| `sketch_payload_bytes` | A+B payload, the memory objective for one comparison |
| `build_a/b_seconds` | Constructor, native ingestion, and explicit finalization |
| `setup_seconds` | Per-pair access to the already initialized KSSD parameters |
| `build_seconds` | Per-pair setup access + both builds, elapsed under concurrent batch load |
| `query_cold_seconds` | First native Jaccard call, including any lazy estimator cache |
| `query_seconds` | Warm repeated native Jaccard calls; batched until the registered minimum duration |
| `process_peak_rss_before_build_bytes` | Shared process high-water RSS observed after this pair's input generation |
| `process_peak_rss_after_build_bytes` | Shared process high-water RSS observed after this pair's build/finalization |
| `process_peak_rss_bytes` | Shared process high-water RSS observed after this pair's queries |
| batch `shared_setup_seconds` | One KSSD dictionary initialization before the parallel loop |
| batch `wall_seconds` | Entire process duration, including generation, building and repeated queries |

Payload excludes object headers, seed/parameter metadata, spare capacity,
construction scratch, caches, and shared tables. ProbMinHash's unweighted query
payload is its `m` doubles; its retained winner arrays/tree are resident memory,
not an imaginary compressed serialization. BinDash payload is the packed array;
KSSD uses the actual 32/64-bit vector. No hypothetical bit packing is credited.

Linux `ru_maxrss` includes all threads, generation, carrier strings and allocator
retention across pairs. It is not a per-pair/backend-only allocation counter,
and subtracting two high-water marks does not fix that. RSS is reported
separately and is not a payload objective. The analyzer rejects combining
different execution profiles or thread counts in a frontier.
The N-separated workload includes native parsing/hashing; build GB/s on real
contiguous genomes needs a separate experiment. Shared KSSD setup is included
in total batch time and kept separate from the individual pair timing.

The query loop consumes return values and includes the same small function-call
and batch-loop overhead for all backends. Each thread compares its resident
pair with warm caches while other pairs are being processed. Neither these
timings nor total batch throughput measure streaming query-only/all-vs-all work.

## Errors, uncertainty and pruning

All planned jobs remain visible, including timeouts, errors, invalid estimates,
missing rows and empty samples. Invalid/nonfinite estimates are not coerced to
zero. Native estimates for empty sampled sets are retained, but those
configurations cannot enter a frontier. Conditional MAE/RMSE for partially
valid configurations are labeled with valid/failed/missing counts.

Per configuration/regime report MAE, RMSE, bias, p95 absolute error, zero-estimate
fraction, cardinality relative errors where supported, mean measured payload,
median build/query time, and a paired-replicate bootstrap 95% MAE interval.

The primary strict Pareto vector is `(MAE, bytes, build_seconds, query_seconds)`:
no worse in every coordinate and strictly better in one. A separate
accuracy-space frontier uses `(MAE, bytes)`. An algorithm is on the frontier in
a regime if **any** of its eligible configurations is on it. Regimes receive
equal weight; no pooling over the number of genome pairs in a dense region.

A **general-set pruning review candidate** must satisfy all of:

1. The entire registered grid is present without failed workers, the algorithm
   has an eligible configuration in every regime, there are at least 20
   independent replicates, and this is not a smoke run.
2. It is off the strict frontier in at least 90% of feasible regimes.
3. In at least 90% of regimes, every eligible configuration has a strict
   dominator from another algorithm. The paired bootstrap upper bound of
   `MAE(dominator)-MAE(candidate)` must be <= 0.0001, and a competitor must save
   at least 10% in one resource coordinate. Resource comparisons here use
   replicate medians; they are **not** significance tests for runtime.
4. It never has a >=10% resource advantage over all competing feasible
   configurations under the preregistered MAE upper-bound constraints
   0.001/0.005/0.01/0.05 (nor is it uniquely feasible).
5. The decision is restricted to ordinary set Jaccard. Abundance-aware
   metrics and cardinality roles cannot be removed by this decision.

The bootstrap is exploratory, conditional on this fixed grid; it is not a
family-wise multiple-testing guarantee. Independent genome validation, timing
replication, and held-out parameter selection are required before changing the
final advisor. No library removal or definitive algorithm elimination is
performed by these scripts. In particular, a 3-replicate smoke run can never
produce an elimination claim.

## Query reuse

The report computes `T(Q)=build_seconds + Q*query_seconds` at registered Q values
and positive pairwise crossovers. Equal query times give no finite crossing;
negative crossings are irrelevant for Q>=0. Estimates describe the stated cold
setup access / warm resident pair workload under batch load. Shared dictionary
construction is excluded from this pair model. For a collection, shared setup, sketch
build multiplicity, query distribution and caching must be modeled separately.

## Public datasets and subsequent experiments

[Mash supporting data](https://mash.readthedocs.io/en/latest/data.html) supplies
genome identifiers and alignment-derived ANI observations, with sketch-based
Jaccard tables. It does not supply FASTA or exact Jaccard. The downloaded archive
must be audited before joining on pair IDs. The included preparation script
preserves conflicting duplicates with original row numbers and refuses to
silently pick one truth value. Exact canonical k-mer sets must be recomputed
from the original, checksummed assembly versions at each k.

[Dashing's experimental repository](https://github.com/langmead-lab/dashing-experiments)
provides pair selection, pair manifests and exact hash-set comparison code.
The preparation script pins a commit and downloads the small reference files
and their license outside this repository. Its FASTA paths still need assembly
resolution; these files are methodology references, not completed RabbitSketch
measurements. Do not claim to have reproduced Dashing from downloading them.

Follow-up suite: genome accuracy with exact k-mer truth and separately measured
alignment ANI; independent containment with known inclusion; abundance with
the correct probability/weighted truth; cardinality; and
construction/query/scaling on retained candidates. Large eukaryotic inputs
and AVX2/AVX-512 machines are planned coverage, not measured results here.
Mutation probability, realized substitutions and alignment ANI must be stored
separately, especially when indels or rearrangements are introduced. Candidate
datasets such as Podar, CAMI II and Dashing2 large cohorts need their own pinned
manifests and truth audit before use.
