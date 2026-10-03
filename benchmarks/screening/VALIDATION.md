# Six-candidate follow-up, 2026-10-01

The active candidates are MinHash (= current FastKMV), FracMinHash, HLL,
SetSketch, unweighted ProbMinHash and BinDash. KSSD and OrderMinHash leave the
scope at the user's request, not because of a statistical pruning result.
Do not duplicate FastKMV and MinHash in comparisons.

## Original pilot and analysis populations

The complete pilot has 90,720 successful rows. Removing the KSSD measurements
leaves 84,000 rows, 50 native configurations, 84 feasible regimes, and 20 paired
data replicates. `systematic_analysis.py --baseline-only` recalculates Pareto
frontiers and the original pruning rule from those rows. It does not reuse
seven-candidate frontier labels. Original raw files and snapshots are immutable.

Three primary payload budgets, 2/8/32 KiB **for A+B together**, are filled
by the fixed-size candidates' existing native grids (MinHash uses fewer bytes
when a set has fewer than K distinct elements):

| A+B budget | MinHash K / ProbMinHash m | BinDash bins, 16 bits | HLL / SetSketch p |
|---|---:|---:|---:|
| 2 KiB | 128 | 512 | 10 |
| 8 KiB | 512 | 2,048 | 12 |
| 32 KiB | 2,048 | 8,192 | 14 |

For FracMinHash, choose the largest mean payload among configurations whose
**maximum observed payload** fits the budget and that have no empty sample.
No parameter is chosen using its error. The original scaled grid ended at
10,000, leaving some large-cardinality/small-budget cells untested. A declared
extension measures scaled=20,000/50,000/100,000/200,000 at N_large=10^7,
using the same 28 regimes and 20 paired data replicates (2,240 additional rows).
These values also enter real validation. This fills the budget range rather
than attributing a grid limitation to the algorithm. A coarse native grid can
still leave some budgets underused.
Maximum observed bytes do not guarantee a future realization's hard memory cap.

Supplementary precision coverage uses upper 95% bootstrap MAE / true J <= 25%.
The reuse diagnostic also requires MAE <= .005. These are post-hoc diagnostic
criteria, not replacements for the original registered pruning rule. Pointwise
paired bootstrap intervals are not corrected for the many comparisons. Synthetic
data replicates share hash seed 42; they are independent data realizations, not
20 independent hash seeds. Feasible J is limited by 1 / size ratio; comparisons
of imbalance effects must match J, k and larger cardinality.

## SetSketch implementation correction

The pinned RabbitSketch commit is `8db783ea75f597068b72705272ae207c94194e83`.
Its native update sends an item to one prefix-selected register. The original
joint estimator used a dense SetSketch rank-order-only likelihood. The native
cardinality estimator already describes the partitioned state. Mixing these
models causes a sparse-state Jaccard bias. A static inclusion-exclusion path
also used a harmonic union approximation inconsistent with the corrected
marginal cardinalities, and could return J > 1.

The isolated worktree `source/RabbitSketch-validated` applies the reviewable
`patches/setsketch_partitioned_joint.patch`. This is a correction of the
RabbitSketch implementation's statistical model; it is **not** a claim that
the original SetSketch paper's likelihood is wrong. The primary reference
and its independent implementation are available from the
[SetSketch authors](https://github.com/dynatrace-research/set-sketch-paper).

For register loads lambda_A = n_A/m, lambda_B = n_B/m, shared load z = I/m,
and native constants a and b, define

```
C_r = 1 - exp(-a b^(-r)),   r = 0,...,q
C_(q+1) = 0               (censored final rank)
F(r,s) = exp(-lambda_A C_r - lambda_B C_s + z min(C_r,C_s)).
```

Negative-rank CDFs are zero. The joint cell probability is the four-corner
finite difference of F. Rank zero includes empty registers and items below
the first rank threshold. Stable factorizations of the cell probabilities
avoid catastrophic cancellation. Optimize z over [0,min(lambda_A,lambda_B)],
including exact endpoints; convert to J = z/(lambda_A+lambda_B-z).
The two marginal cardinalities are fitted separately and then held fixed.
This is a **conditional likelihood under a Poissonized register model**,
not a simultaneous three-parameter MLE or a finite-sample optimality claim.

The patch also makes static inclusion-exclusion estimate union cardinality
using the same partitioned model. The old harmonic tail bound is invalid for
that model, so threshold functions currently evaluate the full estimate.
Their inclusion-exclusion screen is not a proven bound on the separate joint
estimate. No index speedup or no-false-negative search claim is made.

Validation includes 450 independent long-double four-CDF-corner probability
checks, 48 native API cases, symmetry/bounds/identity, flat/OO agreement and
20 hashed sparse-state repetitions. The latter gives MAE 0.000466766 on the
chosen regression case. These tests alone do not establish accuracy on all
genomes. The correction is therefore rerun over **all original 11,760
SetSketch synthetic jobs**, with exactly matched input digests, cardinalities,
payloads, native resolutions, seeds and replicate IDs.

The new accuracy view replaces only SetSketch estimates and adds the declared
FracMinHash grid extension. Timing comes from three executions on the same host,
so its four-objective frontier is labeled
descriptive. The correction is not silently merged into the original run or
used to claim a newly preregistered algorithm-level elimination.

## Real-genome validation

`run_real.py` freezes the local 200-pair manifest from
`examples/real_genome_eval/pairs.tsv`. Its former `/ref/233` paths no longer exist;
`prepare_real.py` resolves all 396 matching filenames/accessions from local
`/ref/bacteria/*.fna.gz`, decompresses into the experiment's own data directory,
and records both compressed and uncompressed SHA-256 values. Historical
uncompressed checksums were unavailable; historical measurements are not reused.
The runner records the source/binary identity and hashes every referenced FASTA.
It computes new exact canonical 2-bit k-mer sets
at k=20 and 24, separating contigs with N and treating non-ACGT symbols as
invalid. It does not use the historical k=32 exact-J script (which concatenates
contigs), historical sketch estimates, or ANI as Jaccard truth.

All six backends use the same parsed sequences. Three coordinated hash seeds
(42,101,202), the three fixed budgets above, and twelve FracMinHash scaled
values produce **32,400** native measurements plus 400 exact-truth jobs.
The adapter is reached through the existing `examples/main.cpp` entry and
`--screening-manifest ... -t 64`. A filesystem lock serializes the synthetic
correction, grid extension and real validation: one benchmark process at a time.

Each real pair builds its two sketches privately. Construction includes native
hashing/update/finalization but excludes FASTA parsing and exact truth. Payload
is query-array bytes, not a serialized file or peak process memory. First query
and repeated warm queries are measured separately. Per-pair times are elapsed
times under concurrent 64-thread sample-list load; this is not an isolated
one-core latency benchmark or cross-pair cached construction benchmark.

Real genome pairs reuse assemblies. Results are descriptive pair-level
comparisons after averaging the three seeds; do not treat all 32,400 rows as
independent biological replicates or quote an IID pair confidence interval.
Real and synthetic times are never pooled. The real conditional reuse table
models B+cold+(Q-1)*warm for the same pair; it is not measured all-vs-all or
indexed search performance.

ANI, directed containment, abundance, eukaryotic scale, multihost architecture
and large all-vs-all workloads remain distinct validation targets. In particular,
[ProbMinHash targets probability Jaccard](https://github.com/oertl/probminhash),
which must not be substituted for min/max weighted Jaccard. The unweighted
results here do not establish weighted accuracy or an abundance use boundary.

## Reproduction and review

The original and follow-up runs live under
`/home/bigssd/zt/rabbitsketch2_benchmark/runs/`. Each executed batch has a
manifest, command, row IDs, source/binary identity and raw JSON rows. Frozen
`code_snapshot/` directories preserve the code used by each run.
The archived `real_validation_20261001_stale_paths` directory contains only
the canceled preflight that referenced obsolete paths; it has no measurements
and must not be pooled with the resolved run.

Apply the SetSketch patch only in a separate checkout of the pinned revision.
Build with the existing `main_screening.patch` and screening CMake project.
For a fresh setup (the current environment already has these paths):

```bash
git -C /home/bigssd/zt/rabbitsketch2_benchmark/source/RabbitSketch worktree add \
  --detach /home/bigssd/zt/rabbitsketch2_benchmark/source/RabbitSketch-validated \
  8db783ea75f597068b72705272ae207c94194e83
git -C /home/bigssd/zt/rabbitsketch2_benchmark/source/RabbitSketch-validated apply \
  /home/user_home/tongz/Finish/RabbitSketch/benchmarks/screening/patches/main_screening.patch
git -C /home/bigssd/zt/rabbitsketch2_benchmark/source/RabbitSketch-validated apply \
  /home/user_home/tongz/Finish/RabbitSketch/benchmarks/screening/patches/setsketch_partitioned_joint.patch
cmake -S benchmarks/screening \
  -B /home/bigssd/zt/rabbitsketch2_benchmark/build_validated \
  -DRABBITSKETCH_SOURCE=/home/bigssd/zt/rabbitsketch2_benchmark/source/RabbitSketch-validated \
  -DRABBITSKETCH_SCREENING_MAIN=ON -DCMAKE_BUILD_TYPE=Release
cmake --build /home/bigssd/zt/rabbitsketch2_benchmark/build_validated \
  --target rabbitsketch_pilot synthetic_pilot test_partitioned_setsketch --parallel 16
```

When the partitioned-model header is present, CMake exposes
`test_partitioned_setsketch`; `ctest --test-dir <build> --output-on-failure`
runs the independent model/API regression checks. The Python interface and
input/truth checks use `SCREENING_BINARY` and `SCREENING_MAIN_BINARY` to select
the built executables.

After all three follow-up runs finish:

```bash
OPENBLAS_NUM_THREADS=1 \
  /home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python \
  benchmarks/screening/systematic_analysis.py
```

`results/systematic_20261001/` holds the six-candidate evidence tables, separate
real/synthetic plots and manifests. Every plot is exported as PNG and standalone
PDF. Neither a color in an observed-winner map nor an exact-looking crossover
number establishes a universal decision boundary outside the tested regimes.
