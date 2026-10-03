# RabbitSketch2 metric and workload validation — 2026-10-01

Stage I is frozen. User scope: MinHash is FastKMV; restore KSSD; exclude
OrderMinHash. One native main process, 64 OpenMP threads, one configuration
at a time. Use the isolated validated SetSketch build and identify its patch.

## Endpoints and truth

1. ANI: original Mash archive has two conflicting observations per pair ID.
   Do not choose one arbitrarily. Use accession-version-matched local genomes
   as a **142-genome Mash accession subset**, with newly computed FastANI
   references (FastANI is an approximate, alignment-free ANI estimator, not
   exact nucleotide alignment truth). Also use existing 200 diverse pairs and
   controlled substitutions of real references with known positional identity.
   No-output FastANI pairs remain missing; never assign them zero identity.
   Record both directions and alignment fractions. Primary natural-pair analysis
   requires both directions and both alignment fractions >= 0.2.
   Separate sampling error against transformed exact canonical Jaccard from
   total error against the independent identity reference. Report RabbitSketch
   native exp(-D) and Mash 1-D separately, in fraction units and percentage points.
   Use k=20,24, fixed native grids, three sketch seeds. Report MAE, RMSE, p95,
   bias and zero estimates, ANI bins 80–90,90–95,95–99,99–100%, and threshold
   errors within 95±1 and 99±0.5 percentage points.
2. Cardinality: exact distinct canonical 2-bit words; real FASTA windows at
   10 kb,100 kb,1 Mb,5 Mb,50 Mb. HLL vs validated SetSketch at p=8,10,12,14,16;
   report relative error, single-object query payload, build and card-query costs.
3. Directed containment: exact intersection / query cardinality, both directions.
   Real reference fragments of 5 Mb,2 Mb,1 Mb,500 kb,100 kb,50 kb,10 kb,
   nominal source overlap 1,.8,.5,.2,.1. Exact k-mer truth supersedes nominal
   sequence fractions. Run MinHash, BinDash, FracMinHash, SetSketch, KSSD.
   Record native containment separately from containment derived from native
   joint Jaccard and native cardinalities; no free exact-cardinality oracle.
   Empty samples and saturated/nonfinite cardinalities are failures, not good
   estimates. Report availability as well as error.
4. Workload: build every distinct input once through existing main, prepare
   immutable query state before parallel queries, and evaluate actual distinct
   all-vs-all pairs. Separate parsing, build, query, and output costs. Repeated
   passes are labelled cached passes, never extrapolated as measured larger N.
   Only query payload is algorithm space; peak RSS includes resident FASTA,
   construction state, and shared dictionaries and is labelled accordingly.
   Q=P/N denotes charged comparisons per object: T/N=B+Q*t; incident reuse
   R=2P/N gives T/N=B+(R/2)*t. For all-vs-all Q=(N-1)/2, R=N-1.
   Crossovers require positive finite solutions and compatible accuracy/memory
   constraints. Do not treat all-vs-all throughput as isolated query latency.
5. Weighted ProbMinHash: audit winner comparison versus priority comparison.
   Truth is probability Jaccard J_P, not sum(min)/sum(max). Insert one aggregated
   weight per key. Test uniform/lognormal/Zipf, proportional rescaling, and CAMI
   II Marine official per-genome abundance profiles. CAMI profile results are
   profile-level validation, not automatically read/k-mer-level validation.

## Selection and uncertainty

Compare configurations under measured pair-total payload 2/8/32 KiB (single
objects 1/4/16 KiB when symmetric); variable-size sketches may not fit every
pair. Configuration selection may use cardinality/bytes but not held-out error.
Keep biological replicates distinct from sketch-seed replicates. For controlled
experiments split reference IDs deterministically into calibration/validation
halves. Natural Mash subset also splits by genome, excluding cross-split pairs
from selection validation. Selection maps use task-specific error constraints;
mark infeasible cells rather than name a winner. Empirical boundaries apply to
tested regimes and this implementation/hardware. No new global pruning claim
is implied by a missing configuration or a limited native parameter grid.

## Task-specific supplements (declared before inspecting containment outcomes)

The native BinDash occupancy cardinality saturates at genome-scale cardinalities.
Measure a separate **BinDash+HLL** composition, not a silent correction: bins
384/1536/6144, plus HLL p=8/10/12 respectively, charging both payloads and both
construction passes. Pair payloads are 2/8/32 KiB. No exact-cardinality oracle.
Also include KSSD drlevel=4 with half_subk=7 in ANI, containment, and workload
to cover its native compressed regime, recording the much larger shared
dictionary setup. The main KSSD dr=2,3 uses half_subk=6. KSSD's current seed
field does not change its reduction table or sequence-update sampling, so
do not count identical seed reruns as independent statistical replicates.
These additions use a separate runner, binary snapshot and output directory.

## Containment validation correction

The initial fragment generator paired each reference with a donor from the
opposite split. That run remains an exploratory diagnostic and is excluded
from final containment summaries. `containment_holdout_20261001` replaces the
donor with the next genome **within** the same split, then rebuilds all 700
pairs, exact truths, and all configurations including the two supplements.
The two source-genome sets are now disjoint. No ANI or performance input changes.
Final containment primary comparison uses J-derived, [0,1]-bounded estimates
with native cardinalities, labelled as such. Native outputs and out-of-range
rates remain separately available; neither path uses exact cardinality.

Natural ANI is additionally checked with MUMmer 4.0.1 NUCmer `--maxmatch -c 65
-l 20 -t 64` and dnadiff's 1-to-1 AvgIdentity on a deterministic stratified
subset selected without sketch errors. This is a separate alignment-derived
reference; do not silently pool it with FastANI. Report aligned fractions and
the dnadiff report's 0.0001-percentage-point rounding. Threshold confusion tables
include both all pairs in the window and a separately labelled view excluding
exact-threshold ties (fraction tolerance 1e-6).

The high-compression KSSD comparison adds BinDash 128 bins (256 query bytes
per object) for ANI and actual all-pairs, and all-pairs timing for the existing
MinHash K=128 / BinDash 512-bin configurations. Otherwise a compressed KSSD
would be compared only with much larger signatures. This targeted supplement
is isolated in `compressed_validation_20261001`; it is not another ordinary
Jaccard parameter sweep.
