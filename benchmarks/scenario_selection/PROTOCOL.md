# Observable scenario selection: registered 2026-11-05

This extension tests **when a precision-oriented choice works**, with no runtime,
scalar/SIMD or construction-RSS ranking. It does not assume that a selector wins.
The earlier Oct30 suite remains frozen. The native binary and isolated corrected
SetSketch source are exactly those of that suite; SetSketch results therefore
remain conditional on that experimental correction.

## Independence and order

All 60 previously studied genera (900 genomes) are development data in this
extension, including their earlier validation half. They are **not fresh evidence**.
Select 30 previously unused genera, three species each, five genomes per species
from the local NCBI assembly snapshot. Exclude previous accession roots and paired
GenBank/RefSeq aliases, BioSamples, and sequence fingerprints invariant to record
order and reverse complement. Selection uses taxonomy and sequence QC only:
1–10 Mb, at most 300 contigs, at least 99% ACGT, N50 at least 10 kb. Select with
deterministic hash ordering, balancing phyla/families before hash ties. No selection
on sketch errors, species separability or exact Jaccard. If fewer than 30 genera
pass, report the deficit rather than recycle development genera.

Within each species, a deterministic split assigns two reference genomes and three
query source genomes. The independent unit for intervals is the genus, not seeds,
pairs, reads, or multiple derivatives of a source. Development and validation use
disjoint hash seeds for seed-sensitive backends. KSSD uses its existing fixed
parameter construction, one seed, and is not represented as three repetitions.

Freeze this protocol/config and data manifests before native outcomes. Finish all
development measurements, fit thresholds and the advisor, and hash their artifacts
**before any validation algorithm/exact measurements**. Validation pilot results
may then be measured; freeze every decision before running the validation outcome
grid. No validation retuning. Bugs require documented amendments and reruns of
affected artifacts; no silent replacement of a failed selector.

## Workloads and estimands

The ten fixed conditions in config.json comprise full genomes, 1/10/100 kb
fragments, incomplete assemblies with/without 20% same-genus other-species
contamination, and simulated 150-base reads at two coverages and three substitution
settings. These are controlled derivatives of real assemblies, **not real sequenced
read datasets**. Reads sample valid positions within contigs with replacement;
errors are uniform substitutions. Short contigs truncate a requested single
fragment, with actual lengths reported. Incomplete assemblies sample nonoverlapping
10 kb tiles. Record boundaries are retained; no artificial junction k-mers.

Use k=20/24/32 and pair-query-payload budgets 512/2048/8192/32768 bytes. Compare
FastKMV, FracMinHash, HLL, corrected SetSketch, KSSD, unweighted ProbMinHash and
BinDash, using the same native grid as Oct30. Legacy MinHash is an implementation
baseline, and order similarity and weighted probability Jaccard are different
estimands: do not mix them into this seven-backend set-Jaccard ranking. The earlier
weighted results remain the separate evidence for weighted ProbMinHash.

For each genus/condition, score 9 queries against the 6 reference genomes. This is
a **known-genus candidate bank**, not unrestricted database classification. Also
measure all query-query edges for clustering and full-genome within-genus edges.
Evaluate query→reference Jaccard, directed query containment, query cardinality
relative error, exact-Jaccard nearest-neighbor recovery, source-species retrieval,
and query clustering ARI. Exact canonical k-mer sets provide set truth, not aligned
ANI. Uniform fractional tie credit; self and query source genomes are absent from
the reference bank. Each source contributes equally.

Open-set stress test: remove the query's species from its candidate references,
then measure false acceptance. Fit acceptance thresholds on development only, by
observable input kind/k/backend/parameter, at a 90% development known-query recall
target. Report achieved known recall and unknown false acceptance, including zero
score ties; do not claim an open-set detector from a closed-set success.

Clustering uses single linkage and thresholds {0.01,...,0.99}, fitted on development
genera only by input kind/k/backend/parameter, with smallest threshold breaking
ties. Report exact-set task scores as an empirical task ceiling/reference, not as
proof that the k-mer metric matches biological labels.

## Observable advisor and baselines

One compatible backend and parameter is selected for an entire query workload.
The advisor never compares sketches from different algorithms, nor changes an
already-built index's backend per query. Rebuilding a bank and the pilot's time
cost are explicitly outside this precision-only validation.

A seed-919 FastKMV pilot of K=32 consumes at most 256 payload bytes per object,
512 per pair. It is discarded before final sketches are built. Report its total
temporary payload separately; the main budget is final query payload, not RSS or
pilot-plus-final persistent memory. A deployment retaining both must add their
payloads. The pilot is not a free operation or a calibrated confidence interval.

Metadata features: input kind, k, log query/reference bases, query record count and
record length. Pilot features additionally summarize estimated cardinalities,
size ratios, Jaccard, bounded directed containment, zero estimates and top-score
margins. Exact counts/Jaccard, taxonomic identifiers, nominal coverage, true error
rate, contamination proportion, and correctness labels are **forbidden inputs**.

Use a fixed nearest-evidence model, no hyperparameter search: robust median/IQR
scaling fitted on development; nearest eight **distinct genera** with matching
input kind and k; equal weights. Each development genus supplies at most one
neighbor, preventing repeated derivatives from dominating. Fit separate losses
for each objective. A metadata-only variant is an ablation. Fixed baselines are
the best development configuration per objective/k/budget/input kind, plus every
individual backend with its development-selected configuration. Native variable
payload configurations must pass an observable size guard (1.5× pilot cardinality,
or input-bases upper bound for metadata). Actual budget failures remain failures.
Report uncertainty as candidate sets/empirical neighboring risks, not guarantees.

Loss is the mean over the declared workload; unavailable or over-budget pair/
query outputs contribute loss 1. Jaccard/containment use absolute error, cardinality
relative error is clipped to 1 for selection (unclipped errors also reported),
retrieval uses 1-credit, and clustering (1-ARI)/2. Conditional MAE/P95, availability,
zero estimates and budget violations are reported separately. Predictions and
feasibility gates must not silently hide failed outputs.

Accuracy profiles are fixed in config.json. Report recommendation coverage, false
feasible rate, actual target attainment, loss/regret against the same-budget
post-hoc oracle, and comparisons against the strongest development-fixed baseline.
The oracle is an offline diagnostic and never an input. Thresholds, candidate
sets and models are fixed before validation. Out-of-support features are flagged;
unsupported estimands are rejected explicitly. Retain unfavorable results.

## Evidence and publication

Paired genus bootstrap (4,000 draws) compares advisor vs fixed choices. Intervals
are pointwise unless explicitly adjusted; numerous scenario cells cannot be
read as simultaneous winner claims. Aggregate seeds within genera first. Publish
scenario/budget/metric tables, complete frozen decisions, calibration choices,
failure counts, taxonomy/provenance and hashes. Preserve all native archives and
exact counts externally. Figures must show failures and uncertainty as well as
point estimates. Claims are limited to this taxonomy snapshot, known-genus banks,
controlled conditions, three k values and the frozen native implementations.
