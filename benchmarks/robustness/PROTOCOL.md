# RabbitSketch2 final robustness validation — frozen scope, 2026-10-01

User scope: add independent multi-species genomes and mixed substitution/indel
validation. Do not seek another CPU or extend ordinary screening, HLL, weighted,
containment grids, or OrderMinHash. One measured process at a time, 64 threads,
existing main and isolated validated SetSketch binary. No ML selector training.

## Frozen configurations

Primary k=24, using only previously measured parameters:
MinHash K=128/512/2048; BinDash bins=128/512/2048/8192 (16 bits);
SetSketch-fixed p=12; FracMinHash scaled=5000/10000;
KSSD dr=3,half_subk=6 and dr=4,half_subk=7.
Secondary k=20: MinHash512, BinDash2048, SetSketch12, Frac10000, KSSD3.
Seeds 42/101/202, except KSSD42 because its seed does not alter reduction.
Exact canonical sets at both k. Each panel has 47 jobs, including two exact jobs.
Budgets are pair-total query payload 2/8/32 KiB; variable-size configs are scored
at actual bytes, with feasibility shown. No new parameter is picked from test error.

## Independent natural panel

Use the existing RefSeq assembly snapshot and taxonomy files under
`/home/bigssd/zt/ref`; record source SHA256 and snapshot dates. Target 60 genera,
three species per genus, five distinct BioSamples per species (180 species,
900 genomes). Select by deterministic hashes and balance phyla/families, before
examining sketch results. Eligible metadata: bacteria, full representation,
not excluded from RefSeq, 1–10 Mb, <=300 contigs, known BioSample and lineage.
Sequence QC: >=99% ACGT, N50>=10 kb. Missing lineage/QC failures are logged.

Exclude assembly roots (all versions and known GCA/GCF pairs), BioSamples, and
identical sequence fingerprints from ALL previous benchmark inputs, including
the 10k performance panel. Fingerprints ignore contig order and strand, retain
contig boundaries, and normalize non-ACGT to N. Exclude duplicate selected
BioSamples/fingerprints. Closely related distinct isolates are not claimed to
be independent biological draws; taxonomy/genome clusters are retained.

Candidate pairs: all within each selected genus, plus 120 different-genus
same-family and 120 different-family pairs, chosen deterministically. First
compute FastANI in both directions, minFraction=.2, on the candidate genomes.
Freeze final pairs without any sketch estimates: four hash-selected within-
species pairs per species (coverage panel), add up to 600 reference-qualified
pairs per ANI interval [90,95),[95,99),[99,100] using genus-balanced selection,
and retain the 240 cross-genus/family controls. Union is at most 2760 pairs.
Report any bin shortfall; do not silently expand the population to obtain wins.
No-output ANI remains missing, not zero. All pairs still have exact J truth.

Primary endpoints: equal-byte BD/MH J error, SS/BD high-J (>0.5) error, and
MH/BD build time on identical genomes. Report pair-weighted and genus-balanced
summaries, MAE/P95, actual payload and zero estimates. Bootstrap at genus level
(dyadic weights for cross-genus pairs), conditional on tested sketch seeds.
Report named/unnamed species composition and overlap with earlier taxonomy.
Use NUCmer/dnadiff on 20 hash-selected qualified pairs per ANI bin, without
sketch-error selection, as a separate alignment-reference audit.

Apply prior frozen ANI advisor requests/configurations to the corresponding
new ANI bins. Report accuracy/budget failures; do not refit or select alternatives
on these data. Retain hardware cost model from Stage II and label it as such.

## Mixed mutations

Reuse the 20 prior 5 Mb references only for paired mutation-model sensitivity;
this panel is explicitly NOT new biological validation. Four targets .90/.95/
.97/.99 and substitution:indel affected-base compositions 100:0,90:10,75:25,
50:50. Include a newly generated matched substitution control in the same run.
Insertion and deletion bases share the indel budget equally (integer rounding).
Indel event lengths follow truncated geometric(p=.5), cap10 bases. Place
nonoverlapping events only in ACGT runs; substitutions avoid deleted bases.
Store event scripts and exactly count matches M, substitutions S, inserted I,
deleted D. The target is edit-alignment column identity M/(M+S+I+D), NOT assumed
to equal biological ANI. Also report ungapped homologous identity M/(M+S), edit
burden, event counts, achieved composition, and sequence lengths.

Freeze deterministic inputs and edit scripts before sketch scoring. Validate
event replay, coordinates, output lengths, substitutions and identity counts.
For valid reference length L and indel share q, choose total edited bases near
E=(1-target)*L/[1-(1-target)*q/2], then round I,D,S and record achieved identity.
Total 320 pairs / 340 objects. Exact k-mer truth is primary for sketch sampling.
ANI model residuals use both known edit-alignment definitions, with FastANI
reference reported separately. NUCmer audit: first five deterministically chosen
references × four targets × four compositions =80 pairs, selected in advance.

Compare mutation compositions at matched target and reference, and additionally
stratify by observed exact J. Bootstrap reference clusters; genotype seed
replication is not independent biology. A robust direction must be supported
by observed effect and uncertainty; no automatic robustness claim if it reverses.

## Delivery

Preserve Stage I/II and frozen advice. New runs/data/analysis use
`robustness_20261001` under `/home/bigssd/zt/rabbitsketch2_benchmark`.
Produce raw results, data/command/code hashes, coverage audits, plots/PDF,
paired intervals and an MD addendum in the current repository. Update the
existing analysis document with a linked status/addendum; final paper claims
must reflect successes, failures and reference-definition sensitivity.
