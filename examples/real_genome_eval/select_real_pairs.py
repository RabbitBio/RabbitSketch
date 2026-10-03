#!/usr/bin/env python3
"""
Select ~200 bacterial-genome pairs spanning the full ANI range
for the real-genome validation experiment.

Sampling strategy (deterministic via fixed seed):
  Tier A  (same species)              60 pairs  → expected ANI 95-100%
  Tier B  (same genus, diff species)  60 pairs  → expected ANI 80-95%
  Tier C  (different genus)           60 pairs  → expected ANI 70-85%
  Tier D  (random cross-genus pool)   20 pairs  → expected ANI <75% / no align

Output:  real_genome_eval/pairs.tsv
  columns: pair_id, tier, accA, accB, fileA, fileB, taxA, taxB, organism_A, organism_B
"""

import csv, os, random, sys
from pathlib import Path

GENOME_DIR  = "/home/bigssd/zt/ref/233"
SPECIES_TSV = "/home/bigssd/zt/ref/233_species_groundtruth"
GENUS_TSV   = "/home/bigssd/zt/ref/233_genus_groundtruth"
OUT_TSV     = "/home/user_home/tongz/Finish/RabbitSketch/examples/real_genome_eval/pairs.tsv"

N_TIER_A = 60
N_TIER_B = 60
N_TIER_C = 60
N_TIER_D = 20
SEED     = 20260519

random.seed(SEED)

# ── Load annotations ──────────────────────────────────────────────────────────
print("Loading annotations…")
acc2species = {}
acc2name    = {}
with open(SPECIES_TSV) as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3: continue
        acc, sp, name = parts[0], parts[1], parts[2]
        acc2species[acc] = sp
        acc2name[acc]    = name

acc2genus = {}
with open(GENUS_TSV) as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3: continue
        acc, gn = parts[0], parts[1]
        acc2genus[acc] = gn

print(f"  annotations:  species={len(acc2species)}  genus={len(acc2genus)}")

# ── Index actual files on disk (some accessions may not have FNA) ────────────
print("Indexing genome files…")
file_map = {}   # accession → full path
for entry in os.listdir(GENOME_DIR):
    if not entry.endswith(".fna"): continue
    # GCF_000003135.1_ASM313v1_genomic.fna → accession "GCF_000003135.1"
    acc = "_".join(entry.split("_")[:2])
    if acc not in file_map:
        file_map[acc] = os.path.join(GENOME_DIR, entry)
print(f"  found {len(file_map)} genome FNA files")

# Filter annotations to accessions that exist on disk
def have(acc): return acc in file_map

# Build species → [acc...] and genus → [acc...]
sp2accs = {}
for acc, sp in acc2species.items():
    if not have(acc): continue
    sp2accs.setdefault(sp, []).append(acc)

gn2accs = {}
for acc, gn in acc2genus.items():
    if not have(acc): continue
    gn2accs.setdefault(gn, []).append(acc)

valid_accs = list(file_map.keys())
print(f"  species with ≥2 genomes: {sum(1 for v in sp2accs.values() if len(v)>=2)}")
print(f"  genera with ≥2 species:  "
      f"{sum(1 for g, accs in gn2accs.items() if len({acc2species.get(a) for a in accs}) >= 2)}")

# ── Tier A: same-species pairs ────────────────────────────────────────────────
print(f"\nSelecting Tier A ({N_TIER_A} same-species pairs)…")
candidate_species = [sp for sp, accs in sp2accs.items()
                     if 2 <= len(accs) <= 20]
random.shuffle(candidate_species)
tier_a = []
for sp in candidate_species:
    if len(tier_a) >= N_TIER_A: break
    accs = sp2accs[sp]
    random.shuffle(accs)
    a, b = accs[0], accs[1]
    tier_a.append(("A", a, b))
print(f"  picked {len(tier_a)} pairs")

# ── Tier B: same-genus, different-species pairs ──────────────────────────────
print(f"\nSelecting Tier B ({N_TIER_B} same-genus diff-species pairs)…")
genus_diverse = [g for g, accs in gn2accs.items()
                 if len({acc2species.get(a) for a in accs}) >= 2]
random.shuffle(genus_diverse)
tier_b = []
for g in genus_diverse:
    if len(tier_b) >= N_TIER_B: break
    accs = gn2accs[g]
    by_sp = {}
    for a in accs:
        by_sp.setdefault(acc2species[a], []).append(a)
    sp_keys = list(by_sp.keys())
    if len(sp_keys) < 2: continue
    random.shuffle(sp_keys)
    a = random.choice(by_sp[sp_keys[0]])
    b = random.choice(by_sp[sp_keys[1]])
    tier_b.append(("B", a, b))
print(f"  picked {len(tier_b)} pairs")

# ── Tier C: different-genus, same-domain (i.e., both bacteria) ───────────────
print(f"\nSelecting Tier C ({N_TIER_C} different-genus pairs)…")
tier_c = []
genus_list = list(gn2accs.keys())
random.shuffle(genus_list)
attempts = 0
while len(tier_c) < N_TIER_C and attempts < 5000:
    attempts += 1
    g1, g2 = random.sample(genus_list, 2)
    a = random.choice(gn2accs[g1])
    b = random.choice(gn2accs[g2])
    tier_c.append(("C", a, b))
print(f"  picked {len(tier_c)} pairs")

# ── Tier D: random pairs (probably very distant; many will fail FastANI) ──────
print(f"\nSelecting Tier D ({N_TIER_D} random pairs)…")
tier_d = []
for _ in range(N_TIER_D):
    a, b = random.sample(valid_accs, 2)
    tier_d.append(("D", a, b))

# ── Write output ──────────────────────────────────────────────────────────────
all_pairs = tier_a + tier_b + tier_c + tier_d
out_path = Path(OUT_TSV)
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w") as f:
    f.write("pair_id\ttier\taccA\taccB\tfileA\tfileB\ttaxA\ttaxB\torgA\torgB\n")
    for i, (tier, a, b) in enumerate(all_pairs):
        f.write("\t".join([
            str(i), tier, a, b, file_map[a], file_map[b],
            acc2species.get(a, ""), acc2species.get(b, ""),
            acc2name.get(a, ""), acc2name.get(b, "")
        ]) + "\n")

print(f"\nWrote {len(all_pairs)} pairs → {out_path}")
print("Tier summary:")
for t, n in [("A", N_TIER_A), ("B", N_TIER_B), ("C", N_TIER_C), ("D", N_TIER_D)]:
    print(f"  Tier {t}: {n} pairs")
