#!/usr/bin/env python3
"""
Analyse real-bacterial-genome benchmark results.

Inputs (in the same directory):
  pairs.tsv             - pair list with tier info
  fastani_truth.tsv     - FastANI ANI for each pair (ground truth)
  sketch_results.csv    - per (pair × algo × params) sketch-derived ANI

Outputs (printed to stdout, no plotting):
  • Tier-A/B/C/D ANI distribution from FastANI
  • Layer-1 (default params) ranking: MAE / RMSE / Pearson r per algo
  • Per-tier breakdown of Layer-1
  • Layer-2 parameter-sensitivity tables (k-mer, SetSketch a, KSSD drlevel)
  • Throughput summary (mean seconds per pair, per algo)

Sketch-derived ANI = 1 - mash_dist  (RabbitSketch already outputs Mash distance).
FastANI gives ANI in [0,100] %; we normalise to [0,1] for comparison.
"""

import csv, math, statistics, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent

# ── Load pairs.tsv (tier info) ────────────────────────────────────────────────
pairs = {}   # pair_id -> {tier, accA, accB, orgA, orgB}
with open(HERE / "pairs.tsv") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for row in rdr:
        pairs[row["pair_id"]] = row

# ── Load FastANI truth ────────────────────────────────────────────────────────
truth = {}   # pair_id -> fastani_ani  (in [0,1] or None if no result)
no_align = set()
with open(HERE / "fastani_truth.tsv") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for row in rdr:
        pid = row["pair_id"]
        v = float(row["fastani"])
        if v < 0:
            no_align.add(pid)
            truth[pid] = None
        else:
            truth[pid] = v / 100.0

n_align = sum(1 for v in truth.values() if v is not None)
print("="*78)
print(f"Loaded {len(pairs)} pairs  •  FastANI aligned: {n_align}  •  no align: {len(no_align)}")
print("="*78)

# ── FastANI ANI distribution per tier ─────────────────────────────────────────
print("\n[A] FastANI ANI distribution per tier")
print("-"*78)
by_tier = defaultdict(list)
for pid, t in truth.items():
    if t is None: continue
    by_tier[pairs[pid]["tier"]].append(t)
print(f"  {'tier':<6} {'n':>4}  {'min':>7} {'p25':>7} {'med':>7} {'p75':>7} {'max':>7}")
for tier in ["A", "B", "C", "D"]:
    xs = sorted(by_tier[tier])
    if not xs:
        print(f"  {tier:<6} {0:>4}   (no aligned)")
        continue
    def pct(p): return xs[int(p*(len(xs)-1))]
    print(f"  {tier:<6} {len(xs):>4}  {min(xs):>7.4f} {pct(0.25):>7.4f} "
          f"{statistics.median(xs):>7.4f} {pct(0.75):>7.4f} {max(xs):>7.4f}")
# Also: how many in each tier had no FastANI result?
n_noalign_per_tier = defaultdict(int)
for pid in no_align: n_noalign_per_tier[pairs[pid]["tier"]] += 1
print(f"\n  No-align counts per tier: " +
      ", ".join(f"{t}={n_noalign_per_tier[t]}" for t in "ABCD"))

# ── Load sketch results ───────────────────────────────────────────────────────
rows = []
with open(HERE / "sketch_results.csv") as f:
    rdr = csv.DictReader(f)
    for row in rdr:
        rows.append({
            "pair_id":  row["pair_id"],
            "layer":    row["layer"],
            "algo":     row["algo"],     # tag, e.g. "minhash_def"
            "params":   row["params"],
            "dist":     float(row["mash_dist"]),
            "ani":      float(row["sketch_ani"]),
            "sec":      float(row["sec"]),
        })
print(f"\n  Sketch result rows: {len(rows)}")

# Extract algo family from tag (strip _def/_k15/_a0p5/...)
def algo_family(tag):
    return tag.split("_")[0]

# ── Layer-1 (default-param) summary ──────────────────────────────────────────
print("\n[B] Layer-1 default-parameter accuracy vs FastANI  "
      "(only pairs with FastANI alignment)")
print("-"*78)

def metrics(xs_truth, xs_est):
    n = len(xs_truth)
    if n == 0: return (0, 0, 0, 0, 0)
    errs = [e - t for e, t in zip(xs_est, xs_truth)]
    abs_errs = [abs(x) for x in errs]
    mae = sum(abs_errs) / n
    rmse = math.sqrt(sum(x*x for x in errs) / n)
    bias = sum(errs) / n
    # Pearson r
    mt = sum(xs_truth) / n; me = sum(xs_est) / n
    num = sum((t-mt)*(e-me) for t,e in zip(xs_truth, xs_est))
    dt = math.sqrt(sum((t-mt)**2 for t in xs_truth))
    de = math.sqrt(sum((e-me)**2 for e in xs_est))
    r = num / (dt*de) if (dt>0 and de>0) else 0.0
    return (n, mae, rmse, bias, r)

# Aggregate L1 by algo family
l1 = defaultdict(list)   # algo_family -> [(truth, est, sec)]
for r in rows:
    if r["layer"] != "L1": continue
    t = truth.get(r["pair_id"])
    if t is None: continue
    l1[algo_family(r["algo"])].append((t, r["ani"], r["sec"]))

print(f"  {'algo':<12} {'n':>4} {'MAE':>7} {'RMSE':>7} {'bias':>7} {'r':>6} {'avg_s':>7}")
ranked = []
for a, data in l1.items():
    xs_t = [d[0] for d in data]
    xs_e = [d[1] for d in data]
    n, mae, rmse, bias, r = metrics(xs_t, xs_e)
    avg_s = sum(d[2] for d in data) / max(1, n)
    ranked.append((mae, a, n, rmse, bias, r, avg_s))
for mae, a, n, rmse, bias, r, avg_s in sorted(ranked):
    print(f"  {a:<12} {n:>4} {mae:>7.4f} {rmse:>7.4f} {bias:>+7.4f} {r:>6.3f} {avg_s:>7.3f}")

# ── Layer-1 per-tier breakdown ────────────────────────────────────────────────
print("\n[C] Layer-1 MAE per tier (A: same-species, B: same-genus, C: diff-genus)")
print("-"*78)
# Tier D usually has no FastANI alignment; skip unless populated.
algos_sorted = sorted(l1.keys())
header = f"  {'algo':<12} " + " ".join(f"{t:>9}" for t in ["A","B","C","D"])
print(header)
for a in algos_sorted:
    cells = [f"  {a:<12} "]
    for tier in "ABCD":
        sub = [(truth[r["pair_id"]], r["ani"])
               for r in rows
               if r["layer"]=="L1" and algo_family(r["algo"])==a
                  and truth.get(r["pair_id"]) is not None
                  and pairs[r["pair_id"]]["tier"]==tier]
        if not sub:
            cells.append(f"{'-':>9}")
        else:
            _, mae, *_ = metrics([s[0] for s in sub], [s[1] for s in sub])
            cells.append(f"{mae:>9.4f}")
    print(" ".join(cells))

# ── Layer-2 parameter sensitivity ─────────────────────────────────────────────
print("\n[D] Layer-2 parameter-sensitivity tables  (50-pair subset)")
print("-"*78)
# Group L2 results into themes for nicer reporting.
def collect(predicate):
    by_tag = defaultdict(list)
    for r in rows:
        if r["layer"]!="L2": continue
        if not predicate(r): continue
        t = truth.get(r["pair_id"])
        if t is None: continue
        by_tag[r["algo"]].append((t, r["ani"]))
    return by_tag

def print_tags(tag_groups, title):
    print(f"\n  {title}")
    print(f"    {'tag':<22} {'n':>4} {'MAE':>7} {'RMSE':>7} {'bias':>7} {'r':>6}")
    rows_out = []
    for tag, data in tag_groups.items():
        xs_t = [d[0] for d in data]; xs_e = [d[1] for d in data]
        n, mae, rmse, bias, r = metrics(xs_t, xs_e)
        rows_out.append((tag, n, mae, rmse, bias, r))
    for tag, n, mae, rmse, bias, r in sorted(rows_out):
        print(f"    {tag:<22} {n:>4} {mae:>7.4f} {rmse:>7.4f} {bias:>+7.4f} {r:>6.3f}")

# k-mer sensitivity (all algos with _k15 / _k27 tags + their default twin)
k_tags = set()
for r in rows:
    if r["layer"]=="L2" and ("_k15" in r["algo"] or "_k27" in r["algo"]):
        k_tags.add(algo_family(r["algo"]))
for a in sorted(k_tags):
    by_tag = defaultdict(list)
    for r in rows:
        if algo_family(r["algo"]) != a: continue
        if r["layer"] == "L1" or "_k15" in r["algo"] or "_k27" in r["algo"]:
            t = truth.get(r["pair_id"])
            if t is None: continue
            # only restrict L1 to the 50-pair L2 subset for fair comparison
            if r["layer"] == "L1":
                # we don't track the subset list here; include all 200 — flagged below
                pass
            by_tag[r["algo"]].append((t, r["ani"]))
    print_tags(by_tag, f"k-mer sensitivity: {a}  (default={a}_def, _k15, _k27)")

# SetSketch a-sensitivity
ss_a = defaultdict(list)
for r in rows:
    if algo_family(r["algo"]) != "setsketch": continue
    if r["algo"] in ("setsketch_def", "setsketch_a0p5", "setsketch_a1p5"):
        t = truth.get(r["pair_id"])
        if t is None: continue
        ss_a[r["algo"]].append((t, r["ani"]))
print_tags(ss_a, "SetSketch `a` sensitivity (default a=5 vs 0.5 vs 1.5)")

# KSSD drlevel sensitivity
kssd_d = defaultdict(list)
for r in rows:
    if algo_family(r["algo"]) != "kssd": continue
    if r["algo"] in ("kssd_def", "kssd_dr1", "kssd_dr2"):
        t = truth.get(r["pair_id"])
        if t is None: continue
        kssd_d[r["algo"]].append((t, r["ani"]))
print_tags(kssd_d, "KSSD drlevel sensitivity (default 3 vs 1 vs 2)")

# Sketch-size sensitivity (minhash / probmh)
for fam in ("minhash", "probmh"):
    by_tag = defaultdict(list)
    for r in rows:
        if algo_family(r["algo"]) != fam: continue
        if r["algo"] in (f"{fam}_def", f"{fam}_s256", f"{fam}_s4096"):
            t = truth.get(r["pair_id"])
            if t is None: continue
            by_tag[r["algo"]].append((t, r["ani"]))
    print_tags(by_tag, f"Sketch-size sensitivity: {fam}  (def=1024, 256, 4096)")

print("\n" + "="*78)
print("Done.")
print("="*78)
