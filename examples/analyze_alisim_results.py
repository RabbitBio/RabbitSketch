#!/usr/bin/env python3
"""
Comprehensive analysis of results/alisim_sweep.csv.

Sections:
  A. Data inventory
  B. Default-parameter accuracy across models (algo ranking)
  C. Model sensitivity – does HKY/GTR/AT-rich change ranking?
  D. Sequence-length effect on accuracy
  E. k-mer size sensitivity (per-algo curves)
  F. Sketch size sensitivity
  G. SetSketch `a` × L cross-tab (the "low fill" problem)
  H. BinDash bbits/s64 trade-off
  I. Bias direction (positive = over-estimate)
  J. Robustness across models (within-algo stddev across models)
"""

import csv, math, statistics
from collections import defaultdict

PATH = "results/alisim_sweep.csv"

rows = []
with open(PATH) as f:
    r = csv.DictReader(f)
    for row in r:
        try:
            rmse_v = row.get("RMSE")
            rmse_f = float(rmse_v) if rmse_v not in (None, "", "None") else None
            rows.append({
                "dir":     row["dir"],
                "model":   row["model"],
                "bl":      float(row["branch_len"]),
                "algo":    row["algo"],
                "kmer":    int(row["kmer"]),
                "size":    int(row["size"]),
                "bits":    int(row["bits"]),
                "ss_a":    float(row["ss_a"]),
                "ss_base": float(row["ss_base"]),
                "bd_bbits":int(row["bd_bbits"]),
                "bd_s64":  int(row["bd_s64"]),
                "drlevel": int(row["drlevel"]),
                "n":       int(row["n"]),
                "mean_ej": float(row["mean_exact_J"]),
                "mean_est":float(row["mean_est_J"]),
                "mae":     float(row["MAE"]),
                "rmse":    rmse_f,
            })
        except Exception:
            continue

# Parse dir → (model, mode, L_int)
def parse_L(d):
    parts = d.strip("/").split("/")
    md_L = parts[-2]  # fixed_L500k or varlen_L500k
    mode, Ls = md_L.split("_", 1)
    s = Ls.lstrip("L")
    if   s.endswith("M"): L = int(float(s[:-1]) * 1_000_000)
    elif s.endswith("k"): L = int(float(s[:-1]) * 1_000)
    else: L = int(s)
    return mode, L, Ls

for r in rows:
    r["mode"], r["L"], r["Lstr"] = parse_L(r["dir"])

def is_default(r):
    a = r["algo"]
    if a == "minhash":   return r["kmer"]==21 and r["size"]==1024
    if a == "fastkmv":   return r["kmer"]==21 and r["size"]==1024
    if a == "hll":       return r["kmer"]==21 and r["bits"]==13
    if a == "setsketch": return (r["kmer"]==21 and r["bits"]==13
                                  and r["ss_a"]==5.0 and r["ss_base"]==2.0)
    if a == "bindash":   return r["kmer"]==21 and r["bd_bbits"]==16 and r["bd_s64"]==16
    if a == "kssd":      return r["drlevel"]==3
    if a == "probmh":    return r["kmer"]==21 and r["size"]==1024
    return False

ALGOS = ["minhash","fastkmv","hll","setsketch","bindash","kssd","probmh"]

def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs)/len(xs) if xs else None

def std(xs):
    xs = [x for x in xs if x is not None]
    return statistics.pstdev(xs) if len(xs)>=2 else None

# ════════════════════════════════════════════════════════════════════════════
# A. Inventory
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("A. Data inventory")
print("="*80)
models = sorted(set(r["model"] for r in rows))
Ls     = sorted(set(r["L"] for r in rows if r["mode"]=="fixed"))
bls    = sorted(set(r["bl"] for r in rows))
print(f"  Models  ({len(models)}):  {', '.join(models)}")
print(f"  Lengths ({len(Ls)}):     {Ls}")
print(f"  BLs     ({len(bls)}):     {bls}")
print(f"  Rows: {len(rows)}, fixed: {sum(1 for r in rows if r['mode']=='fixed')}, "
      f"varlen: {sum(1 for r in rows if r['mode']=='varlen')}")
print()

# ════════════════════════════════════════════════════════════════════════════
# B. Default-parameter ranking, averaged over all combos & BLs
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("B. Default-parameter MAE by algorithm (averaged over all model/L/BL)")
print("="*80)
def_rows = [r for r in rows if is_default(r) and r["mode"]=="fixed"]
print(f"\n   algo         mean_MAE   median_MAE   n_combos    rank")
ranked = []
for a in ALGOS:
    maes = [r["mae"] for r in def_rows if r["algo"]==a]
    if not maes: continue
    m  = sum(maes)/len(maes)
    md = statistics.median(maes)
    ranked.append((m, md, len(maes), a))
ranked.sort()
for rank,(m,md,n,a) in enumerate(ranked,1):
    print(f"   {a:<12}  {m:.4f}     {md:.4f}      {n:>4}        #{rank}")
print()

# ════════════════════════════════════════════════════════════════════════════
# C. Model sensitivity at L=100k, BL=0.05
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("C. Model sensitivity: MAE at L=100k, BL=0.05, default params")
print("   (Does HKY / GTR+G / AT-rich change ranking?)")
print("="*80)
sub = [r for r in def_rows if r["L"]==100_000 and abs(r["bl"]-0.05)<1e-6]
print(f"\n   {'algo':<10} " + "".join(f"{m[:8]:>9} " for m in models))
for a in ALGOS:
    cells=[]
    for m in models:
        row=[r for r in sub if r["algo"]==a and r["model"]==m]
        cells.append(f"{row[0]['mae']:>9.4f} " if row else "    -     ")
    print(f"   {a:<10} " + "".join(cells))
print()

# ════════════════════════════════════════════════════════════════════════════
# D. Length effect: averaged over models, at BL=0.05
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("D. Sequence-length effect (BL=0.05, default params, avg over 8 models)")
print("="*80)
sub = [r for r in def_rows if abs(r["bl"]-0.05)<1e-6]
print(f"\n   {'algo':<10} " + "".join(f"L={L:>9,} " for L in Ls))
for a in ALGOS:
    cells=[]
    for L in Ls:
        ms = [r["mae"] for r in sub if r["algo"]==a and r["L"]==L]
        cells.append(f"{mean(ms):>9.4f}   " if ms else "      -      ")
    print(f"   {a:<10} " + "".join(cells))
print()

# ════════════════════════════════════════════════════════════════════════════
# E. k-mer sensitivity: avg over models, at L=100k BL=0.05
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("E. k-mer sweep at L=100k, BL=0.05 (avg over 8 substitution models)")
print("="*80)
sub = [r for r in rows if r["mode"]=="fixed" and r["L"]==100_000
       and abs(r["bl"]-0.05)<1e-6]
def kmer_sweep(algo):
    out=[]
    for r in sub:
        if r["algo"]!=algo: continue
        ok=False
        if algo=="minhash"   and r["size"]==1024: ok=True
        if algo=="fastkmv"   and r["size"]==1024: ok=True
        if algo=="hll"       and r["bits"]==13:   ok=True
        if algo=="setsketch" and r["bits"]==13 and r["ss_a"]==5 and r["ss_base"]==2: ok=True
        if algo=="bindash"   and r["bd_bbits"]==16 and r["bd_s64"]==16: ok=True
        if algo=="probmh"    and r["size"]==1024: ok=True
        if ok: out.append(r)
    return out

for a in ["minhash","fastkmv","probmh","bindash","hll","setsketch"]:
    sw = kmer_sweep(a)
    by_k = defaultdict(list)
    for r in sw: by_k[r["kmer"]].append(r["mae"])
    ks = sorted(by_k.keys())
    if not ks: continue
    print(f"\n   {a}:")
    print(f"     k   = " + "".join(f"{k:>6d} " for k in ks))
    print(f"     MAE = " + "".join(f"{mean(by_k[k]):>6.4f} " for k in ks))
print()

# ════════════════════════════════════════════════════════════════════════════
# F. Sketch size sensitivity (BL=0.05, L=100k)
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("F. Sketch size / register sensitivity, L=100k, BL=0.05")
print("="*80)
sub = [r for r in rows if r["mode"]=="fixed" and r["L"]==100_000
       and abs(r["bl"]-0.05)<1e-6]
for algo, param, filt in [
    ("minhash","size",   lambda r: r["kmer"]==21),
    ("fastkmv","size",   lambda r: r["kmer"]==21),
    ("probmh","size",    lambda r: r["kmer"]==21),
    ("hll","bits",       lambda r: r["kmer"]==21),
    ("setsketch","bits", lambda r: r["kmer"]==21 and r["ss_a"]==5 and r["ss_base"]==2),
    ("bindash","bd_s64", lambda r: r["kmer"]==21 and r["bd_bbits"]==16),
]:
    sw = [r for r in sub if r["algo"]==algo and filt(r)]
    by_p = defaultdict(list)
    for r in sw: by_p[r[param]].append(r["mae"])
    ps = sorted(by_p.keys())
    if not ps: continue
    print(f"\n   {algo} (varying {param}):")
    print(f"     {param:<6} = " + "".join(f"{p:>6} " for p in ps))
    print(f"     MAE    = " + "".join(f"{mean(by_p[p]):>6.4f} " for p in ps))
print()

# ════════════════════════════════════════════════════════════════════════════
# G. SetSketch a × L cross-tab (THE big finding)
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("G. SetSketch `a` × Length cross-tab (BL=0.05, bits=13, base=2, kmer=21)")
print("   Shows how `a` must be tuned per L (low fill ↔ small a)")
print("="*80)
sub = [r for r in rows if r["mode"]=="fixed" and abs(r["bl"]-0.05)<1e-6
       and r["algo"]=="setsketch" and r["kmer"]==21
       and r["bits"]==13 and r["ss_base"]==2.0]
# average over models
by_aL = defaultdict(list)
for r in sub: by_aL[(r["ss_a"], r["L"])].append(r["mae"])
a_vals = sorted(set(k[0] for k in by_aL.keys()))
L_vals = sorted(set(k[1] for k in by_aL.keys()))
print(f"\n   {'a\\L':<7} " + "".join(f"{L:>9,} " for L in L_vals))
for a in a_vals:
    cells=[]
    for L in L_vals:
        ms = by_aL.get((a,L))
        cells.append(f"{mean(ms):>9.4f} " if ms else "    -      ")
    print(f"   {a:<7} " + "".join(cells))
print()

# Best a per L
print("\n   Best `a` for each L (lowest MAE):")
print(f"   {'L':>10}   best_a   MAE_at_best   MAE_at_default(a=5)   ratio")
for L in L_vals:
    best=(None, float('inf'))
    a_default=None
    for a in a_vals:
        ms = by_aL.get((a,L))
        if ms is None: continue
        m = mean(ms)
        if m < best[1]: best=(a,m)
        if a==5.0: a_default=m
    if best[0] is None: continue
    ratio = a_default/best[1] if (a_default and best[1]>0) else None
    print(f"   {L:>10,}   {best[0]:>4}    {best[1]:.4f}        "
          f"{a_default:.4f}              {ratio:.1f}×" if ratio else
          f"   {L:>10,}   {best[0]:>4}    {best[1]:.4f}")
print()

# ════════════════════════════════════════════════════════════════════════════
# H. BinDash bbits sensitivity (L=100k, BL=0.05)
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("H. BinDash bbits & s64 trade-off, L=100k, BL=0.05")
print("="*80)
sub = [r for r in rows if r["mode"]=="fixed" and r["L"]==100_000
       and abs(r["bl"]-0.05)<1e-6 and r["algo"]=="bindash" and r["kmer"]==21]

bb = defaultdict(list)
for r in sub:
    if r["bd_s64"]==16: bb[r["bd_bbits"]].append(r["mae"])
bbs=sorted(bb.keys())
print(f"\n   bbits sweep (s64=16):")
print(f"     bb = " + "".join(f"{b:>6} " for b in bbs))
print(f"     MAE= " + "".join(f"{mean(bb[b]):>6.4f} " for b in bbs))

s64 = defaultdict(list)
for r in sub:
    if r["bd_bbits"]==16: s64[r["bd_s64"]].append(r["mae"])
s64s = sorted(s64.keys())
print(f"\n   s64 sweep (bbits=16):  # regs = s64×64")
print(f"     s64= " + "".join(f"{s:>6} " for s in s64s))
print(f"     reg= " + "".join(f"{s*64:>6} " for s in s64s))
print(f"     MAE= " + "".join(f"{mean(s64[s]):>6.4f} " for s in s64s))
print()

# ════════════════════════════════════════════════════════════════════════════
# I. Bias direction at default params
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("I. Bias = mean_est − mean_exact at L=100k (avg over models)")
print("   Positive = over-estimate; near 0 = unbiased")
print("="*80)
sub = [r for r in def_rows if r["L"]==100_000]
bls_here = sorted(set(r["bl"] for r in sub))
print(f"\n   {'algo':<10} " + "".join(f"BL={b:.4f} " for b in bls_here))
for a in ALGOS:
    cells=[]
    for b in bls_here:
        rs=[r for r in sub if r["algo"]==a and r["bl"]==b]
        if rs:
            biases=[r["mean_est"]-r["mean_ej"] for r in rs]
            cells.append(f"  {mean(biases):>+7.4f} ")
        else:
            cells.append("    -      ")
    print(f"   {a:<10} " + "".join(cells))
print()

# ════════════════════════════════════════════════════════════════════════════
# J. Cross-model stability (default params, BL=0.05, L=100k)
# ════════════════════════════════════════════════════════════════════════════
print("="*80)
print("J. Cross-model stability: std of MAE across 8 models")
print("   (Smaller = algo less affected by evolutionary model)")
print("="*80)
sub = [r for r in def_rows if r["L"]==100_000 and abs(r["bl"]-0.05)<1e-6]
print(f"\n   {'algo':<10}   mean_MAE    std_MAE    cv%")
for a in ALGOS:
    rs=[r for r in sub if r["algo"]==a]
    if len(rs)<2: continue
    maes=[r["mae"] for r in rs]
    m=mean(maes); s=std(maes)
    cv=100*s/m if m else 0
    print(f"   {a:<10}   {m:.4f}      {s:.4f}     {cv:.1f}%")
print()
