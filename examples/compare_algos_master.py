#!/usr/bin/env python3
"""Cross-algorithm comparison on alisim_master.csv (SetSketch now = joint-MLE).

Only uses rows with is_default=1, so we are comparing each algorithm's
*recommended* configuration, not their cherry-picked best parameters.
"""
from __future__ import annotations
import csv
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

MASTER = Path("results/alisim_master.csv")
ALGOS  = ["minhash","fastkmv","hll","setsketch","bindash","kssd","probmh"]

COMMON_LS = {1000, 10000, 50000, 100000, 500000}   # other algos only cover these

def load(path, restrict_L=True):
    rows=[]
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                if int(r["is_default"]) != 1: continue
                r["BL"]   = float(r["BL"])
                r["L"]    = int(r["L"])
                if restrict_L and r["L"] not in COMMON_LS: continue
                r["MAE"]  = float(r["MAE_mean"])
                r["RMSE"] = float(r["RMSE_mean"]) if r["RMSE_mean"] not in ("","nan") else float("nan")
                r["bias"] = float(r["bias_mean"])
                r["eJ"]   = float(r["exactJ_mean"])
                r["estJ"] = float(r["estJ_mean"])
                rows.append(r)
            except (ValueError, KeyError):
                continue
    return rows

def mean(xs):
    xs=list(xs); return sum(xs)/len(xs) if xs else float("nan")

def fmt_row(d, header):
    out=[]
    for h in header:
        v = d.get(h, float("nan"))
        if isinstance(v, float):
            out.append("%.5f" % v if not math.isnan(v) else "  -    ")
        else:
            out.append(str(v))
    return out

def table(title, data, row_label, col_labels):
    print()
    print("=" * 80)
    print(title)
    print("=" * 80)
    w = max(8, max(len(str(c)) for c in col_labels)+1)
    print(f"  {row_label:<10}" + "".join(f"{c:>{w}}" for c in col_labels))
    for r in data:
        print(f"  {r['label']:<10}" + "".join(
            (("%9.5f" % r[c]) if isinstance(r.get(c, float("nan")), float)
                                 and not math.isnan(r.get(c, float("nan")))
             else f"{'  -':>{w}}")
            for c in col_labels))

def main():
    print("[Fair comparison: restricted to L ∈ {1k, 10k, 50k, 100k, 500k},\n",
          " the range covered by ALL 7 algorithms.\n",
          " SetSketch has additional 1M-8M data that other algos don't.]")
    rows = load(MASTER, restrict_L=True)
    print(f"Loaded {len(rows)} default-config rows")
    print(f"Per-algo row counts (default × 5 Ls × 8 BLs × n_models):")
    by_algo = defaultdict(list)
    for r in rows: by_algo[r["algo"]].append(r)
    for a in ALGOS:
        print(f"  {a:<10} {len(by_algo[a]):>4}")

    # ── overall summary
    print()
    print("=" * 80)
    print("OVERALL: default params, all (L × BL × model) datasets")
    print("(NOTE: mean MAE is dominated by L=1k where SetSketch/KSSD defaults")
    print(" don't fit; median is the robust comparison)")
    print("=" * 80)
    print(f"  {'algo':<10}  {'rows':>4}  {'mean MAE':>10}  {'med MAE':>10}  {'mean bias':>11}  {'true J':>8}")
    summary = []
    for a in ALGOS:
        rs = by_algo[a]
        if not rs: continue
        m  = mean(r["MAE"] for r in rs)
        md = st.median(r["MAE"] for r in rs)
        bi = mean(r["bias"] for r in rs)
        tj = mean(r["eJ"] for r in rs)
        summary.append((a, len(rs), m, md, bi, tj))
        print(f"  {a:<10}  {len(rs):>4}  {m:>10.5f}  {md:>10.5f}  {bi:>+11.5f}  {tj:>8.4f}")

    # Rank
    print()
    print("  Rank by median MAE (robust to L=1k outliers):")
    for i,(a,n,m,md,*_) in enumerate(sorted(summary, key=lambda x: x[3]), 1):
        print(f"    {i}. {a:<10}  median MAE={md:.5f}   (mean={m:.5f})")

    # Restricted to L ≥ 10k (params fit)
    print()
    print("=" * 80)
    print("FAIR RANK: L ∈ {10k, 50k, 100k, 500k}  (all algo defaults fit)")
    print("=" * 80)
    print(f"  {'algo':<10}  {'rows':>4}  {'mean MAE':>10}  {'med MAE':>10}  {'mean bias':>11}")
    fair_summary = []
    for a in ALGOS:
        rs = [r for r in by_algo[a] if r["L"] >= 10000]
        if not rs: continue
        m  = mean(r["MAE"] for r in rs)
        md = st.median(r["MAE"] for r in rs)
        bi = mean(r["bias"] for r in rs)
        fair_summary.append((a, len(rs), m, md, bi))
        print(f"  {a:<10}  {len(rs):>4}  {m:>10.5f}  {md:>10.5f}  {bi:>+11.5f}")

    print()
    print("  Rank by mean MAE (fair, L ≥ 10k):")
    for i,(a,n,m,md,bi) in enumerate(sorted(fair_summary, key=lambda x: x[2]), 1):
        print(f"    {i}. {a:<10}  MAE={m:.5f}")

    # ── BL breakdown
    print()
    print("=" * 80)
    print("MAE by branch length (default params, averaged over L and models)")
    print("=" * 80)
    BLs = sorted({r["BL"] for r in rows})
    print(f"  {'algo':<10}" + "".join(f"  d={bl:<6.4f}" for bl in BLs))
    for a in ALGOS:
        line = f"  {a:<10}"
        for bl in BLs:
            rs = [r for r in by_algo[a] if abs(r["BL"]-bl)<1e-6]
            if rs:
                line += "  %8.5f" % mean(r["MAE"] for r in rs)
            else:
                line += "     -    "
        print(line)

    # ── L breakdown
    print()
    print("=" * 80)
    print("MAE by sequence length L (default params, averaged over BL and models)")
    print("=" * 80)
    Ls = sorted({r["L"] for r in rows})
    L_tags = {1000:"1k",10000:"10k",50000:"50k",100000:"100k",
              500000:"500k",1000000:"1M",2000000:"2M",4000000:"4M",
              8000000:"8M"}
    print(f"  {'algo':<10}" + "".join(f"  L={L_tags.get(L,'?'):<5}" for L in Ls))
    for a in ALGOS:
        line = f"  {a:<10}"
        for L in Ls:
            rs = [r for r in by_algo[a] if r["L"]==L]
            if rs:
                line += "  %8.5f" % mean(r["MAE"] for r in rs)
            else:
                line += "     -    "
        print(line)

    # ── J-bin breakdown
    print()
    print("=" * 80)
    print("MAE by true-J bin (default params)")
    print("=" * 80)
    bins=[(0.95,1.01),(0.85,0.95),(0.70,0.85),(0.50,0.70),
          (0.30,0.50),(0.15,0.30),(0.05,0.15),(0.0,0.05)]
    print(f"  {'algo':<10}" + "".join(f"  [{lo:.2f},{hi:.2f})" for lo,hi in bins))
    for a in ALGOS:
        line = f"  {a:<10}"
        for lo,hi in bins:
            rs = [r for r in by_algo[a] if lo <= r["eJ"] < hi]
            if rs:
                line += "  %10.5f" % mean(r["MAE"] for r in rs)
            else:
                line += "      -     "
        print(line)

    # ── per-(L, BL) winner table
    print()
    print("=" * 80)
    print("WINNER per (L × BL) cell  (default params, lowest MAE)")
    print("=" * 80)
    print(f"  {'L':<6} | " + "".join(f"  BL={bl:<6.4f}" for bl in BLs))
    print("  " + "-" * (8 + 13*len(BLs)))
    for L in Ls:
        line = f"  {L_tags.get(L,'?'):<6} | "
        for bl in BLs:
            best_a, best_m = None, float("inf")
            for a in ALGOS:
                rs = [r for r in by_algo[a] if r["L"]==L and abs(r["BL"]-bl)<1e-6]
                if rs:
                    m = mean(r["MAE"] for r in rs)
                    if m < best_m:
                        best_m, best_a = m, a
            line += f"  {best_a or '-':<11}"
        print(line)

if __name__ == "__main__":
    main()
