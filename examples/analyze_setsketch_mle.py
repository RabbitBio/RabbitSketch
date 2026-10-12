#!/usr/bin/env python3
"""Compare SetSketch inclusion-exclusion vs joint-MLE estimator.

Reads results/alisim_setsketch_mle.csv produced by run_setsketch_mle_eval.sh
(both `setsketch` and `setsketch_mle` rows are mirrored on identical params),
then emits a self-contained text report covering:

  Q1  Overall MAE/RMSE/bias across all combos
  Q2  Default-parameter (k=21, b=13, a=5, base=2) vs everything else
  Q3  By branch length: who wins where
  Q4  By sequence length: short→long scaling
  Q5  By bits: does MLE fix the "MAE↑ with bits" anomaly?
  Q6  By a parameter: which a is now best for MLE?
  Q7  By base parameter
  Q8  By k-mer size
  Q9  Bias by true J
  Q10 Per-model robustness
  Q11 ANI translation behaviour
"""
from __future__ import annotations
import csv
import math
import statistics as st
from collections import defaultdict


CSV_PATH = "results/alisim_setsketch_mle.csv"


def load(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                r["BL"]     = float(r["branch_len"])
                r["kmer"]   = int(r["kmer"])
                r["bits"]   = int(r["bits"])
                r["ss_a"]   = float(r["ss_a"])
                r["ss_base"]= float(r["ss_base"])
                r["eJ"]     = float(r["mean_exact_J"])
                r["estJ"]   = float(r["mean_est_J"])
                r["MAE"]    = float(r["MAE"])
                r["RMSE"]   = float(r["RMSE"]) if r["RMSE"] not in ("", None) else float("nan")
                r["bias"]   = r["estJ"] - r["eJ"]
                rows.append(r)
            except (ValueError, KeyError):
                continue
    return rows


def L_from_dir(d):
    for tag in ("L1k","L10k","L50k","L100k","L500k","L1M","L2M","L4M","L8M"):
        if "/fixed_%s/" % tag in d:
            return tag
    return "?"


def model_from_dir(d):
    return d.split("/")[-3] if "/fixed_" in d else "?"


def is_default(r):
    return (r["kmer"]==21 and r["bits"]==13
            and abs(r["ss_a"]-5.0)<0.1 and abs(r["ss_base"]-2.0)<0.05)


def pair_rows(rows):
    """Pair (setsketch, setsketch_mle) rows on identical (dir, k, b, a, base)."""
    key = lambda r: (r["dir"], r["kmer"], r["bits"], r["ss_a"], r["ss_base"])
    inc = {key(r): r for r in rows if r["algo"]=="setsketch"}
    mle = {key(r): r for r in rows if r["algo"]=="setsketch_mle"}
    pairs = []
    for k, ri in inc.items():
        if k in mle:
            pairs.append((ri, mle[k]))
    return pairs


def mean(xs):
    xs = list(xs)
    return sum(xs)/len(xs) if xs else float("nan")


def fmtimprov(a, b):
    if b <= 1e-12: return "inf"
    return "%.2fx" % (a/b)


# ─────────────────────────────────────────────────────────────────────────────
def main():
    rows = load(CSV_PATH)
    print(f"Loaded {len(rows)} rows from {CSV_PATH}")
    print(f"  setsketch (incl-excl)   : {sum(1 for r in rows if r['algo']=='setsketch')}")
    print(f"  setsketch_mle (joint-MLE): {sum(1 for r in rows if r['algo']=='setsketch_mle')}")

    pairs = pair_rows(rows)
    print(f"  paired combos for comparison: {len(pairs)}")
    print()

    # ── Q1: overall ─────────────────────────────────────────────────────────
    print("=" * 78)
    print("Q1. Overall accuracy across ALL parameter combos & datasets")
    print("=" * 78)
    inc_mae = [a["MAE"] for a,_ in pairs]
    mle_mae = [b["MAE"] for _,b in pairs]
    inc_bias= [a["bias"] for a,_ in pairs]
    mle_bias= [b["bias"] for _,b in pairs]
    print(f"  metric           incl-excl       MLE      MLE wins")
    print(f"  mean MAE         {mean(inc_mae):.5f}       {mean(mle_mae):.5f}    {fmtimprov(mean(inc_mae), mean(mle_mae))}")
    print(f"  median MAE       {st.median(inc_mae):.5f}       {st.median(mle_mae):.5f}    {fmtimprov(st.median(inc_mae), st.median(mle_mae))}")
    print(f"  mean bias        {mean(inc_bias):+.5f}      {mean(mle_bias):+.5f}")
    wins = sum(1 for a,b in pairs if b["MAE"] < a["MAE"])
    ties = sum(1 for a,b in pairs if abs(a["MAE"]-b["MAE"]) < 1e-7)
    print(f"  combos MLE strictly better : {wins}/{len(pairs)} ({100*wins/len(pairs):.1f}%)")
    print(f"  combos tied (within 1e-7) : {ties}/{len(pairs)}")

    # ── Q2: default vs other ────────────────────────────────────────────────
    print()
    print("=" * 78)
    print("Q2. Default-parameter comparison (k=21, bits=13, a=5, base=2)")
    print("=" * 78)
    def_pairs = [(a,b) for a,b in pairs if is_default(a)]
    inc_mae = [a["MAE"] for a,_ in def_pairs]
    mle_mae = [b["MAE"] for _,b in def_pairs]
    inc_bias= [a["bias"] for a,_ in def_pairs]
    mle_bias= [b["bias"] for _,b in def_pairs]
    print(f"  combos (default params, all BL × L × model): {len(def_pairs)}")
    print(f"  metric           incl-excl       MLE      MLE wins")
    print(f"  mean MAE         {mean(inc_mae):.5f}       {mean(mle_mae):.5f}    {fmtimprov(mean(inc_mae), mean(mle_mae))}")
    print(f"  mean bias        {mean(inc_bias):+.5f}      {mean(mle_bias):+.5f}")

    # ── Q3: by branch length ────────────────────────────────────────────────
    print()
    print("=" * 78)
    print("Q3. Default params, MAE by branch length (truer J → lower BL)")
    print("=" * 78)
    print(f"  {'BL':>8}  {'true J':>8}  {'incl MAE':>10}  {'MLE MAE':>10}  {'gain':>6}")
    BLs = sorted({a["BL"] for a,_ in def_pairs})
    for bl in BLs:
        ip = [(a,b) for a,b in def_pairs if abs(a["BL"]-bl)<1e-6]
        if not ip: continue
        tj = mean(a["eJ"] for a,_ in ip)
        mi = mean(a["MAE"] for a,_ in ip)
        mm = mean(b["MAE"] for _,b in ip)
        print(f"  {bl:>8.4f}  {tj:>8.4f}  {mi:>10.5f}  {mm:>10.5f}  {fmtimprov(mi,mm):>6}")

    # ── Q4: by sequence length ──────────────────────────────────────────────
    print()
    print("=" * 78)
    print("Q4. Default params, MAE by sequence length L")
    print("=" * 78)
    print(f"  {'L':>6}  {'incl MAE':>10}  {'MLE MAE':>10}  {'gain':>6}")
    for L in ("L1k","L10k","L50k","L100k","L500k","L1M","L2M","L4M","L8M"):
        ip = [(a,b) for a,b in def_pairs if L_from_dir(a["dir"])==L]
        if not ip: continue
        mi = mean(a["MAE"] for a,_ in ip)
        mm = mean(b["MAE"] for _,b in ip)
        print(f"  {L:>6}  {mi:>10.5f}  {mm:>10.5f}  {fmtimprov(mi,mm):>6}")

    # ── Q5: by bits — does anomaly disappear? ───────────────────────────────
    print()
    print("=" * 78)
    print("Q5. MAE by `bits` at k=21, a=5, base=2  (was: incl MAE↑ with bits)")
    print("=" * 78)
    print(f"  {'bits':>4}  {'regs':>7}  {'incl MAE':>10}  {'MLE MAE':>10}  {'gain':>6}")
    bits_pairs = [(a,b) for a,b in pairs
                  if a["kmer"]==21 and abs(a["ss_a"]-5.0)<0.1 and abs(a["ss_base"]-2)<0.05]
    bits_set = sorted({a["bits"] for a,_ in bits_pairs})
    for bts in bits_set:
        ip = [(a,b) for a,b in bits_pairs if a["bits"]==bts]
        mi = mean(a["MAE"] for a,_ in ip)
        mm = mean(b["MAE"] for _,b in ip)
        print(f"  {bts:>4}  {1<<bts:>7}  {mi:>10.5f}  {mm:>10.5f}  {fmtimprov(mi,mm):>6}")

    # ── Q6: by a parameter ──────────────────────────────────────────────────
    print()
    print("=" * 78)
    print("Q6. MAE by `a` parameter at k=21, bits=13, base=2")
    print("=" * 78)
    print(f"  {'a':>6}  {'incl MAE':>10}  {'MLE MAE':>10}  {'gain':>6}")
    a_pairs = [(a,b) for a,b in pairs
               if a["kmer"]==21 and a["bits"]==13 and abs(a["ss_base"]-2)<0.05]
    a_set = sorted({a["ss_a"] for a,_ in a_pairs})
    for av in a_set:
        ip = [(a,b) for a,b in a_pairs if abs(a["ss_a"]-av)<1e-6]
        mi = mean(a["MAE"] for a,_ in ip)
        mm = mean(b["MAE"] for _,b in ip)
        print(f"  {av:>6.2f}  {mi:>10.5f}  {mm:>10.5f}  {fmtimprov(mi,mm):>6}")

    # ── Q7: by base ─────────────────────────────────────────────────────────
    print()
    print("=" * 78)
    print("Q7. MAE by `base` parameter at k=21, bits=13, a=5")
    print("=" * 78)
    print(f"  {'base':>6}  {'incl MAE':>10}  {'MLE MAE':>10}  {'gain':>6}")
    b_pairs = [(a,b) for a,b in pairs
               if a["kmer"]==21 and a["bits"]==13 and abs(a["ss_a"]-5)<0.1]
    base_set = sorted({a["ss_base"] for a,_ in b_pairs})
    for bv in base_set:
        ip = [(a,b) for a,b in b_pairs if abs(a["ss_base"]-bv)<1e-6]
        mi = mean(a["MAE"] for a,_ in ip)
        mm = mean(b["MAE"] for _,b in ip)
        print(f"  {bv:>6.2f}  {mi:>10.5f}  {mm:>10.5f}  {fmtimprov(mi,mm):>6}")

    # ── Q8: by k-mer ────────────────────────────────────────────────────────
    print()
    print("=" * 78)
    print("Q8. MAE by k-mer at bits=13, a=5, base=2")
    print("=" * 78)
    print(f"  {'k':>4}  {'incl MAE':>10}  {'MLE MAE':>10}  {'gain':>6}")
    k_pairs = [(a,b) for a,b in pairs
               if a["bits"]==13 and abs(a["ss_a"]-5)<0.1 and abs(a["ss_base"]-2)<0.05]
    k_set = sorted({a["kmer"] for a,_ in k_pairs})
    for kv in k_set:
        ip = [(a,b) for a,b in k_pairs if a["kmer"]==kv]
        mi = mean(a["MAE"] for a,_ in ip)
        mm = mean(b["MAE"] for _,b in ip)
        print(f"  {kv:>4}  {mi:>10.5f}  {mm:>10.5f}  {fmtimprov(mi,mm):>6}")

    # ── Q9: bias by true J  (where was incl-excl most biased?) ──────────────
    print()
    print("=" * 78)
    print("Q9. Default params, BIAS (est − true) by true-J bin")
    print("=" * 78)
    bins = [(0.0,0.05),(0.05,0.15),(0.15,0.30),(0.30,0.50),
            (0.50,0.70),(0.70,0.85),(0.85,0.95),(0.95,1.01)]
    print(f"  {'J-bin':>14}  {'n':>4}  {'incl bias':>10}  {'MLE bias':>10}  {'incl MAE':>10}  {'MLE MAE':>10}")
    for lo,hi in bins:
        ip = [(a,b) for a,b in def_pairs if lo<=a["eJ"]<hi]
        if not ip: continue
        bi = mean(a["bias"] for a,_ in ip)
        bm = mean(b["bias"] for _,b in ip)
        mi = mean(a["MAE"] for a,_ in ip)
        mm = mean(b["MAE"] for _,b in ip)
        print(f"  [{lo:.2f},{hi:.2f})  {len(ip):>4}  {bi:>+10.5f}  {bm:>+10.5f}  {mi:>10.5f}  {mm:>10.5f}")

    # ── Q10: per-model robustness ───────────────────────────────────────────
    print()
    print("=" * 78)
    print("Q10. Default params, MAE by evolutionary model")
    print("=" * 78)
    print(f"  {'model':>10}  {'incl MAE':>10}  {'MLE MAE':>10}  {'gain':>6}")
    for m in sorted({model_from_dir(a["dir"]) for a,_ in def_pairs}):
        ip = [(a,b) for a,b in def_pairs if model_from_dir(a["dir"])==m]
        mi = mean(a["MAE"] for a,_ in ip)
        mm = mean(b["MAE"] for _,b in ip)
        print(f"  {m:>10}  {mi:>10.5f}  {mm:>10.5f}  {fmtimprov(mi,mm):>6}")

    # ── Q11: ANI behaviour at default params ────────────────────────────────
    # Mash distance d_mash = -ln(2J/(1+J)) / k; ANI ≈ 1 − d_mash
    print()
    print("=" * 78)
    print("Q11. ANI translation: estimated ANI vs true ANI (default params, k=21)")
    print("=" * 78)
    def mash_dist(j, k=21):
        if j<=0 or j>=1: return 0.0
        return -math.log(2.0*j/(1.0+j))/k
    print(f"  {'BL':>8}  {'true ANI':>9}  {'incl ANI':>9}  {'MLE ANI':>9}  {'incl Δ':>8}  {'MLE Δ':>8}")
    for bl in BLs:
        ip = [(a,b) for a,b in def_pairs if abs(a["BL"]-bl)<1e-6]
        if not ip: continue
        tj  = mean(a["eJ"] for a,_ in ip)
        eji = mean(a["estJ"] for a,_ in ip)
        ejm = mean(b["estJ"] for _,b in ip)
        ani_t  = 1.0 - mash_dist(tj)
        ani_i  = 1.0 - mash_dist(eji)
        ani_m  = 1.0 - mash_dist(ejm)
        print(f"  {bl:>8.4f}  {ani_t:>9.5f}  {ani_i:>9.5f}  {ani_m:>9.5f}  {ani_i-ani_t:>+8.5f}  {ani_m-ani_t:>+8.5f}")

    print()
    print("=" * 78)
    print("Done.")
    print("=" * 78)


if __name__ == "__main__":
    main()
