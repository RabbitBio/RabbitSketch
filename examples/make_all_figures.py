#!/usr/bin/env python3
"""
Generate the full 12-figure set from results/alisim_master.csv.

Each figure is independent; run individual functions or `python3 make_all_figures.py`
to produce them all. Outputs go to results/figures/figNN_<slug>.{png,pdf}.

Figure inventory:
  fig01  Algorithm landscape  — MAE vs L, 7 lines, BL∈{0.005,0.05} panels
  fig02  BL × L heatmaps per algo (default params)               — 7 small heatmaps
  fig03  k-mer dependence per algo, lines = L                    — 6 small panels
  fig04  Sketch size dependence per algo, lines = L              — 5 small panels
  fig05  SetSketch a × L heatmap (BL=0.05)                       — 1 panel
  fig06  KSSD drlevel × L, BL panels                             — 3 lines × 3 panels
  fig07  Bias direction, x=BL, lines=algo (L=500k)               — 1 panel
  fig08  Fixed vs varlen, x=L per algo                           — 7 small panels
  fig09  Per-evolutionary-model robustness heatmap               — 1 panel
  fig10  BinDash bbits × s64 heatmap                             — 1 panel
  fig11  k × BL heatmaps per algo (low-vs-high mutation regime)  — 6 small panels
  fig12  SetSketch base & bits sweeps                            — 2 panels
"""

import csv, math, os
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

HERE = Path(__file__).resolve().parent
CSV  = HERE / "results/alisim_master.csv"
OUT  = HERE / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "figure.dpi":      120,
    "savefig.dpi":     180,
    "font.size":       9,
    "axes.titlesize":  10,
    "axes.labelsize":  9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 1.8,
    "lines.markersize":4,
})

# ── Algorithm colours and order ───────────────────────────────────────────────
ALGOS = ["minhash", "fastkmv", "probmh", "bindash", "hll", "setsketch", "kssd"]
ALGO_COLOR = {
    "minhash":   "#1f77b4",  # blue
    "fastkmv":   "#ff7f0e",  # orange
    "probmh":    "#2ca02c",  # green
    "bindash":   "#d62728",  # red
    "hll":       "#9467bd",  # purple
    "setsketch": "#8c564b",  # brown
    "kssd":      "#e377c2",  # pink
}

# ── Loader ────────────────────────────────────────────────────────────────────
def load():
    rows = []
    with open(CSV) as f:
        for r in csv.DictReader(f):
            try:
                rows.append({
                    "mode":     r["mode"],
                    "L":        int(r["L"]),
                    "Lstr":     r["Lstr"],
                    "BL":       float(r["BL"]),
                    "algo":     r["algo"],
                    "kmer":     int(r["kmer"]),
                    "size":     int(r["size"]),
                    "bits":     int(r["bits"]),
                    "ss_a":     float(r["ss_a"]),
                    "ss_base":  float(r["ss_base"]),
                    "bd_bbits": int(r["bd_bbits"]),
                    "bd_s64":   int(r["bd_s64"]),
                    "drlevel":  int(r["drlevel"]),
                    "is_def":   int(r["is_default"]),
                    "n_models": int(r["n_models"]),
                    "ej":       float(r["exactJ_mean"]),
                    "est":      float(r["estJ_mean"]),
                    "mae":      float(r["MAE_mean"]),
                    "mae_std":  float(r["MAE_std"]),
                    "rmse":     float(r["RMSE_mean"]) if r["RMSE_mean"] else None,
                    "bias":     float(r["bias_mean"]),
                })
            except Exception:
                continue
    print(f"loaded {len(rows)} rows from {CSV}")
    return rows

def save(name, fig, also_pdf=False):
    png = OUT / f"{name}.png"
    fig.savefig(png, bbox_inches="tight", dpi=180)
    if also_pdf:
        fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    size = os.path.getsize(png) // 1024
    print(f"  wrote {png.name}  ({size} KB)")

def L_fmt(L):
    return f"{L//1000}k" if L < 1_000_000 else f"{L//1_000_000}M" if L >= 1_000_000 else f"{L}"

# ════════════════════════════════════════════════════════════════════════════
# Fig 1 — Algorithm landscape
# ════════════════════════════════════════════════════════════════════════════
def fig01(rows):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=False)
    for ax, BL in zip(axes, [0.005, 0.05]):
        Ls = sorted(set(r["L"] for r in rows if r["mode"]=="fixed"))
        for a in ALGOS:
            ys = []
            for L in Ls:
                ms = [r["mae"] for r in rows
                      if r["mode"]=="fixed" and r["is_def"] and r["algo"]==a
                         and r["L"]==L and abs(r["BL"]-BL)<1e-7]
                ys.append(np.mean(ms) if ms else np.nan)
            ax.plot(Ls, ys, "-o", color=ALGO_COLOR[a], label=a)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Sequence length L (bp)")
        ax.set_ylabel("MAE")
        ax.set_title(f"BL = {BL}   ({'close-related' if BL<0.02 else 'moderate-divergence'})")
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xticks(Ls); ax.set_xticklabels([L_fmt(L) for L in Ls], rotation=30)
    axes[0].legend(ncol=2, loc="upper right", framealpha=0.9)
    fig.suptitle("Fig 1. Algorithm landscape: MAE vs sequence length at default parameters",
                 y=1.01, fontsize=11)
    save("fig01_algo_landscape", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 2 — BL × L heatmaps per algorithm
# ════════════════════════════════════════════════════════════════════════════
def fig02(rows):
    Ls  = sorted(set(r["L"] for r in rows if r["mode"]=="fixed"))
    BLs = sorted(set(r["BL"] for r in rows if r["mode"]=="fixed"))
    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    axes = axes.flatten()
    for ax, a in zip(axes, ALGOS):
        M = np.full((len(BLs), len(Ls)), np.nan)
        for i, bl in enumerate(BLs):
            for j, L in enumerate(Ls):
                ms = [r["mae"] for r in rows
                      if r["mode"]=="fixed" and r["is_def"] and r["algo"]==a
                         and r["L"]==L and abs(r["BL"]-bl)<1e-7]
                if ms: M[i,j] = np.mean(ms)
        im = ax.imshow(M, aspect="auto", origin="lower",
                       norm=mcolors.LogNorm(vmin=0.0005, vmax=1.0),
                       cmap="YlOrRd")
        ax.set_xticks(range(len(Ls)))
        ax.set_xticklabels([L_fmt(L) for L in Ls], rotation=45, fontsize=7)
        ax.set_yticks(range(len(BLs)))
        ax.set_yticklabels([f"{b:.4f}" for b in BLs], fontsize=7)
        ax.set_title(a)
        for i in range(len(BLs)):
            for j in range(len(Ls)):
                v = M[i,j]
                if not np.isnan(v):
                    txt = f"{v:.3f}" if v >= 0.001 else f"{v:.0e}"
                    ax.text(j, i, txt, ha="center", va="center",
                            color="black" if v<0.05 else "white", fontsize=5.5)
    axes[-1].axis("off")  # 8th cell unused
    cbar = fig.colorbar(im, ax=axes.tolist(), shrink=0.7, pad=0.02)
    cbar.set_label("MAE (log)")
    fig.suptitle("Fig 2. BL × L MAE heatmap per algorithm (default parameters)",
                 fontsize=11, y=0.99)
    save("fig02_heatmap_bl_L", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 3 — k-mer dependence per algorithm, lines = L
# ════════════════════════════════════════════════════════════════════════════
def fig03(rows):
    Ls_pick = [10_000, 100_000, 1_000_000, 8_000_000]
    Lcolors = ["#bbbbbb", "#888888", "#444444", "#000000"]
    BL = 0.05
    algos = ["minhash","fastkmv","probmh","bindash","hll","setsketch"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, a in zip(axes, algos):
        for L, c in zip(Ls_pick, Lcolors):
            ks=[]; ms=[]
            cand = [r for r in rows
                    if r["mode"]=="fixed" and r["algo"]==a and r["L"]==L
                       and abs(r["BL"]-BL)<1e-7]
            # 1-D sweep: vary kmer at default-other-params
            if a in ("minhash","fastkmv","probmh"):
                sub=[r for r in cand if r["size"]==1024]
            elif a in ("hll",):
                sub=[r for r in cand if r["bits"]==13]
            elif a=="setsketch":
                sub=[r for r in cand if r["bits"]==13 and r["ss_a"]==5.0 and r["ss_base"]==2.0]
            elif a=="bindash":
                sub=[r for r in cand if r["bd_bbits"]==16 and r["bd_s64"]==16]
            by_k=defaultdict(list)
            for r in sub: by_k[r["kmer"]].append(r["mae"])
            ks=sorted(by_k.keys())
            ms=[np.mean(by_k[k]) for k in ks]
            if ks:
                ax.plot(ks, ms, "-o", color=c, label=f"L={L_fmt(L)}")
        ax.set_yscale("log")
        ax.set_title(a)
        ax.grid(True, alpha=0.3)
        ax.set_xticks([11,15,19,21,23,27,31])
    for ax in axes[:3]: ax.set_xlabel("")
    for ax in axes[3:]: ax.set_xlabel("k-mer size")
    for ax in axes[::3]: ax.set_ylabel("MAE")
    axes[0].legend(loc="best", fontsize=7)
    fig.suptitle(f"Fig 3. k-mer effect at BL = {BL} (avg over 11 models)", fontsize=11, y=1.00)
    save("fig03_kmer_effect", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 4 — Sketch-size dependence per algorithm
# ════════════════════════════════════════════════════════════════════════════
def fig04(rows):
    Ls_pick = [10_000, 100_000, 1_000_000, 8_000_000]
    Lcolors = ["#bbbbbb", "#888888", "#444444", "#000000"]
    BL = 0.05
    specs = [
        ("minhash",   "size",     "sketch size",
            lambda r: r["kmer"]==21),
        ("fastkmv",   "size",     "sketch size",
            lambda r: r["kmer"]==21),
        ("probmh",    "size",     "sketch size",
            lambda r: r["kmer"]==21),
        ("hll",       "bits",     "log2(registers)",
            lambda r: r["kmer"]==21),
        ("setsketch", "bits",     "log2(registers)",
            lambda r: r["kmer"]==21 and r["ss_a"]==5.0 and r["ss_base"]==2.0),
        ("bindash",   "bd_s64",   "s64 (×64 regs)",
            lambda r: r["kmer"]==21 and r["bd_bbits"]==16),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey=True)
    axes = axes.flatten()
    for ax, (a, pkey, xlabel, filt) in zip(axes, specs):
        for L, c in zip(Ls_pick, Lcolors):
            sub=[r for r in rows
                 if r["mode"]=="fixed" and r["algo"]==a and r["L"]==L
                    and abs(r["BL"]-BL)<1e-7 and filt(r)]
            by_p=defaultdict(list)
            for r in sub: by_p[r[pkey]].append(r["mae"])
            pp=sorted(by_p.keys())
            if not pp: continue
            ms=[np.mean(by_p[p]) for p in pp]
            ax.plot(pp, ms, "-o", color=c, label=f"L={L_fmt(L)}")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(a)
        ax.set_xlabel(xlabel)
        ax.grid(True, which="both", alpha=0.3)
    for ax in axes[::3]: ax.set_ylabel("MAE")
    axes[0].legend(loc="best", fontsize=7)
    fig.suptitle(f"Fig 4. Sketch-size effect at BL = {BL}, k = 21 (avg over 11 models)",
                 fontsize=11, y=1.00)
    save("fig04_size_effect", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 5 — SetSketch a × L heatmap
# ════════════════════════════════════════════════════════════════════════════
def fig05(rows):
    BL = 0.05
    sub = [r for r in rows if r["mode"]=="fixed" and r["algo"]=="setsketch"
                              and r["kmer"]==21 and r["bits"]==13
                              and r["ss_base"]==2.0 and abs(r["BL"]-BL)<1e-7]
    aL = defaultdict(list)
    for r in sub: aL[(r["ss_a"], r["L"])].append(r["mae"])
    Ls = sorted(set(L for _,L in aL))
    avals = sorted(set(a for a,_ in aL))
    M = np.full((len(avals), len(Ls)), np.nan)
    for i,a in enumerate(avals):
        for j,L in enumerate(Ls):
            if (a,L) in aL: M[i,j] = np.mean(aL[(a,L)])
    fig, ax = plt.subplots(figsize=(9.5, 6))
    im = ax.imshow(M, aspect="auto", origin="lower",
                   norm=mcolors.LogNorm(vmin=0.005, vmax=1.0), cmap="YlOrRd")
    ax.set_xticks(range(len(Ls))); ax.set_xticklabels([L_fmt(L) for L in Ls], rotation=30)
    ax.set_yticks(range(len(avals))); ax.set_yticklabels([str(a) for a in avals])
    ax.set_xlabel("Sequence length L"); ax.set_ylabel("SetSketch  a")
    # highlight default a=5
    if 5.0 in avals:
        ax.axhline(avals.index(5.0), color="cyan", lw=1.0, ls="--", alpha=0.6)
        ax.text(-0.5, avals.index(5.0), "default", color="cyan", va="center",
                fontsize=8, ha="right")
    # annotate best-a per L
    best_idx = np.nanargmin(M, axis=0)
    for j,b in enumerate(best_idx):
        ax.scatter([j],[b], color="black", s=30, marker="x")
    for i in range(len(avals)):
        for j in range(len(Ls)):
            v = M[i,j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2g}",
                        ha="center", va="center", fontsize=6,
                        color="black" if v<0.1 else "white")
    fig.colorbar(im, ax=ax, label="MAE (log)").ax.yaxis.set_label_position("left")
    ax.set_title(f"Fig 5. SetSketch a × L  (BL={BL}, bits=13, base=2, k=21)"
                 "\nCyan dashed = default a=5;  ✕ = optimal a per L")
    save("fig05_setsketch_aL", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 6 — KSSD drlevel × L
# ════════════════════════════════════════════════════════════════════════════
def fig06(rows):
    BLs_panels = [0.005, 0.05, 0.15]
    Ls = sorted(set(r["L"] for r in rows if r["mode"]=="fixed"))
    drcolor = {1:"#2ca02c", 2:"#ff7f0e", 3:"#d62728"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), sharey=True)
    for ax, BL in zip(axes, BLs_panels):
        for dr in (1,2,3):
            ys=[]
            for L in Ls:
                ms = [r["mae"] for r in rows
                      if r["mode"]=="fixed" and r["algo"]=="kssd"
                         and r["drlevel"]==dr and r["L"]==L
                         and abs(r["BL"]-BL)<1e-7]
                ys.append(np.mean(ms) if ms else np.nan)
            ax.plot(Ls, ys, "-o", color=drcolor[dr],
                    label=f"drlevel={dr}{'  (default)' if dr==3 else ''}")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xticks(Ls); ax.set_xticklabels([L_fmt(L) for L in Ls], rotation=45)
        ax.set_xlabel("L"); ax.grid(True, which="both", alpha=0.3)
        ax.set_title(f"BL = {BL}")
    axes[0].set_ylabel("MAE")
    axes[0].legend(loc="upper right", framealpha=0.9)
    fig.suptitle("Fig 6. KSSD drlevel × L  (smaller drlevel = more registers, better accuracy)",
                 fontsize=11, y=1.02)
    save("fig06_kssd_drlevel", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 7 — Bias direction
# ════════════════════════════════════════════════════════════════════════════
def fig07(rows):
    L = 500_000
    BLs = sorted(set(r["BL"] for r in rows if r["mode"]=="fixed"))
    fig, ax = plt.subplots(figsize=(11, 4.5))
    x = np.arange(len(BLs))
    width = 0.11
    for k, a in enumerate(ALGOS):
        ys=[]
        for bl in BLs:
            ms=[r["bias"] for r in rows
                if r["mode"]=="fixed" and r["is_def"] and r["algo"]==a
                   and r["L"]==L and abs(r["BL"]-bl)<1e-7]
            ys.append(np.mean(ms) if ms else 0.0)
        ax.bar(x + (k-3)*width, ys, width, color=ALGO_COLOR[a], label=a)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels([f"{b:.4f}" for b in BLs], rotation=20)
    ax.set_xlabel("Branch length"); ax.set_ylabel("Bias (mean_est − mean_exact_J)")
    ax.set_title(f"Fig 7. Bias direction at default params, L = {L:,}\n"
                 "Positive = over-estimate; KSSD ≈ 0 (unbiased)")
    ax.legend(ncol=4, fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)
    save("fig07_bias", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 8 — Fixed vs Varlen, default parameters
# ════════════════════════════════════════════════════════════════════════════
def fig08(rows):
    BL = 0.05
    Ls = sorted(set(r["L"] for r in rows if r["mode"]=="fixed"))
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.5), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, a in zip(axes, ALGOS):
        for mode, ls, mk in [("fixed","-","o"), ("varlen","--","s")]:
            ys=[]
            for L in Ls:
                ms=[r["mae"] for r in rows
                    if r["mode"]==mode and r["is_def"] and r["algo"]==a
                       and r["L"]==L and abs(r["BL"]-BL)<1e-7]
                ys.append(np.mean(ms) if ms else np.nan)
            ax.plot(Ls, ys, ls+mk, color=ALGO_COLOR[a], label=mode, alpha=0.85)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(a)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xticks(Ls); ax.set_xticklabels([L_fmt(L) for L in Ls], rotation=45)
        ax.legend(loc="best", fontsize=7)
    axes[-1].axis("off")
    fig.suptitle(f"Fig 8. Fixed vs Variable-length comparison at BL = {BL} (default params)",
                 fontsize=11, y=0.99)
    save("fig08_fixed_vs_varlen", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 9 — Per-evolutionary-model heatmap (CSV stores 11-model average per row;
#   so we need per-model breakdown from raw sweep CSV).
# ════════════════════════════════════════════════════════════════════════════
def fig09():
    # use raw sweep csv (per-model), aggregate just over the L=500k BL=0.05 default rows
    RAW = HERE / "results/alisim_sweep.csv"
    if not RAW.exists():
        print("  fig09 skipped: results/alisim_sweep.csv missing")
        return
    cells = defaultdict(list)  # (algo, model) -> mae
    with open(RAW) as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            try:
                if "/fixed_L500k/" not in r["dir"]: continue
                if abs(float(r["branch_len"]) - 0.05) > 1e-7: continue
                a = r["algo"]; k=int(r["kmer"]); s=int(r["size"]); b=int(r["bits"])
                ssa=float(r["ss_a"]); ssb=float(r["ss_base"])
                bb=int(r["bd_bbits"]); s64=int(r["bd_s64"]); dr=int(r["drlevel"])
                ok=False
                if a=="minhash"   and k==21 and s==1024: ok=True
                if a=="fastkmv"   and k==21 and s==1024: ok=True
                if a=="probmh"    and k==21 and s==1024: ok=True
                if a=="hll"       and k==21 and b==13:   ok=True
                if a=="setsketch" and k==21 and b==13 and ssa==20 and ssb==2.0: ok=True
                if a=="bindash"   and k==21 and bb==16 and s64==16: ok=True
                if a=="kssd"      and dr==3: ok=True
                if not ok: continue
                cells[(a, r["model"])].append(float(r["MAE"]))
            except Exception:
                continue
    models = sorted(set(m for _,m in cells))
    if not models:
        print("  fig09 skipped: no rows found")
        return
    M = np.full((len(ALGOS), len(models)), np.nan)
    for i,a in enumerate(ALGOS):
        for j,m in enumerate(models):
            vs = cells.get((a,m))
            if vs: M[i,j] = np.mean(vs)
    fig, ax = plt.subplots(figsize=(13, 4))
    im = ax.imshow(M, aspect="auto",
                   norm=mcolors.LogNorm(vmin=0.005, vmax=1.0), cmap="YlOrRd")
    ax.set_xticks(range(len(models))); ax.set_xticklabels(models, rotation=30, ha="right")
    ax.set_yticks(range(len(ALGOS))); ax.set_yticklabels(ALGOS)
    for i in range(len(ALGOS)):
        for j in range(len(models)):
            v = M[i,j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.3f}",
                        ha="center", va="center", fontsize=7,
                        color="black" if v<0.05 else "white")
    fig.colorbar(im, ax=ax, label="MAE", shrink=0.85)
    ax.set_title("Fig 9. Per-evolutionary-model robustness (L=500k, BL=0.05, default params)")
    save("fig09_per_model", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 10 — BinDash bbits × s64 heatmap
# ════════════════════════════════════════════════════════════════════════════
def fig10(rows):
    L=500_000; BL=0.05
    sub=[r for r in rows if r["mode"]=="fixed" and r["algo"]=="bindash"
                            and r["L"]==L and r["kmer"]==21
                            and abs(r["BL"]-BL)<1e-7]
    cells=defaultdict(list)
    for r in sub: cells[(r["bd_bbits"], r["bd_s64"])].append(r["mae"])
    bbs = sorted(set(b for b,_ in cells))
    s64s = sorted(set(s for _,s in cells))
    M = np.full((len(bbs), len(s64s)), np.nan)
    for i,b in enumerate(bbs):
        for j,s in enumerate(s64s):
            if (b,s) in cells: M[i,j] = np.mean(cells[(b,s)])
    fig, ax = plt.subplots(figsize=(8, 5))
    im=ax.imshow(M, aspect="auto", origin="lower",
                 norm=mcolors.LogNorm(vmin=0.005, vmax=0.05), cmap="YlOrRd")
    ax.set_xticks(range(len(s64s))); ax.set_xticklabels(s64s)
    ax.set_yticks(range(len(bbs)));  ax.set_yticklabels(bbs)
    ax.set_xlabel("s64 (each = 64 regs)"); ax.set_ylabel("bbits")
    for i in range(len(bbs)):
        for j in range(len(s64s)):
            v = M[i,j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                        fontsize=7, color="black" if v<0.025 else "white")
    fig.colorbar(im, ax=ax, label="MAE")
    ax.set_title(f"Fig 10. BinDash bbits × s64  (L={L:,}, BL={BL}, k=21)")
    save("fig10_bindash_bbits_s64", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 11 — k-mer × BL heatmaps per algorithm
# ════════════════════════════════════════════════════════════════════════════
def fig11(rows):
    L=500_000
    algos=["minhash","fastkmv","probmh","bindash","hll","setsketch"]
    BLs=sorted(set(r["BL"] for r in rows if r["mode"]=="fixed"))
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    axes=axes.flatten()
    for ax, a in zip(axes, algos):
        cand=[r for r in rows
              if r["mode"]=="fixed" and r["algo"]==a and r["L"]==L]
        if a in ("minhash","fastkmv","probmh"):
            cand=[r for r in cand if r["size"]==1024]
        if a=="hll":      cand=[r for r in cand if r["bits"]==13]
        if a=="setsketch":cand=[r for r in cand if r["bits"]==13 and r["ss_a"]==5.0 and r["ss_base"]==2.0]
        if a=="bindash":  cand=[r for r in cand if r["bd_bbits"]==16 and r["bd_s64"]==16]
        ks = sorted(set(r["kmer"] for r in cand))
        if not ks:
            ax.axis("off"); continue
        M = np.full((len(BLs), len(ks)), np.nan)
        for i,bl in enumerate(BLs):
            for j,k in enumerate(ks):
                ms=[r["mae"] for r in cand if r["kmer"]==k and abs(r["BL"]-bl)<1e-7]
                if ms: M[i,j] = np.mean(ms)
        im = ax.imshow(M, aspect="auto", origin="lower",
                       norm=mcolors.LogNorm(vmin=0.0005, vmax=0.5),
                       cmap="YlOrRd")
        ax.set_xticks(range(len(ks))); ax.set_xticklabels(ks)
        ax.set_yticks(range(len(BLs))); ax.set_yticklabels([f"{b:.4f}" for b in BLs], fontsize=6)
        ax.set_title(a); ax.set_xlabel("k-mer"); ax.set_ylabel("BL")
    fig.colorbar(im, ax=axes.tolist(), shrink=0.6, label="MAE (log)")
    fig.suptitle(f"Fig 11. k-mer × BL heatmap (L={L:,}, default non-k params)",
                 fontsize=11, y=0.99)
    save("fig11_kmer_x_BL", fig)

# ════════════════════════════════════════════════════════════════════════════
# Fig 12 — SetSketch base & bits sweeps
# ════════════════════════════════════════════════════════════════════════════
def fig12(rows):
    Ls_pick=[10_000, 100_000, 1_000_000]
    Lcolors=["#888888", "#444444", "#000000"]
    BL=0.05
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # 12-a: ss_base sweep (fixed bits=13 ss_a=5 k=21)
    ax = axes[0]
    for L, c in zip(Ls_pick, Lcolors):
        sub=[r for r in rows
             if r["mode"]=="fixed" and r["algo"]=="setsketch" and r["L"]==L
                and r["kmer"]==21 and r["bits"]==13 and r["ss_a"]==5.0
                and abs(r["BL"]-BL)<1e-7]
        by_b=defaultdict(list)
        for r in sub: by_b[r["ss_base"]].append(r["mae"])
        bb=sorted(by_b.keys())
        if not bb: continue
        ms=[np.mean(by_b[b]) for b in bb]
        ax.plot(bb, ms, "-o", color=c, label=f"L={L_fmt(L)}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("ss_base"); ax.set_ylabel("MAE")
    ax.set_title("SetSketch ss_base sweep (a=5, bits=13, k=21)")
    ax.legend(loc="best", fontsize=8); ax.grid(True, which="both", alpha=0.3)

    # 12-b: bits sweep (fixed ss_a=5 ss_base=2 k=21)
    ax = axes[1]
    for L, c in zip(Ls_pick, Lcolors):
        sub=[r for r in rows
             if r["mode"]=="fixed" and r["algo"]=="setsketch" and r["L"]==L
                and r["kmer"]==21 and r["ss_a"]==5.0 and r["ss_base"]==2.0
                and abs(r["BL"]-BL)<1e-7]
        by_b=defaultdict(list)
        for r in sub: by_b[r["bits"]].append(r["mae"])
        bb=sorted(by_b.keys())
        if not bb: continue
        ms=[np.mean(by_b[b]) for b in bb]
        ax.plot(bb, ms, "-o", color=c, label=f"L={L_fmt(L)}")
    ax.set_yscale("log")
    ax.set_xlabel("bits (log2 regs)"); ax.set_ylabel("MAE")
    ax.set_title("SetSketch bits sweep (a=5, base=2, k=21)")
    ax.legend(loc="best", fontsize=8); ax.grid(True, which="both", alpha=0.3)

    fig.suptitle("Fig 12. SetSketch ss_base & bits sweeps", fontsize=11, y=1.02)
    save("fig12_setsketch_base_bits", fig)

# ════════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════════
def main():
    rows = load()
    print()
    print("Generating figures…")
    fig01(rows)
    fig02(rows)
    fig03(rows)
    fig04(rows)
    fig05(rows)
    fig06(rows)
    fig07(rows)
    fig08(rows)
    fig09()                 # uses raw sweep csv
    fig10(rows)
    fig11(rows)
    fig12(rows)
    print()
    print(f"All figures saved to: {OUT}")

if __name__ == "__main__":
    main()
