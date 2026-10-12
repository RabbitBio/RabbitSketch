#!/usr/bin/env python3
"""
Parameter sensitivity plots for MinHash, FastKMV, ProbMinHash, KSSD, BinDash.
Each algorithm gets its own 4-panel figure tailored to its parameter space.
"""

import csv, re, math, statistics
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
CSV  = HERE / "results" / "alisim_sweep.csv"

LENGTHS = [1000,10000,50000,100000,500000,1000000,2000000,4000000,8000000]
LBL = {1000:"1k",10000:"10k",50000:"50k",100000:"100k",500000:"500k",
       1000000:"1M",2000000:"2M",4000000:"4M",8000000:"8M"}
BL_LABELS = {0.0005:"0.0005\n(≈0.1%)",0.0025:"0.0025\n(≈0.5%)",
             0.005:"0.005\n(≈1%)",0.01:"0.01\n(≈2%)",
             0.025:"0.025\n(≈5%)",0.05:"0.05\n(≈10%)",
             0.1:"0.10\n(≈18%)",0.15:"0.15\n(≈26%)"}

CMAP = plt.cm.viridis_r
L_COLORS = {L: CMAP(i/(len(LENGTHS)-1)) for i, L in enumerate(LENGTHS)}

def parse_L(d):
    m = re.search(r'fixed_L(\d+)(k|M)?', d)
    if not m: return None
    v,u = int(m.group(1)), m.group(2)
    return v*1000 if u=="k" else v*1_000_000 if u=="M" else v

def load_algo(algo):
    rows = []
    with open(CSV) as f:
        for r in csv.DictReader(f):
            if r["algo"] != algo: continue
            L = parse_L(r["dir"])
            if L is None: continue
            rows.append({
                "k":    int(r["kmer"]) if r["kmer"]!="-1" else -1,
                "size": int(r["size"]) if r["size"]!="-1" else -1,
                "bits": int(r["bits"]) if r["bits"]!="-1" else -1,
                "dr":   int(r["drlevel"]) if r["drlevel"]!="-1" else -1,
                "bbits":int(r["bd_bbits"]) if r["bd_bbits"]!="-1" else -1,
                "s64":  int(r["bd_s64"])   if r["bd_s64"]!="-1"   else -1,
                "L":    L,
                "BL":   float(r["branch_len"]),
                "MAE":  float(r["MAE"]),
                "model":r["model"],
            })
    return rows

def med_by_x(rows, xfn):
    d = defaultdict(lambda: defaultdict(list))
    for r in rows:
        d[xfn(r)][r["L"]].append(r["MAE"])
    return {x: {L: statistics.median(v) for L,v in Ld.items()} for x,Ld in d.items()}

def draw_k_panel(ax, rows_sel, title_note=""):
    xvals = sorted({r["k"] for r in rows_sel})
    data = med_by_x(rows_sel, lambda r: r["k"])
    for L in LENGTHS:
        y = [data[x].get(L) for x in xvals]
        if any(v is not None for v in y):
            ax.plot(xvals, y, marker="o", markersize=4,
                    color=L_COLORS[L], linewidth=1.6, label=LBL[L])
    ax.axvline(21, color="dimgray", linestyle="--", linewidth=1, alpha=0.7)
    ylims = ax.get_ylim()
    ax.text(21.4, (ylims[0]*ylims[1])**0.5 if ylims[0]>0 else ylims[1]*0.5,
            "← k=21\n  (default)", fontsize=7, color="dimgray", va="center")
    ax.set_xlabel("k-mer size  k", fontsize=10)
    ax.set_ylabel("Median MAE", fontsize=10)
    ax.set_title(f"(A) k-mer size{title_note}", fontsize=9, fontweight="bold")
    ax.set_xticks(xvals); ax.set_xticklabels([str(k) for k in xvals], fontsize=8)
    ax.set_yscale("log"); ax.grid(True, alpha=0.3)
    ax.legend(title="Seq length", fontsize=7, title_fontsize=7, ncol=2, loc="upper right")

def draw_size_panel(ax, rows_sel, def_size, title="(B) Sketch size",
                    slope_label="Theory (∝1/√size)"):
    xvals = sorted({r["size"] for r in rows_sel})
    data  = med_by_x(rows_sel, lambda r: r["size"])
    for L in LENGTHS:
        y = [data[x].get(L) for x in xvals]
        if any(v is not None for v in y):
            ax.plot(xvals, y, marker="s", markersize=4,
                    color=L_COLORS[L], linewidth=1.6, label=LBL[L])
    # theoretical slope
    anchor = statistics.median([data[def_size][L] for L in LENGTHS
                                 if L in data.get(def_size,{})])
    xarr = np.array(xvals, dtype=float)
    ax.plot(xarr, anchor * np.sqrt(def_size/xarr), "k--",
            linewidth=1.2, alpha=0.6, label=slope_label)
    ax.axvline(def_size, color="dimgray", linestyle=":", linewidth=1, alpha=0.7)
    ylims = ax.get_ylim()
    ax.text(def_size*1.1, (ylims[0]*ylims[1])**0.5 if ylims[0]>0 else ylims[1]*0.5,
            f"← s={def_size}\n  (default)", fontsize=7, color="dimgray", va="center")
    ax.set_xlabel("Sketch size  s", fontsize=10)
    ax.set_ylabel("Median MAE", fontsize=10)
    ax.set_title(title, fontsize=9, fontweight="bold")
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks(xvals); ax.set_xticklabels([str(s) for s in xvals], fontsize=7, rotation=40)
    ax.grid(True, alpha=0.3)
    ax.legend(title="Seq length", fontsize=7, title_fontsize=7, ncol=2, loc="upper right")

def draw_reduction_panel(ax, rows_sel, xfn, xvals,
                         xlabel="Size doubling step", title="(C) MAE reduction"):
    """Bar chart: % MAE reduction between consecutive xvals."""
    transitions = [f"{xvals[i]}→{xvals[i+1]}" for i in range(len(xvals)-1)]
    x = np.arange(len(transitions))
    n_L = len(LENGTHS)
    width = 0.08
    offsets = np.linspace(-(n_L-1)/2, (n_L-1)/2, n_L) * width
    data = med_by_x(rows_sel, xfn)
    for i, L in enumerate(LENGTHS):
        reds = []
        for j in range(len(xvals)-1):
            m0 = data.get(xvals[j],{}).get(L)
            m1 = data.get(xvals[j+1],{}).get(L)
            if m0 and m1 and m0>0:
                reds.append((m0-m1)/m0*100)
            else:
                reds.append(0)
        ax.bar(x + offsets[i], reds, width=width*0.9,
               color=L_COLORS[L], alpha=0.85, label=LBL[L])
    ax.axhline(29.3, color="red", linestyle="--", linewidth=1.4,
               label="Theory 29.3%")
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("MAE reduction (%)", fontsize=10)
    ax.set_title(title, fontsize=9, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(transitions, fontsize=8, rotation=30)
    ax.set_ylim(0, 60); ax.axhline(0, color="black", linewidth=0.5)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(title="Seq length", fontsize=6, title_fontsize=7,
              ncol=3, loc="upper right")

def draw_heatmap(ax, rows_sel, title="(D) MAE heatmap (default params)"):
    BLS = sorted({r["BL"] for r in rows_sel})
    heat = np.full((len(BLS), len(LENGTHS)), np.nan)
    for i, bl in enumerate(BLS):
        for j, L in enumerate(LENGTHS):
            v = [r["MAE"] for r in rows_sel if abs(r["BL"]-bl)<1e-9 and r["L"]==L]
            heat[i,j] = statistics.median(v) if v else np.nan
    im = ax.imshow(np.log10(heat+1e-8), aspect="auto", origin="lower",
                   cmap="RdYlGn_r", vmin=-3, vmax=0)
    ax.set_xticks(range(len(LENGTHS)))
    ax.set_xticklabels([LBL[L] for L in LENGTHS], fontsize=8, rotation=45)
    ax.set_yticks(range(len(BLS)))
    ax.set_yticklabels([BL_LABELS.get(bl,f"{bl}") for bl in BLS], fontsize=7)
    ax.set_xlabel("Sequence length", fontsize=10)
    ax.set_ylabel("Branch length (≈ Mash distance)", fontsize=10)
    ax.set_title(title, fontsize=9, fontweight="bold")
    for i in range(len(BLS)):
        for j in range(len(LENGTHS)):
            val = heat[i,j]
            if not np.isnan(val):
                txt = f"{val:.3f}" if val>=0.001 else f"{val:.4f}"
                color = "white" if np.log10(val+1e-8)<-1.5 else "black"
                ax.text(j, i, txt, ha="center", va="center", fontsize=6, color=color)
    cbar = plt.gcf().colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("log₁₀(MAE)", fontsize=8)
    cbar.set_ticks([-3,-2,-1,0])
    cbar.set_ticklabels(["0.001","0.01","0.1","1.0"], fontsize=7)

def save_fig(fig, name):
    fig.savefig(HERE/"results"/f"{name}.pdf", dpi=150, bbox_inches="tight")
    fig.savefig(HERE/"results"/f"{name}.png", dpi=150, bbox_inches="tight")
    print(f"Saved: results/{name}.png")
    plt.close(fig)

# ══════════════════════════════════════════════════════════════════════════
# 1. MinHash
# ══════════════════════════════════════════════════════════════════════════
print("Plotting MinHash...")
mh = load_algo("minhash")
SIZES_MH = sorted({r["size"] for r in mh})

fig = plt.figure(figsize=(17,14))
fig.suptitle("MinHash Parameter Sensitivity  (median MAE, 11 models)",
             fontsize=12, fontweight="bold", y=0.99)
gs = fig.add_gridspec(2,2, hspace=0.46, wspace=0.34, left=0.07, right=0.97, top=0.94, bottom=0.07)

ax = fig.add_subplot(gs[0,0])
draw_k_panel(ax, [r for r in mh if r["size"]==1024], "\n(size=1024)")
ax.set_title("(A) k-mer size  (size=1024)\nSmaller k helps near-identical pairs; larger k improves specificity", fontsize=9, fontweight="bold")

ax = fig.add_subplot(gs[0,1])
draw_size_panel(ax, [r for r in mh if r["k"]==21], def_size=1024,
    title="(B) Sketch size  (k=21)\nMAE ∝ 1/√s: doubling size reduces MAE ~29%")

ax = fig.add_subplot(gs[1,0])
draw_reduction_panel(ax, [r for r in mh if r["k"]==21],
    lambda r: r["size"], SIZES_MH,
    xlabel="Sketch-size doubling step  (s → 2s)",
    title="(C) MAE reduction per sketch doubling  (k=21)\nMinHash follows 1/√s law like HLL")

ax = fig.add_subplot(gs[1,1])
draw_heatmap(ax, [r for r in mh if r["k"]==21 and r["size"]==1024],
    "(D) MAE heatmap  (default k=21, size=1024)")

save_fig(fig, "minhash_params")

# ══════════════════════════════════════════════════════════════════════════
# 2. FastKMV
# ══════════════════════════════════════════════════════════════════════════
print("Plotting FastKMV...")
fk = load_algo("fastkmv")
SIZES_FK = sorted({r["size"] for r in fk})

fig = plt.figure(figsize=(17,14))
fig.suptitle("FastKMV Parameter Sensitivity  (median MAE, 11 models)",
             fontsize=12, fontweight="bold", y=0.99)
gs = fig.add_gridspec(2,2, hspace=0.46, wspace=0.34, left=0.07, right=0.97, top=0.94, bottom=0.07)

ax = fig.add_subplot(gs[0,0])
draw_k_panel(ax, [r for r in fk if r["size"]==1024], "\n(size=1024)")
ax.set_title("(A) k-mer size  (size=1024)\nKMV (exact min-hash) vs MinHash: nearly identical k-sensitivity", fontsize=9, fontweight="bold")

ax = fig.add_subplot(gs[0,1])
draw_size_panel(ax, [r for r in fk if r["k"]==21], def_size=1024,
    title="(B) Sketch size  (k=21)\nKMV stores exact minima → same 1/√s law as MinHash")

ax = fig.add_subplot(gs[1,0])
draw_reduction_panel(ax, [r for r in fk if r["k"]==21],
    lambda r: r["size"], SIZES_FK,
    xlabel="Sketch-size doubling step  (s → 2s)",
    title="(C) MAE reduction per sketch doubling  (k=21)")

ax = fig.add_subplot(gs[1,1])
draw_heatmap(ax, [r for r in fk if r["k"]==21 and r["size"]==1024],
    "(D) MAE heatmap  (default k=21, size=1024)")

save_fig(fig, "fastkmv_params")

# ══════════════════════════════════════════════════════════════════════════
# 3. ProbMinHash
# ══════════════════════════════════════════════════════════════════════════
print("Plotting ProbMinHash...")
pmh = load_algo("probmh")
SIZES_PMH = sorted({r["size"] for r in pmh})  # 64,96,128,...4096 (√2 steps)

fig = plt.figure(figsize=(17,14))
fig.suptitle("ProbMinHash Parameter Sensitivity  (median MAE, 11 models)",
             fontsize=12, fontweight="bold", y=0.99)
gs = fig.add_gridspec(2,2, hspace=0.46, wspace=0.34, left=0.07, right=0.97, top=0.94, bottom=0.07)

ax = fig.add_subplot(gs[0,0])
draw_k_panel(ax, [r for r in pmh if r["size"]==1024], "\n(size=1024)")
ax.set_title("(A) k-mer size  (size=1024)\nWeight-proportional sampling amplifies short-k effects", fontsize=9, fontweight="bold")

ax = fig.add_subplot(gs[0,1])
draw_size_panel(ax, [r for r in pmh if r["k"]==21], def_size=1024,
    title="(B) Sketch size  (k=21)\nDenser √2 grid confirms MAE ∝ 1/√s for weighted scheme")

# Panel C: MAE vs size, log-log + slope check per L group
ax = fig.add_subplot(gs[1,0])
sel_pmh_b = [r for r in pmh if r["k"]==21]
data_pmh_s = med_by_x(sel_pmh_b, lambda r: r["size"])
for L in LENGTHS:
    y = [data_pmh_s[s].get(L) for s in SIZES_PMH]
    if any(v is not None for v in y):
        ax.plot(SIZES_PMH, y, marker="^", markersize=4,
                color=L_COLORS[L], linewidth=1.6, label=LBL[L])
# theoretical -0.5 slope from size=1024
anchor = statistics.median([data_pmh_s[1024][L] for L in LENGTHS if L in data_pmh_s.get(1024,{})])
sarr = np.array(SIZES_PMH, dtype=float)
ax.plot(sarr, anchor*np.sqrt(1024/sarr), "k--", lw=1.2, alpha=0.6, label="∝1/√s")
ax.set_xlabel("Sketch size  s", fontsize=10)
ax.set_ylabel("Median MAE", fontsize=10)
ax.set_title("(C) MAE vs size, log-log  (k=21)\n√2 grid resolves slope precisely: all lengths parallel at −0.5",
             fontsize=9, fontweight="bold")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xticks(SIZES_PMH)
ax.set_xticklabels([str(s) for s in SIZES_PMH], fontsize=7, rotation=40)
ax.grid(True, alpha=0.3)
ax.legend(title="Seq length", fontsize=7, title_fontsize=7, ncol=2, loc="upper right")

ax = fig.add_subplot(gs[1,1])
draw_heatmap(ax, [r for r in pmh if r["k"]==21 and r["size"]==1024],
    "(D) MAE heatmap  (default k=21, size=1024)")

save_fig(fig, "probmh_params")

# ══════════════════════════════════════════════════════════════════════════
# 4. KSSD
# ══════════════════════════════════════════════════════════════════════════
print("Plotting KSSD...")
ks = load_algo("kssd")
DR_COLORS = {1:"#1a9641", 2:"#fdae61", 3:"#d73027"}
BLS_KS = sorted({r["BL"] for r in ks})

fig = plt.figure(figsize=(17,14))
fig.suptitle("KSSD Parameter Sensitivity  (median MAE, 11 models)\n"
             "Only parameter: drlevel (subsampling ratio ≈ 1/4^drlevel)",
             fontsize=12, fontweight="bold", y=0.99)
gs = fig.add_gridspec(2,2, hspace=0.46, wspace=0.34, left=0.07, right=0.97, top=0.93, bottom=0.07)

# Panel A: MAE vs seq length per drlevel
ax = fig.add_subplot(gs[0,0])
for dr in [1,2,3]:
    sel = [r for r in ks if r["dr"]==dr]
    by_L = defaultdict(list)
    for r in sel: by_L[r["L"]].append(r["MAE"])
    x = sorted(by_L.keys())
    y = [statistics.median(by_L[L]) for L in x]
    lw = 2.2 if dr==3 else 1.8
    ax.plot([LBL[L] for L in x], y, marker="o", markersize=5,
            color=DR_COLORS[dr], linewidth=lw,
            label=f"drlevel={dr}  (keep≈1/{4**dr} k-mers)")
ax.axhline(0.01, color="gray", linestyle=":", linewidth=1, alpha=0.7)
ax.text(8.3, 0.011, "MAE=0.01", fontsize=7, color="gray")
ax.set_xlabel("Sequence length", fontsize=10)
ax.set_ylabel("Median MAE (log scale)", fontsize=10)
ax.set_title("(A) MAE vs sequence length by drlevel\n"
             "drlevel=1: accurate from L=4M;  drlevel=3: needs L>10M",
             fontsize=9, fontweight="bold")
ax.set_yscale("log"); ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper right")

# Panel B: MAE vs BL (all L averaged) per drlevel
ax = fig.add_subplot(gs[0,1])
for dr in [1,2,3]:
    sel = [r for r in ks if r["dr"]==dr]
    by_BL = defaultdict(list)
    for r in sel: by_BL[r["BL"]].append(r["MAE"])
    x = sorted(by_BL.keys())
    y = [statistics.median(by_BL[bl]) for bl in x]
    ax.plot(x, y, marker="s", markersize=5,
            color=DR_COLORS[dr], linewidth=1.8,
            label=f"drlevel={dr}")
ax.set_xlabel("Branch length (≈ Mash distance)", fontsize=10)
ax.set_ylabel("Median MAE (log scale)", fontsize=10)
ax.set_title("(B) MAE vs divergence by drlevel\n"
             "All drlevels converge at BL≥0.05 (high divergence, sketches empty)",
             fontsize=9, fontweight="bold")
ax.set_yscale("log"); ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper right")

# Bottom row: 3 heatmaps side by side (one per drlevel)
from matplotlib.gridspec import GridSpecFromSubplotSpec
gs_bot = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1,:], wspace=0.08)
for pi, dr in enumerate([1,2,3]):
    ax_h = fig.add_subplot(gs_bot[pi])
    sel = [r for r in ks if r["dr"]==dr]
    heat = np.full((len(BLS_KS), len(LENGTHS)), np.nan)
    for i, bl in enumerate(BLS_KS):
        for j, L in enumerate(LENGTHS):
            v = [r["MAE"] for r in sel if abs(r["BL"]-bl)<1e-9 and r["L"]==L]
            heat[i,j] = statistics.median(v) if v else np.nan
    im = ax_h.imshow(np.log10(heat+1e-8), aspect="auto", origin="lower",
                     cmap="RdYlGn_r", vmin=-3, vmax=0)
    ax_h.set_xticks(range(len(LENGTHS)))
    ax_h.set_xticklabels([LBL[L] for L in LENGTHS], fontsize=6, rotation=50)
    ax_h.set_xlabel("Seq length", fontsize=9)
    ax_h.set_title(f"(C{pi+1}) drlevel={dr}\n(keep ≈1/{4**dr} k-mers)",
                   fontsize=9, fontweight="bold", color=DR_COLORS[dr])
    if pi == 0:
        ax_h.set_yticks(range(len(BLS_KS)))
        ax_h.set_yticklabels([BL_LABELS.get(bl, f"{bl}") for bl in BLS_KS], fontsize=6)
        ax_h.set_ylabel("Branch length", fontsize=9)
    else:
        ax_h.set_yticks([])
    for i in range(len(BLS_KS)):
        for j in range(len(LENGTHS)):
            val = heat[i,j]
            if not np.isnan(val):
                txt = f"{val:.2f}" if val >= 0.01 else f"{val:.3f}"
                clr = "white" if np.log10(val+1e-8) < -1.5 else "black"
                ax_h.text(j, i, txt, ha="center", va="center", fontsize=5, color=clr)
    if pi == 2:
        cbar = fig.colorbar(im, ax=ax_h, fraction=0.046, pad=0.04)
        cbar.set_label("log₁₀(MAE)", fontsize=7)
        cbar.set_ticks([-3,-2,-1,0])
        cbar.set_ticklabels(["0.001","0.01","0.1","1.0"], fontsize=6)

save_fig(fig, "kssd_params")

# ══════════════════════════════════════════════════════════════════════════
# 5. BinDash
# ══════════════════════════════════════════════════════════════════════════
print("Plotting BinDash...")
bd = load_algo("bindash")

fig = plt.figure(figsize=(17,14))
fig.suptitle("BinDash Parameter Sensitivity  (median MAE, 11 models)",
             fontsize=12, fontweight="bold", y=0.99)
gs = fig.add_gridspec(2,2, hspace=0.46, wspace=0.34, left=0.07, right=0.97, top=0.94, bottom=0.07)

# Panel A: k-mer size
ax = fig.add_subplot(gs[0,0])
draw_k_panel(ax, [r for r in bd if r["bbits"]==16 and r["s64"]==16], "\n(bbits=16, s64=16)")
ax.set_title("(A) k-mer size  (bbits=16, s64=16)\nBinDash most sensitive to k among all algorithms", fontsize=9, fontweight="bold")

# Panel B: bbits (bits per bin hash)
ax = fig.add_subplot(gs[0,1])
sel_bb = [r for r in bd if r["k"]==21 and r["s64"]==16]
BBITS = sorted({r["bbits"] for r in sel_bb})
data_bb = med_by_x(sel_bb, lambda r: r["bbits"])
for L in LENGTHS:
    y = [data_bb[x].get(L) for x in BBITS]
    if any(v is not None for v in y):
        ax.plot(BBITS, y, marker="D", markersize=4,
                color=L_COLORS[L], linewidth=1.6, label=LBL[L])
ax.axvline(16, color="dimgray", linestyle="--", linewidth=1, alpha=0.7)
ylims = ax.get_ylim()
ax.text(16.3, (ylims[0]*ylims[1])**0.5 if ylims[0]>0 else ylims[1]*0.5,
        "← bbits=16\n  (default)", fontsize=7, color="dimgray", va="center")
ax.set_xlabel("Bits per bin hash  bbits", fontsize=10)
ax.set_ylabel("Median MAE", fontsize=10)
ax.set_title("(B) Hash precision (bbits)  (k=21, s64=16)\nPlateau beyond bbits=10–12; diminishing returns at >16",
             fontsize=9, fontweight="bold")
ax.set_xticks(BBITS); ax.set_xticklabels([str(b) for b in BBITS], fontsize=8)
ax.set_yscale("log"); ax.grid(True, alpha=0.3)
ax.legend(title="Seq length", fontsize=7, title_fontsize=7, ncol=2, loc="upper right")

# Panel C: s64 (number of 64-element bins → sketch size = s64 × 64)
ax = fig.add_subplot(gs[1,0])
sel_s64 = [r for r in bd if r["k"]==21 and r["bbits"]==16]
S64_VALS = sorted({r["s64"] for r in sel_s64})
data_s64 = med_by_x(sel_s64, lambda r: r["s64"])
for L in LENGTHS:
    y = [data_s64[x].get(L) for x in S64_VALS]
    if any(v is not None for v in y):
        ax.plot([v*64 for v in S64_VALS], y, marker="^", markersize=4,
                color=L_COLORS[L], linewidth=1.6, label=LBL[L])
# theoretical 1/sqrt(s) slope
anchor_s = statistics.median([data_s64[16][L] for L in LENGTHS if L in data_s64.get(16,{})])
sarr = np.array([v*64 for v in S64_VALS], dtype=float)
ax.plot(sarr, anchor_s*np.sqrt(16*64/sarr), "k--", lw=1.2, alpha=0.6, label="∝1/√(s64×64)")
ax.axvline(16*64, color="dimgray", linestyle=":", linewidth=1, alpha=0.7)
ylims = ax.get_ylim()
ax.text(16*64*1.1, (ylims[0]*ylims[1])**0.5 if ylims[0]>0 else ylims[1]*0.5,
        "← s64=16\n  (default)", fontsize=7, color="dimgray", va="center")
ax.set_xlabel("Total sketch size  s64 × 64", fontsize=10)
ax.set_ylabel("Median MAE", fontsize=10)
ax.set_title("(C) Sketch density (s64)  (k=21, bbits=16)\nTotal elements = s64×64; MAE ∝ 1/√(s64×64)",
             fontsize=9, fontweight="bold")
ax.set_xscale("log"); ax.set_yscale("log")
ax.grid(True, alpha=0.3)
ax.legend(title="Seq length", fontsize=7, title_fontsize=7, ncol=2, loc="upper right")

ax = fig.add_subplot(gs[1,1])
draw_heatmap(ax, [r for r in bd if r["k"]==21 and r["bbits"]==16 and r["s64"]==16],
    "(D) MAE heatmap  (default k=21, bbits=16, s64=16)")

save_fig(fig, "bindash_params")

print("All done.")
