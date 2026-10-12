#!/usr/bin/env python3
"""
SetSketch parameter sensitivity visualization.

Layout (4-panel figure):
  Panel A  – k-mer size (k=11–31, bits=13, a=5, B=2)
             Lines: one per sequence length, x=k, y=median MAE averaged over BL
  Panel B  – Sketch size / bits (bits=8–16, k=21, a=5, B=2)
             Lines: one per sequence length, x=2^bits, y=median MAE
  Panel C  – Fill parameter a (a=0.3–120, k=21, bits=13, B=2)
             Lines: one per sequence length, x=a (log scale), y=median MAE
  Panel D  – Heatmap: median MAE(L × BL) for default params (k=21,b=13,a=5,B=2)

Also saves individual panels for potential inclusion.
"""

import csv, math, statistics
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

HERE = Path(__file__).parent
CSV  = HERE / "results" / "alisim_setsketch_mle.csv"
OUT  = HERE / "results" / "setsketch_params.pdf"
OUT_PNG = HERE / "results" / "setsketch_params.png"

# ── Load ──────────────────────────────────────────────────────────────────
rows = []
with open(CSV) as f:
    for r in csv.DictReader(f):
        rows.append({
            "k":    int(r["kmer"]),
            "bits": int(r["bits"]),
            "a":    float(r["ss_a"]),
            "B":    float(r["ss_base"]),
            "L":    int(r["n"]) if False else None,  # placeholder
            "MAE":  float(r["MAE"]),
            "BL":   float(r["branch_len"]),
            "dir":  r["dir"],
        })

# Recover L from directory path e.g. .../fixed_L10k/...
def parse_L(d):
    import re
    m = re.search(r'fixed_L(\d+)(k|M)?', d)
    if not m: return None
    v, u = int(m.group(1)), m.group(2)
    if u == "k": return v * 1000
    if u == "M": return v * 1_000_000
    return v

for r in rows:
    r["L"] = parse_L(r["dir"])

rows = [r for r in rows if r["L"] is not None]

LENGTHS = sorted({r["L"] for r in rows})
LBL = {1000:"1k", 10000:"10k", 50000:"50k", 100000:"100k",
       500000:"500k", 1000000:"1M", 2000000:"2M", 4000000:"4M", 8000000:"8M"}

# Color palette: light→dark for short→long
CMAP = plt.cm.viridis_r
L_COLORS = {L: CMAP(i / (len(LENGTHS)-1)) for i, L in enumerate(LENGTHS)}

def med_by(rows, groupfn, xfn):
    """Return dict: x → {L: median MAE}"""
    data = defaultdict(lambda: defaultdict(list))
    for r in rows:
        data[xfn(r)][r["L"]].append(r["MAE"])
    return {x: {L: statistics.median(v) for L, v in Ld.items()}
            for x, Ld in data.items()}

# ─────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(17, 14))
fig.suptitle("SetSketch Parameter Sensitivity  (median MAE averaged over JC69 + GTR_G05, 8 branch lengths)",
             fontsize=12, fontweight="bold", y=0.99)

gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.34,
                      left=0.07, right=0.97, top=0.94, bottom=0.07)

# ─── Panel A: k-mer size ─────────────────────────────────────────────────
ax_k = fig.add_subplot(gs[0, 0])
sel_k = [r for r in rows if r["bits"]==13 and r["a"]==5.0 and r["B"]==2.0]
xvals_k = sorted({r["k"] for r in sel_k})
data_k = med_by(sel_k, None, lambda r: r["k"])

for L in LENGTHS:
    y = [data_k[x].get(L) for x in xvals_k]
    if any(v is not None for v in y):
        ax_k.plot(xvals_k, y, marker="o", markersize=4,
                  color=L_COLORS[L], linewidth=1.6, label=LBL[L])

ax_k.axvline(21, color="dimgray", linestyle="--", linewidth=1.2, alpha=0.8)
ax_k.set_xlabel("k-mer size  k", fontsize=10)
ax_k.set_ylabel("Median MAE", fontsize=10)
ax_k.set_title("(A) k-mer size  (bits=13, a=5, B=2)\n"
               "Long seqs: k barely matters; short/distant: smaller k better",
               fontsize=9, fontweight="bold")
ax_k.set_xticks(xvals_k)
ax_k.set_xticklabels([str(k) for k in xvals_k], fontsize=8)
ax_k.set_yscale("log")
ax_k.grid(True, alpha=0.3)
ax_k.legend(title="Seq length", fontsize=7, title_fontsize=7,
            ncol=2, loc="upper right")
# mark default
ylims = ax_k.get_ylim()
ax_k.annotate("k=21\n(default)", xy=(21, ylims[0]*2), xytext=(23, ylims[0]*3),
              fontsize=7, color="dimgray",
              arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.8))

# ─── Panel B: bits (sketch size) ─────────────────────────────────────────
ax_b = fig.add_subplot(gs[0, 1])
sel_b = [r for r in rows if r["k"]==21 and r["a"]==5.0 and r["B"]==2.0]
xvals_b = sorted({r["bits"] for r in sel_b})
data_b = med_by(sel_b, None, lambda r: r["bits"])

for L in LENGTHS:
    y = [data_b[x].get(L) for x in xvals_b]
    if any(v is not None for v in y):
        lw = 2.0 if L == 1000 else 1.6
        ls = "--" if L == 1000 else "-"
        ax_b.plot([2**x for x in xvals_b], y, marker="s", markersize=4,
                  color=L_COLORS[L], linewidth=lw, linestyle=ls, label=LBL[L])

ax_b.axvline(2**13, color="dimgray", linestyle="--", linewidth=1.2, alpha=0.8)
ax_b.set_xlabel("Register count  m = 2^bits", fontsize=10)
ax_b.set_ylabel("Median MAE", fontsize=10)
ax_b.set_title("(B) Sketch size (bits)  (k=21, a=5, B=2)\n"
               "L≥100k: more bits always better; L=1k: U-shape (fill-rate limited)",
               fontsize=9, fontweight="bold")
ax_b.set_xscale("log", base=2)
ax_b.set_yscale("log")
ax_b.set_xticks([2**b for b in xvals_b])
ax_b.set_xticklabels([f"$2^{{{b}}}$" for b in xvals_b], fontsize=8)
ax_b.grid(True, alpha=0.3)
ax_b.legend(title="Seq length", fontsize=7, title_fontsize=7,
            ncol=2, loc="upper right")
ylims_b = ax_b.get_ylim()
ax_b.annotate("m=8192\n(b=13 default)", xy=(2**13, ylims_b[0]*3),
              xytext=(2**11, ylims_b[0]*5),
              fontsize=7, color="dimgray",
              arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.8))

# ─── Panel C: fill parameter a ───────────────────────────────────────────
ax_a = fig.add_subplot(gs[1, 0])
sel_a = [r for r in rows if r["k"]==21 and r["bits"]==13 and r["B"]==2.0]
xvals_a = sorted({r["a"] for r in sel_a})
data_a = med_by(sel_a, None, lambda r: r["a"])

for L in LENGTHS:
    y = [data_a[x].get(L) for x in xvals_a]
    if any(v is not None for v in y):
        lw = 2.0 if L == 1000 else 1.6
        ax_a.plot(xvals_a, y, marker="^", markersize=4,
                  color=L_COLORS[L], linewidth=lw, label=LBL[L])

ax_a.axvline(5.0, color="dimgray", linestyle="--", linewidth=1.2, alpha=0.8)
ax_a.set_xlabel("Fill parameter  a", fontsize=10)
ax_a.set_ylabel("Median MAE", fontsize=10)
ax_a.set_title("(C) Fill parameter a  (k=21, bits=13, B=2)\n"
               "Large a: saturates long seqs; small a: under-fills short seqs",
               fontsize=9, fontweight="bold")
ax_a.set_xscale("log")
ax_a.set_yscale("log")
ax_a.set_xticks(xvals_a)
ax_a.set_xticklabels([str(a) for a in xvals_a], fontsize=7, rotation=45)
ax_a.grid(True, alpha=0.3)
ax_a.legend(title="Seq length", fontsize=7, title_fontsize=7,
            ncol=2, loc="upper left")
ylims_a = ax_a.get_ylim()
ax_a.annotate("a=5\n(default)", xy=(5, ylims_a[0]*2),
              xytext=(7, ylims_a[0]*5),
              fontsize=7, color="dimgray",
              arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.8))

# ─── Panel D: Heatmap (L × BL) at default params ─────────────────────────
ax_h = fig.add_subplot(gs[1, 1])

sel_def = [r for r in rows if r["k"]==21 and r["bits"]==13
           and r["a"]==5.0 and r["B"]==2.0]
BLS  = sorted({r["BL"] for r in sel_def})
LENS = sorted({r["L"]  for r in sel_def})

heat = np.full((len(BLS), len(LENS)), np.nan)
for i, bl in enumerate(BLS):
    for j, L in enumerate(LENS):
        v = [r["MAE"] for r in sel_def if abs(r["BL"]-bl)<1e-9 and r["L"]==L]
        heat[i, j] = statistics.median(v) if v else np.nan

# Use log scale for color
im = ax_h.imshow(np.log10(heat + 1e-8), aspect="auto", origin="lower",
                 cmap="RdYlGn_r", vmin=-3, vmax=0)
ax_h.set_xticks(range(len(LENS)))
ax_h.set_xticklabels([LBL[L] for L in LENS], fontsize=8, rotation=45)
BL_LABELS = {0.0005:"0.0005\n(≈0.1%)", 0.0025:"0.0025\n(≈0.5%)",
             0.005:"0.005\n(≈1%)",   0.01:"0.01\n(≈2%)",
             0.025:"0.025\n(≈5%)",   0.05:"0.05\n(≈10%)",
             0.1:"0.10\n(≈18%)",     0.15:"0.15\n(≈26%)"}
ax_h.set_yticks(range(len(BLS)))
ax_h.set_yticklabels([BL_LABELS.get(bl, f"{bl:.4f}") for bl in BLS], fontsize=7)
ax_h.set_xlabel("Sequence length", fontsize=10)
ax_h.set_ylabel("Branch length  (≈ Mash distance)", fontsize=10)
ax_h.set_title("(D) MAE heatmap  (default params k=21, bits=13, a=5)\n"
               "Green = accurate; red = inaccurate",
               fontsize=9, fontweight="bold")

# Annotate cells
for i in range(len(BLS)):
    for j in range(len(LENS)):
        val = heat[i, j]
        if not np.isnan(val):
            txt = f"{val:.3f}" if val >= 0.001 else f"{val:.4f}"
            color = "white" if np.log10(val+1e-8) < -1.5 else "black"
            ax_h.text(j, i, txt, ha="center", va="center",
                      fontsize=6, color=color)

cbar = fig.colorbar(im, ax=ax_h, fraction=0.046, pad=0.04)
cbar.set_label("log₁₀(MAE)", fontsize=8)
cbar.set_ticks([-3, -2, -1, 0])
cbar.set_ticklabels(["0.001", "0.01", "0.1", "1.0"], fontsize=7)

# ─── Save ─────────────────────────────────────────────────────────────────
plt.savefig(OUT, dpi=150, bbox_inches="tight")
plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
print(f"Saved: {OUT}")
print(f"Saved: {OUT_PNG}")
