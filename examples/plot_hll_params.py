#!/usr/bin/env python3
"""
HLL parameter sensitivity visualization (3-panel + heatmap).

Layout:
  Panel A  – k-mer size (k=11–31, bits=13)
             Lines: one per seq length; x=k, y=median MAE averaged over all 11 models
  Panel B  – Sketch size / bits (bits=8–16, k=21)
             Lines: one per seq length, log-log; shows 29% MAE-per-doubling law
  Panel C  – Register-doubling MAE reduction ratio (bits=8–16, k=21)
             Bar chart per length: actual reduction % vs theoretical 29.3%
  Panel D  – Heatmap (L × BL) at default params (k=21, bits=13)
"""

import csv, re, math, statistics
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

HERE = Path(__file__).parent
CSV  = HERE / "results" / "alisim_sweep.csv"
OUT_PNG = HERE / "results" / "hll_params.png"
OUT_PDF = HERE / "results" / "hll_params.pdf"

# ── Load HLL rows ─────────────────────────────────────────────────────────
def parse_L(d):
    m = re.search(r'fixed_L(\d+)(k|M)?', d)
    if not m: return None
    v, u = int(m.group(1)), m.group(2)
    return v*1000 if u=="k" else v*1_000_000 if u=="M" else v

rows = []
with open(CSV) as f:
    for r in csv.DictReader(f):
        if r["algo"] != "hll": continue
        L = parse_L(r["dir"])
        if L is None: continue
        rows.append({
            "k":    int(r["kmer"]),
            "bits": int(r["bits"]),
            "L":    L,
            "BL":   float(r["branch_len"]),
            "MAE":  float(r["MAE"]),
            "model": r["model"],
        })

LENGTHS = sorted({r["L"] for r in rows})
BLS     = sorted({r["BL"] for r in rows})
MODELS  = sorted({r["model"] for r in rows})
LBL = {1000:"1k", 10000:"10k", 50000:"50k", 100000:"100k",
       500000:"500k", 1000000:"1M", 2000000:"2M", 4000000:"4M", 8000000:"8M"}

CMAP = plt.cm.viridis_r
L_COLORS = {L: CMAP(i/(len(LENGTHS)-1)) for i, L in enumerate(LENGTHS)}

def med_group(rows, xfn):
    data = defaultdict(lambda: defaultdict(list))
    for r in rows:
        data[xfn(r)][r["L"]].append(r["MAE"])
    return {x: {L: statistics.median(v) for L, v in Ld.items()}
            for x, Ld in data.items()}

# ─────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(17, 14))
fig.suptitle("HLL Parameter Sensitivity  (median MAE, averaged over 11 evolutionary models)",
             fontsize=12, fontweight="bold", y=0.99)
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.34,
                      left=0.07, right=0.97, top=0.94, bottom=0.07)

# ─── Panel A: k-mer size ─────────────────────────────────────────────────
ax_k = fig.add_subplot(gs[0, 0])
sel_k = [r for r in rows if r["bits"]==13]
xvals_k = sorted({r["k"] for r in sel_k})
data_k = med_group(sel_k, lambda r: r["k"])

for L in LENGTHS:
    y = [data_k[x].get(L) for x in xvals_k]
    if any(v is not None for v in y):
        ax_k.plot(xvals_k, y, marker="o", markersize=4,
                  color=L_COLORS[L], linewidth=1.7, label=LBL[L])

ax_k.axvline(21, color="dimgray", linestyle="--", linewidth=1.2, alpha=0.8)
ax_k.set_xlabel("k-mer size  k", fontsize=10)
ax_k.set_ylabel("Median MAE", fontsize=10)
ax_k.set_title("(A) k-mer size  (bits=13)\n"
               "Optimal k shifts with length; HLL less k-sensitive than MinHash",
               fontsize=9, fontweight="bold")
ax_k.set_xticks(xvals_k)
ax_k.set_xticklabels([str(k) for k in xvals_k], fontsize=8)
ax_k.set_yscale("log")
ax_k.grid(True, alpha=0.3)
ax_k.legend(title="Seq length", fontsize=7, title_fontsize=7,
            ncol=2, loc="upper right")
ylims = ax_k.get_ylim()
ax_k.text(21.3, (ylims[0]*ylims[1])**0.5, "← k=21\n  (default)",
          fontsize=7, color="dimgray", va="center")

# ─── Panel B: bits — log-log, showing 1/√2 slope ─────────────────────────
ax_b = fig.add_subplot(gs[0, 1])
sel_b = [r for r in rows if r["k"]==21]
xvals_b = sorted({r["bits"] for r in sel_b})
data_b = med_group(sel_b, lambda r: r["bits"])

for L in LENGTHS:
    y = [data_b[x].get(L) for x in xvals_b]
    if any(v is not None for v in y):
        ax_b.plot([2**x for x in xvals_b], y, marker="s", markersize=4,
                  color=L_COLORS[L], linewidth=1.7, label=LBL[L])

# Draw theoretical slope: 1/sqrt(m) = 1/2^(bits/2)  → slope -0.5 on log2-log2
m_ref = np.array([2**b for b in xvals_b], dtype=float)
# Anchor to the median across lengths at bits=13
anchor_mae = statistics.median([data_b[13][L] for L in LENGTHS if L in data_b[13]])
theory = anchor_mae * np.sqrt(2**13 / m_ref)
ax_b.plot(m_ref, theory, "k--", linewidth=1.2, alpha=0.6, label="Theory (∝1/√m)")

ax_b.axvline(2**13, color="dimgray", linestyle=":", linewidth=1.2, alpha=0.7)
ax_b.set_xlabel("Register count  m = 2^bits", fontsize=10)
ax_b.set_ylabel("Median MAE", fontsize=10)
ax_b.set_title("(B) Sketch size (bits)  (k=21)\n"
               "MAE ∝ 1/√m: every register doubling reduces MAE by ~29%",
               fontsize=9, fontweight="bold")
ax_b.set_xscale("log", base=2)
ax_b.set_yscale("log")
ax_b.set_xticks([2**b for b in xvals_b])
ax_b.set_xticklabels([f"$2^{{{b}}}$" for b in xvals_b], fontsize=8)
ax_b.grid(True, alpha=0.3)
handles, labels = ax_b.get_legend_handles_labels()
ax_b.legend(handles, labels, title="Seq length", fontsize=7, title_fontsize=7,
            ncol=2, loc="upper right")
ylims_b = ax_b.get_ylim()
ax_b.text(2**13 * 1.15, (ylims_b[0]*ylims_b[1])**0.5,
          "← m=8192\n  (b=13, default)",
          fontsize=7, color="dimgray", va="center")

# ─── Panel C: register-doubling reduction ratio per length ────────────────
ax_r = fig.add_subplot(gs[1, 0])

# For each L, compute actual % MAE reduction per bit doubling (bits b → b+1)
xvals_bits = sorted({r["bits"] for r in rows if r["k"]==21})
reduction_data = {}  # L → list of (bits_step_label, reduction%)
for L in LENGTHS:
    reductions = []
    prev_mae = None
    for b in xvals_bits:
        v = [r["MAE"] for r in rows if r["k"]==21 and r["bits"]==b and r["L"]==L]
        mae = statistics.median(v) if v else None
        if mae is not None and prev_mae is not None and prev_mae > 0:
            pct = (prev_mae - mae) / prev_mae * 100
            reductions.append(pct)
        prev_mae = mae
    reduction_data[L] = reductions

# x-axis: bit transitions (8→9, 9→10, ..., 15→16)
transitions = [f"{b}→{b+1}" for b in xvals_bits[:-1]]
x = np.arange(len(transitions))
n_L = len(LENGTHS)
width = 0.08
offsets = np.linspace(-(n_L-1)/2, (n_L-1)/2, n_L) * width

for i, L in enumerate(LENGTHS):
    reds = reduction_data[L]
    ax_r.bar(x + offsets[i], reds, width=width*0.9,
             color=L_COLORS[L], alpha=0.85, label=LBL[L])

# Theoretical line
ax_r.axhline(29.3, color="red", linestyle="--", linewidth=1.5,
             label="Theory 29.3%  (1−1/√2)")
ax_r.set_xlabel("Register-count doubling step  (bits b → b+1)", fontsize=10)
ax_r.set_ylabel("MAE reduction  (%)", fontsize=10)
ax_r.set_title("(C) MAE reduction per register doubling  (k=21)\n"
               "Bars cluster around 29% — confirms HLL scaling law",
               fontsize=9, fontweight="bold")
ax_r.set_xticks(x)
ax_r.set_xticklabels(transitions, fontsize=8)
ax_r.set_ylim(0, 55)
ax_r.axhline(0, color="black", linewidth=0.5)
ax_r.grid(True, axis="y", alpha=0.3)
ax_r.legend(title="Seq length", fontsize=6, title_fontsize=7,
            ncol=3, loc="upper right")

# ─── Panel D: Heatmap (L × BL) at default params ─────────────────────────
ax_h = fig.add_subplot(gs[1, 1])
sel_def = [r for r in rows if r["k"]==21 and r["bits"]==13]

heat = np.full((len(BLS), len(LENGTHS)), np.nan)
for i, bl in enumerate(BLS):
    for j, L in enumerate(LENGTHS):
        v = [r["MAE"] for r in sel_def if abs(r["BL"]-bl)<1e-9 and r["L"]==L]
        heat[i, j] = statistics.median(v) if v else np.nan

im = ax_h.imshow(np.log10(heat + 1e-8), aspect="auto", origin="lower",
                 cmap="RdYlGn_r", vmin=-3, vmax=0)
ax_h.set_xticks(range(len(LENGTHS)))
ax_h.set_xticklabels([LBL[L] for L in LENGTHS], fontsize=8, rotation=45)
BL_LABELS = {0.0005:"0.0005\n(≈0.1%)", 0.0025:"0.0025\n(≈0.5%)",
             0.005:"0.005\n(≈1%)",   0.01:"0.01\n(≈2%)",
             0.025:"0.025\n(≈5%)",   0.05:"0.05\n(≈10%)",
             0.1:"0.10\n(≈18%)",     0.15:"0.15\n(≈26%)"}
ax_h.set_yticks(range(len(BLS)))
ax_h.set_yticklabels([BL_LABELS.get(bl, f"{bl}") for bl in BLS], fontsize=7)
ax_h.set_xlabel("Sequence length", fontsize=10)
ax_h.set_ylabel("Branch length  (≈ Mash distance)", fontsize=10)
ax_h.set_title("(D) MAE heatmap  (default params k=21, bits=13)\n"
               "Green = accurate; red = inaccurate",
               fontsize=9, fontweight="bold")
for i in range(len(BLS)):
    for j in range(len(LENGTHS)):
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

# ── Save ──────────────────────────────────────────────────────────────────
plt.savefig(OUT_PDF, dpi=150, bbox_inches="tight")
plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
print(f"Saved: {OUT_PNG}")
print(f"Saved: {OUT_PDF}")
