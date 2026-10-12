#!/usr/bin/env python3
"""
Variable-length accuracy landscape (default parameters).

Each cell: MAE at nominal L with per-sequence lengths drawn from [L/4, L×4].
Averaged over 3 substitution models: JC69, GTR_G05, GTR_ATrich.

Output: results/figures/accuracy_landscape_varlen_default.png
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from plot_common import (
    LANDSCAPE_ALGO_COLORS,
    LANDSCAPE_FIG_RC,
    LANDSCAPE_MAE_CMAP,
    landscape_mae_norm,
    varlen_landscape_bls,
    varlen_landscape_xtick_labels,
)

HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE / "results" / "alisim_master.csv"
OUT = HERE / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

ALGOS = ["hll", "fastkmv", "probmh", "bindash", "setsketch", "kssd"]
ALGO_LABEL = {
    "hll": "HLL",
    "fastkmv": "MinHash",
    "probmh": "ProbMH",
    "bindash": "BinDash",
    "setsketch": "SetSketch",
    "kssd": "KSSD",
}
ALGO_COLOR = LANDSCAPE_ALGO_COLORS
MAE_CMAP = LANDSCAPE_MAE_CMAP

BINDASH_S64 = 16  # 2 KB default

L_ORDER = ["1k", "10k", "50k", "100k", "500k", "1M", "2M", "4M", "8M"]
L_DISPLAY = {s: s for s in L_ORDER}

mpl.rcParams.update({**LANDSCAPE_FIG_RC, "xtick.major.width": 0.6, "ytick.major.width": 0.6})


def text_color_for_bg(hex_color: str) -> str:
    r, g, b = mcolors.to_rgb(hex_color)
    lum = 0.299 * r + 0.587 * g + 0.114 * b
    return "white" if lum < 0.58 else "#1a1a1a"


def style_heatmap_axes(ax, nrows: int, ncols: int):
    ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nrows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6, alpha=0.85)
    ax.tick_params(which="minor", size=0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color("#333333")


def load_varlen_default(metric: str = "MAE_mean", master_csv: Path | None = None) -> pd.DataFrame:
    df = pd.read_csv(master_csv or DEFAULT_CSV)
    d = df[(df["is_default"] == 1) & (df["mode"] == "varlen")].copy()
    bd = df[
        (df["algo"] == "bindash")
        & (df["mode"] == "varlen")
        & (df["kmer"] == 21)
        & (df["bd_bbits"] == 16)
        & (df["bd_s64"] == BINDASH_S64)
    ].copy()
    d = d[d["algo"] != "bindash"]
    d = pd.concat([d, bd], ignore_index=True)
    d["metric"] = d[metric].astype(float)
    return d[d["algo"].isin(ALGOS)]


def pivot_metric(d: pd.DataFrame, algo: str, bls: list[float]) -> np.ndarray:
    sub = d[d["algo"] == algo]
    M = np.full((len(L_ORDER), len(bls)), np.nan)
    for i, lstr in enumerate(L_ORDER):
        for j, bl in enumerate(bls):
            cell = sub[(sub["Lstr"] == lstr) & (np.isclose(sub["BL"], bl))]
            if len(cell):
                M[i, j] = cell["metric"].mean()
    return M


def build_winner_tables(d: pd.DataFrame, bls: list[float]) -> np.ndarray:
    nL, nB = len(L_ORDER), len(bls)
    win_idx = np.full((nL, nB), -1, dtype=int)
    for i, lstr in enumerate(L_ORDER):
        for j, bl in enumerate(bls):
            cell = d[(d["Lstr"] == lstr) & (np.isclose(d["BL"], bl))]
            if cell.empty:
                continue
            ranked = cell.groupby("algo")["metric"].mean().sort_values()
            win_idx[i, j] = ALGOS.index(ranked.index[0])
    return win_idx


def plot_varlen_landscape(
    metric: str = "MAE_mean",
    slug: str = "accuracy_landscape_varlen_default",
    *,
    master_csv: Path | None = None,
    subtitle: str = "",
):
    d = load_varlen_default(metric, master_csv=master_csv)
    metric_name = "MAE" if metric == "MAE_mean" else "RMSE"
    bls = varlen_landscape_bls()
    xlabels = varlen_landscape_xtick_labels()

    mats = {algo: pivot_metric(d, algo, bls) for algo in ALGOS}
    all_vals = np.concatenate([m[np.isfinite(m)] for m in mats.values()])
    norm = landscape_mae_norm(all_vals)
    win_idx = build_winner_tables(d, bls)

    ncols = len(bls)
    fig = plt.figure(figsize=(max(14, ncols * 1.55), 11.5))
    gs_top = GridSpec(2, 3, figure=fig, top=0.98, bottom=0.52,
                      left=0.05, right=0.90, hspace=0.20, wspace=0.10)
    gs_bot = GridSpec(1, 1, figure=fig, top=0.40, bottom=0.07,
                      left=0.05, right=0.90)

    im_ref = None
    for k, algo in enumerate(ALGOS):
        row, col = divmod(k, 3)
        ax = fig.add_subplot(gs_top[row, col])
        M = mats[algo]
        im_ref = ax.imshow(M, aspect="auto", origin="upper", norm=norm, cmap=MAE_CMAP)
        style_heatmap_axes(ax, M.shape[0], M.shape[1])
        ax.set_title(ALGO_LABEL[algo], pad=4)
        ax.set_xticks(range(len(bls)))
        if row == 1:
            ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=12)
            if col == 1:
                ax.set_xlabel("Mash distance", labelpad=6)
        else:
            ax.set_xticklabels([])
        ax.set_yticks(range(len(L_ORDER)))
        ax.set_yticklabels([L_DISPLAY[l] for l in L_ORDER])
        if col == 0:
            ax.set_ylabel("Nominal length")
        else:
            ax.set_yticklabels([])

    cax = fig.add_axes([0.915, 0.52, 0.012, 0.44])
    cb = fig.colorbar(im_ref, cax=cax)
    cb.set_label(f"{metric_name} (log scale)")
    cb.ax.tick_params(labelsize=12)
    cb.outline.set_linewidth(0.6)

    ax_w = fig.add_subplot(gs_bot[0])
    winner_rgb = np.zeros((len(L_ORDER), len(bls), 3))
    for i in range(len(L_ORDER)):
        for j in range(len(bls)):
            idx = win_idx[i, j]
            if idx >= 0:
                winner_rgb[i, j] = mcolors.to_rgb(ALGO_COLOR[ALGOS[idx]])
            else:
                winner_rgb[i, j] = (0.9, 0.9, 0.9)
    ax_w.imshow(winner_rgb, aspect="auto", origin="upper")
    style_heatmap_axes(ax_w, len(L_ORDER), len(bls))
    ax_w.set_xticks(range(len(bls)))
    ax_w.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=13)
    ax_w.set_yticks(range(len(L_ORDER)))
    ax_w.set_yticklabels([L_DISPLAY[l] for l in L_ORDER], fontsize=13)
    ax_w.set_ylabel("Nominal length", fontsize=13)
    ax_w.set_xlabel("Mash distance", fontsize=13)

    for i in range(len(L_ORDER)):
        for j in range(len(bls)):
            idx = win_idx[i, j]
            if idx < 0:
                continue
            algo = ALGOS[idx]
            fg = text_color_for_bg(ALGO_COLOR[algo])
            ax_w.text(j, i, ALGO_LABEL[algo], ha="center", va="center",
                      fontsize=14, color=fg)

    sub = (
        "Variable length: each sequence in [L/4, L×4]; "
        "3 models (JC69, GTR_G05, GTR_ATrich); SetSketch=MLE a=5"
    )
    fig.suptitle(subtitle or sub, fontsize=12, y=1.01)

    png = OUT / f"{slug}.png"
    fig.savefig(png, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"wrote {png} ({os.path.getsize(png)//1024} KB)")


if __name__ == "__main__":
    plot_varlen_landscape()
