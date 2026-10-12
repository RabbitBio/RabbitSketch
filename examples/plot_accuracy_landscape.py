#!/usr/bin/env python3
"""
Accuracy landscape under default parameters (fixed length).

Outputs (results/figures/):
  accuracy_landscape_default.{png,pdf}
      — all algos at library defaults; BinDash s64=16 (2 KB)
  accuracy_landscape_bindash_8k.{png,pdf}
      — same but BinDash s64=64 (8 KB, equal-memory vs other ~8 KB algos)
  accuracy_landscape_rmse_* — RMSE versions of the above

Rows    : sequence length L
Columns : Mash distance 0.001 … 0.35 (nearest AliSim BL per column)
"""

from __future__ import annotations

import argparse
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
    landscape_mash_bls,
    landscape_mash_xtick_labels,
    setsketch_param_mask,
)

# Four sketch algorithms at ~8 KB (no BinDash / KSSD).
EQ8K_FOUR_ALGOS = ["hll", "setsketch", "fastkmv", "probmh"]
EQ8K_FOUR_LABELS = {
    "hll": "HLL",
    "setsketch": "SetSketch",
    "fastkmv": "FastKMV",
    "probmh": "ProbMinHash",
}

HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE / "results" / "alisim_master.csv"
CSV = DEFAULT_CSV
OUT = HERE / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

ALGOS_ALL = ["hll", "fastkmv", "probmh", "bindash", "setsketch", "kssd"]
ALGOS = ALGOS_ALL  # backward-compatible alias
ALGO_LABEL_BASE = {
    "hll": "HLL",
    "fastkmv": "MinHash",
    "probmh": "ProbMH",
    "bindash": "BinDash",
    "setsketch": "SetSketch",
    "kssd": "KSSD",
}

# BinDash variants (k=21, bbits=16). s64×64 = number of 64-bit bins in sketch.
# Default: s64=16 → 1024 bins, ~2 KB (not the ~8 KB s64=64 equal-memory variant).
BINDASH_VARIANTS = {
    "2k": {"s64": 16, "bins": 1024, "panel_label": "BinDash"},
    "8k": {"s64": 64, "bins": 4096, "panel_label": "BinDash"},
}
ALGO_COLOR = LANDSCAPE_ALGO_COLORS
MAE_CMAP = LANDSCAPE_MAE_CMAP

# Lstr in alisim_master.csv matches build_master_csv.parseL (no leading "L")
L_ORDER = ["1k", "10k", "50k", "100k", "500k", "1M", "2M", "4M", "8M"]
L_DISPLAY = {s: s for s in L_ORDER}

mpl.rcParams.update({**LANDSCAPE_FIG_RC, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.titleweight": "normal"})


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


def load_eq8k_four_sketch(
    metric: str = "MAE_mean",
    *,
    master_csv: Path | None = None,
) -> pd.DataFrame:
    """HLL/SetSketch bits=13 (~8192 regs); FastKMV/ProbMH size=1024 (~8 KB)."""
    df = pd.read_csv(master_csv or CSV)
    fx = df[df["mode"] == "fixed"]
    masks = {
        "hll": (fx["algo"] == "hll") & (fx["kmer"] == 21) & (fx["bits"] == 13),
        "setsketch": setsketch_param_mask(fx, kmer=21, bits=13),
        "fastkmv": (fx["algo"] == "fastkmv") & (fx["kmer"] == 21) & (fx["size"] == 1024),
        "probmh": (fx["algo"] == "probmh") & (fx["kmer"] == 21) & (fx["size"] == 1024),
    }
    parts = []
    for algo in EQ8K_FOUR_ALGOS:
        sub = fx[masks[algo]].copy()
        if sub.empty:
            raise SystemExit(f"no rows for {algo} at ~8 KB params in {master_csv or CSV}")
        sub["metric"] = sub[metric].astype(float)
        parts.append(sub)
    return pd.concat(parts, ignore_index=True)


def load_default(
    metric: str = "MAE_mean",
    *,
    bindash_mem: str = "2k",
    master_csv: Path | None = None,
) -> pd.DataFrame:
    if bindash_mem not in BINDASH_VARIANTS:
        raise ValueError(f"bindash_mem must be one of {list(BINDASH_VARIANTS)}")
    s64 = BINDASH_VARIANTS[bindash_mem]["s64"]

    df = pd.read_csv(master_csv or CSV)
    d = df[(df["is_default"] == 1) & (df["mode"] == "fixed")].copy()
    # BinDash: always pick explicit s64 (default 2k = s64=16, 1024 bins), not sweep is_default alone
    bd = df[
        (df["algo"] == "bindash")
        & (df["mode"] == "fixed")
        & (df["kmer"] == 21)
        & (df["bd_bbits"] == 16)
        & (df["bd_s64"] == s64)
    ].copy()
    d = d[d["algo"] != "bindash"]
    d = pd.concat([d, bd], ignore_index=True)
    d["metric"] = d[metric].astype(float)
    return d[d["algo"].isin(ALGOS)]


def algo_labels(bindash_mem: str) -> dict[str, str]:
    labels = dict(ALGO_LABEL_BASE)
    labels["bindash"] = BINDASH_VARIANTS[bindash_mem]["panel_label"]
    return labels


def algo_panel_positions(n: int) -> list[tuple[int, int]]:
    """Grid (row, col) in the top panel row for n algorithms."""
    if n == 6:
        return [divmod(k, 3) for k in range(6)]
    if n == 5:
        return [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    if n == 4:
        return [(0, 0), (0, 1), (1, 0), (1, 1)]
    raise ValueError(f"unsupported algorithm count: {n}")


def pivot_metric(d: pd.DataFrame, algo: str, bls: list[float]) -> np.ndarray:
    sub = d[d["algo"] == algo]
    M = np.full((len(L_ORDER), len(bls)), np.nan)
    for i, lstr in enumerate(L_ORDER):
        for j, bl in enumerate(bls):
            cell = sub[(sub["Lstr"] == lstr) & (np.isclose(sub["BL"], bl))]
            if len(cell):
                M[i, j] = cell["metric"].mean()
    return M


def build_winner_tables(
    d: pd.DataFrame, bls: list[float], algos: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return winner index, best MAE, margin (2nd - 1st) per (L, BL) cell."""
    nL, nB = len(L_ORDER), len(bls)
    win_idx = np.full((nL, nB), -1, dtype=int)
    best = np.full((nL, nB), np.nan)
    margin = np.full((nL, nB), np.nan)

    for i, lstr in enumerate(L_ORDER):
        for j, bl in enumerate(bls):
            cell = d[(d["Lstr"] == lstr) & (np.isclose(d["BL"], bl))]
            if cell.empty:
                continue
            ranked = cell.groupby("algo")["metric"].mean().sort_values()
            if len(ranked) == 0:
                continue
            best_algo = ranked.index[0]
            win_idx[i, j] = algos.index(best_algo)
            best[i, j] = ranked.iloc[0]
            if len(ranked) >= 2:
                margin[i, j] = ranked.iloc[1] - ranked.iloc[0]
            else:
                margin[i, j] = 0.0
    return win_idx, best, margin


def plot_landscape(
    metric: str = "MAE_mean",
    slug: str = "accuracy_landscape_default",
    *,
    algos: list[str] | None = None,
    bindash_mem: str = "2k",
    subtitle: str = "",
    master_csv: Path | None = None,
    out_dir: Path | None = None,
    png_only: bool = False,
    labels: dict[str, str] | None = None,
    data: pd.DataFrame | None = None,
):
    algos = list(algos or ALGOS_ALL)
    if labels is None:
        labels = algo_labels(bindash_mem)
    if data is None:
        d = load_default(metric, bindash_mem=bindash_mem, master_csv=master_csv)
        d = d[d["algo"].isin(algos)]
    else:
        d = data[data["algo"].isin(algos)].copy()
    metric_name = "MAE" if metric == "MAE_mean" else "RMSE"

    bls = landscape_mash_bls()
    xlabels = landscape_mash_xtick_labels()

    mats = {}
    for algo in algos:
        mats[algo] = pivot_metric(d, algo, bls)
    all_vals = np.concatenate([m[np.isfinite(m)] for m in mats.values()])
    norm = landscape_mae_norm(all_vals)

    win_idx, _, _ = build_winner_tables(d, bls, algos)
    positions = algo_panel_positions(len(algos))
    last_row = max(r for r, _ in positions)

    ncols = len(bls)
    n_algos = len(algos)
    top_cols = 2 if n_algos == 4 else 3
    fig = plt.figure(figsize=(max(14 if n_algos == 4 else 16, ncols * 1.35), 11.5))
    gs_top = GridSpec(
        2, top_cols, figure=fig,
        top=0.98, bottom=0.52,
        left=0.05, right=0.90,
        hspace=0.20, wspace=0.10,
    )
    gs_bot = GridSpec(
        1, 1, figure=fig,
        top=0.40, bottom=0.07,
        left=0.05, right=0.90,
    )

    # ── panels A: per-algorithm MAE heatmaps (2 × 3) ──
    im_ref = None
    for k, algo in enumerate(algos):
        row, col = positions[k]
        ax = fig.add_subplot(gs_top[row, col])
        M = mats[algo]
        im_ref = ax.imshow(
            M, aspect="auto", origin="upper",
            norm=norm, cmap=MAE_CMAP,
        )
        style_heatmap_axes(ax, M.shape[0], M.shape[1])
        ax.set_title(labels[algo], pad=4)
        ax.set_xticks(range(len(bls)))
        if row == last_row:
            ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=12)
            if col == 1:
                ax.set_xlabel("Mash distance", labelpad=6)
        else:
            ax.set_xticklabels([])
        ax.set_yticks(range(len(L_ORDER)))
        ax.set_yticklabels([L_DISPLAY[l] for l in L_ORDER])
        if col == 0:
            ax.set_ylabel("Sequence length")
        else:
            ax.set_yticklabels([])

    # shared colorbar for MAE panels
    cax = fig.add_axes([0.915, 0.52, 0.012, 0.44])
    cb = fig.colorbar(im_ref, cax=cax)
    cb.set_label(f"{metric_name} (log scale)")
    cb.ax.tick_params(labelsize=12)
    cb.outline.set_linewidth(0.6)

    # ── panel B: winner heatmap ──
    ax_w = fig.add_subplot(gs_bot[0])
    # RGB image from winner colors
    winner_rgb = np.zeros((len(L_ORDER), len(bls), 3))
    for i in range(len(L_ORDER)):
        for j in range(len(bls)):
            idx = win_idx[i, j]
            if idx >= 0:
                winner_rgb[i, j] = mcolors.to_rgb(ALGO_COLOR[algos[idx]])
            else:
                winner_rgb[i, j] = (0.9, 0.9, 0.9)
    ax_w.imshow(winner_rgb, aspect="auto", origin="upper")
    style_heatmap_axes(ax_w, len(L_ORDER), len(bls))
    ax_w.set_xticks(range(len(bls)))
    ax_w.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=13)
    ax_w.set_yticks(range(len(L_ORDER)))
    ax_w.set_yticklabels([L_DISPLAY[l] for l in L_ORDER], fontsize=13)
    ax_w.set_ylabel("Sequence length", fontsize=13)
    ax_w.set_xlabel("Mash distance", fontsize=13)

    for i in range(len(L_ORDER)):
        for j in range(len(bls)):
            idx = win_idx[i, j]
            if idx < 0:
                continue
            algo = algos[idx]
            name = labels[algo]
            fg = text_color_for_bg(ALGO_COLOR[algo])
            ax_w.text(j, i, name, ha="center", va="center", fontsize=14,
                      color=fg, fontweight="normal")

    if subtitle:
        fig.suptitle(subtitle, fontsize=12, y=1.01)

    out = out_dir or OUT
    out.mkdir(parents=True, exist_ok=True)
    png = out / f"{slug}.png"
    fig.savefig(png, bbox_inches="tight", dpi=200)
    if not png_only:
        fig.savefig(out / f"{slug}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {png} ({os.path.getsize(png)//1024} KB)")
    if not png_only:
        print(f"wrote {out / f'{slug}.pdf'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="", help="Path to alisim_master.csv")
    ap.add_argument(
        "--only",
        default="default",
        choices=["default", "all"],
        help="default: MAE default landscape only; all: all four figures",
    )
    ap.add_argument(
        "--subtitle",
        default="",
        help="Optional figure suptitle (e.g. model subset)",
    )
    ap.add_argument(
        "--no-kssd",
        action="store_true",
        help="Exclude KSSD from panels and winner",
    )
    ap.add_argument("--slug", default="", help="Output basename (default from flags)")
    ap.add_argument(
        "--eq8k-four",
        action="store_true",
        help="Four sketch algos at ~8 KB (HLL/SetSketch bits=13, FastKMV/ProbMH 1024); no BinDash",
    )
    args = ap.parse_args()
    master_csv = Path(args.master) if args.master else None
    subtitle = args.subtitle
    if args.eq8k_four:
        algos = EQ8K_FOUR_ALGOS
        slug_mae = args.slug or "accuracy_landscape_eq8k_four"
        if not subtitle:
            subtitle = (
                r"$k{=}21$, ~8 KB: HLL/SetSketch $2^{13}$ regs; FastKMV/ProbMinHash $1024$ regs "
                r"(7 substitution models, mean MAE)"
            )
        data = load_eq8k_four_sketch("MAE_mean", master_csv=master_csv)
        plot_landscape(
            "MAE_mean",
            slug_mae,
            algos=algos,
            labels=EQ8K_FOUR_LABELS,
            data=data,
            subtitle=subtitle,
            master_csv=master_csv,
        )
    else:
        algos = [a for a in ALGOS_ALL if not (args.no_kssd and a == "kssd")]
        if args.slug:
            slug_mae = args.slug
        elif args.no_kssd:
            slug_mae = "accuracy_landscape_default_no_kssd"
        else:
            slug_mae = "accuracy_landscape_default"

        plot_landscape(
            "MAE_mean",
            slug_mae,
            algos=algos,
            bindash_mem="2k",
            subtitle=subtitle,
            master_csv=master_csv,
        )
        if args.only == "all":
            plot_landscape(
                "RMSE_mean", "accuracy_landscape_rmse_default",
                bindash_mem="2k", master_csv=master_csv,
            )
            plot_landscape(
                "MAE_mean", "accuracy_landscape_bindash_8k",
                bindash_mem="8k", master_csv=master_csv,
            )
            plot_landscape(
                "RMSE_mean", "accuracy_landscape_rmse_bindash_8k",
                bindash_mem="8k", master_csv=master_csv,
            )
