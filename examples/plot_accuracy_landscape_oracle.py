#!/usr/bin/env python3
"""
Oracle accuracy landscape: best MAE per (L, BL) for each algorithm,
chosen over the full parameter sweep in alisim_master.csv.

SetSketch: joint-MLE only (setsketch_mle in sweep → setsketch in master).

Outputs (results/figures/):
  accuracy_landscape_oracle.{png,pdf}
  accuracy_landscape_rmse_oracle.{png,pdf}

Top panels: best MAE heatmaps only (no cell text).
Winner panel: algorithm name + winning config (compact).
X-axis: Mash distance (k=21).
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
)

HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE / "results/alisim_master.csv"
OUT = HERE / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

ALGOS_ALL = ["hll", "fastkmv", "probmh", "bindash", "setsketch", "kssd"]
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

L_ORDER = ["1k", "10k", "50k", "100k", "500k", "1M", "2M", "4M", "8M"]
L_DISPLAY = {s: s for s in L_ORDER}

PARAM_COLS = [
    "kmer", "size", "bits", "ss_a", "ss_base", "bd_bbits", "bd_s64", "drlevel",
]

mpl.rcParams.update(LANDSCAPE_FIG_RC)


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


def load_fixed(metric: str, algos: list[str], master_csv: Path | None = None) -> pd.DataFrame:
    df = pd.read_csv(master_csv or DEFAULT_CSV)
    d = df[(df["mode"] == "fixed") & (df["algo"].isin(algos))].copy()
    d["metric"] = d[metric].astype(float)
    return d


def algo_panel_positions(n: int) -> list[tuple[int, int]]:
    """Grid (row, col) in a 2×3 top panel for n algorithms."""
    if n == 6:
        return [divmod(k, 3) for k in range(6)]
    if n == 5:
        return [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    raise ValueError(f"unsupported algorithm count: {n}")


def key_params_label(algo: str, row: pd.Series) -> str:
    """Compact winning-parameter string for one oracle cell."""
    if algo in ("fastkmv", "probmh"):
        parts = []
        if int(row["kmer"]) != 21:
            parts.append(f"k={int(row['kmer'])}")
        parts.append(f"N={int(row['size'])}")
        return ", ".join(parts)
    if algo == "hll":
        parts = []
        if int(row["kmer"]) != 21:
            parts.append(f"k={int(row['kmer'])}")
        parts.append(f"2^{int(row['bits'])}")
        return ", ".join(parts)
    if algo == "setsketch":
        parts = []
        if int(row["kmer"]) != 21:
            parts.append(f"k={int(row['kmer'])}")
        parts.append(f"a={float(row['ss_a']):g}")
        b = float(row["ss_base"])
        if abs(b - 2.0) > 1e-9:
            parts.append(f"B={b:g}")
        return ", ".join(parts)
    if algo == "bindash":
        parts = []
        if int(row["kmer"]) != 21:
            parts.append(f"k={int(row['kmer'])}")
        if int(row["bd_bbits"]) != 16:
            parts.append(f"b={int(row['bd_bbits'])}")
        parts.append(f"s64={int(row['bd_s64'])}")
        return ", ".join(parts)
    if algo == "kssd":
        return f"dr={int(row['drlevel'])}"
    return ""


def pick_oracle_rows(d: pd.DataFrame) -> pd.DataFrame:
    """One row per (Lstr, BL, algo): minimum metric, tie-break by param cols."""
    work = d.sort_values(["Lstr", "BL", "algo", "metric"] + PARAM_COLS)
    idx = work.groupby(["Lstr", "BL", "algo"], sort=False)["metric"].idxmin()
    best = work.loc[idx].copy()
    best["param_note"] = [
        key_params_label(a, best.loc[i]) for i, a in best["algo"].items()
    ]
    return best


def pivot_oracle(best: pd.DataFrame, algo: str, bls: list[float]
                 ) -> tuple[np.ndarray, np.ndarray]:
    sub = best[best["algo"] == algo]
    M = np.full((len(L_ORDER), len(bls)), np.nan)
    P = np.empty((len(L_ORDER), len(bls)), dtype=object)
    P[:] = ""
    for i, lstr in enumerate(L_ORDER):
        for j, bl in enumerate(bls):
            cell = sub[(sub["Lstr"] == lstr) & (np.isclose(sub["BL"], bl))]
            if len(cell):
                r = cell.iloc[0]
                M[i, j] = float(r["metric"])
                P[i, j] = str(r["param_note"])
    return M, P


def build_winner_oracle(
    best: pd.DataFrame, bls: list[float], algos: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    nL, nB = len(L_ORDER), len(bls)
    win_idx = np.full((nL, nB), -1, dtype=int)
    win_param = np.empty((nL, nB), dtype=object)
    win_param[:] = ""

    for i, lstr in enumerate(L_ORDER):
        for j, bl in enumerate(bls):
            cell = best[(best["Lstr"] == lstr) & (np.isclose(best["BL"], bl))]
            if cell.empty:
                continue
            cell = cell.sort_values(["metric"] + PARAM_COLS)
            r = cell.iloc[0]
            win_idx[i, j] = algos.index(r["algo"])
            win_param[i, j] = str(r["param_note"])
    return win_idx, win_param


def plot_oracle(
    metric: str = "MAE_mean",
    slug: str = "accuracy_landscape_oracle",
    *,
    algos: list[str] | None = None,
    master_csv: Path | None = None,
    subtitle: str = "",
    out_dir: Path | None = None,
    png_only: bool = False,
    save_configs: bool = True,
):
    algos = list(algos or ALGOS_ALL)
    d = load_fixed(metric, algos, master_csv=master_csv)
    metric_name = "MAE" if metric == "MAE_mean" else "RMSE"
    best = pick_oracle_rows(d)
    bls = landscape_mash_bls()
    xlabels = landscape_mash_xtick_labels()

    mats, params = {}, {}
    for algo in algos:
        mats[algo], params[algo] = pivot_oracle(best, algo, bls)

    all_vals = np.concatenate([m[np.isfinite(m)] for m in mats.values()])
    norm = landscape_mae_norm(all_vals)

    win_idx, win_param = build_winner_oracle(best, bls, algos)
    positions = algo_panel_positions(len(algos))
    last_row = max(r for r, _ in positions)

    fig = plt.figure(figsize=(16, 11.5))
    gs_top = GridSpec(2, 3, figure=fig, top=0.98, bottom=0.52,
                      left=0.05, right=0.90, hspace=0.20, wspace=0.10)
    gs_bot = GridSpec(1, 1, figure=fig, top=0.40, bottom=0.07,
                      left=0.05, right=0.90)

    im_ref = None
    for k, algo in enumerate(algos):
        row, col = positions[k]
        ax = fig.add_subplot(gs_top[row, col])
        M = mats[algo]
        im_ref = ax.imshow(M, aspect="auto", origin="upper", norm=norm, cmap=MAE_CMAP)
        style_heatmap_axes(ax, M.shape[0], M.shape[1])
        ax.set_title(f"{ALGO_LABEL[algo]} (best)", pad=4)
        ax.set_xticks(range(len(bls)))
        if row == last_row:
            ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=11)
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

    cax = fig.add_axes([0.915, 0.52, 0.012, 0.44])
    cb = fig.colorbar(im_ref, cax=cax)
    cb.set_label(f"{metric_name} (log scale)")
    cb.ax.tick_params(labelsize=12)

    ax_w = fig.add_subplot(gs_bot[0])
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
    ax_w.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=12)
    ax_w.set_yticks(range(len(L_ORDER)))
    ax_w.set_yticklabels([L_DISPLAY[l] for l in L_ORDER], fontsize=13)
    ax_w.set_ylabel("Sequence length", fontsize=13)
    ax_w.set_xlabel("Mash distance", fontsize=13)
    ax_w.set_title("Oracle winner (best params per algorithm)", pad=8)

    for i in range(len(L_ORDER)):
        for j in range(len(bls)):
            idx = win_idx[i, j]
            if idx < 0:
                continue
            algo = algos[idx]
            fg = text_color_for_bg(ALGO_COLOR[algo])
            note = win_param[i, j]
            ax_w.text(j, i - 0.15, ALGO_LABEL[algo], ha="center", va="center",
                      fontsize=11, color=fg)
            if note:
                ax_w.text(j, i + 0.22, note, ha="center", va="center",
                          fontsize=7, color=fg)

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

    if save_configs:
        csv_out = out / f"{slug}_configs.csv"
        cols = ["Lstr", "BL", "algo", "metric", "param_note"] + PARAM_COLS
        best[cols].rename(columns={"metric": metric_name}).to_csv(csv_out, index=False)
        print(f"wrote {csv_out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="", help="Path to alisim_master.csv")
    ap.add_argument(
        "--no-kssd",
        action="store_true",
        help="Exclude KSSD from oracle panels and winner",
    )
    ap.add_argument(
        "--only",
        default="mae",
        choices=["mae", "all"],
        help="mae: MAE oracle only; all: MAE + RMSE",
    )
    ap.add_argument("--slug", default="", help="Output basename (default from flags)")
    args = ap.parse_args()

    master_csv = Path(args.master) if args.master else None
    algos = [a for a in ALGOS_ALL if not (args.no_kssd and a == "kssd")]
    if args.slug:
        slug_mae = args.slug
    elif args.no_kssd:
        slug_mae = "accuracy_landscape_oracle_no_kssd"
    else:
        slug_mae = "accuracy_landscape_oracle"

    plot_oracle(
        "MAE_mean",
        slug_mae,
        algos=algos,
        master_csv=master_csv,
    )
    if args.only == "all":
        slug_rmse = slug_mae.replace("oracle", "rmse_oracle", 1)
        if slug_rmse == slug_mae:
            slug_rmse = f"{slug_mae}_rmse"
        plot_oracle(
            "RMSE_mean",
            slug_rmse,
            algos=algos,
            master_csv=master_csv,
        )
