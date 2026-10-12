#!/usr/bin/env python3
"""
k-mer sensitivity grid at default sketch size.

Columns = fixed Mash levels 0.01, 0.05, 0.10, 0.20, 0.40 (nearest AliSim BL);
          panel titles show pooled exact J at that BL (not Mash).
Rows    = sequence length L

Output: results/figures/k_sensitivity_LxJ_grid.{png,pdf}
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

_MPL_CACHE = os.environ.get("MPLCONFIGDIR", "/tmp/mpl-cache-rabbitsketch")
os.environ["MPLCONFIGDIR"] = _MPL_CACHE
Path(_MPL_CACHE).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_common import (
    SENSITIVITY_MASH_COLUMNS,
    SETSKETCH_BITS_DEFAULT,
    SETSKETCH_K_DEFAULT,
    SETSKETCH_PANEL_LABEL,
    format_j_tick,
    j_mean_by_bl,
    sensitivity_mash_bls,
    setsketch_param_mask,
)

HERE = Path(__file__).resolve().parent
MASTER = HERE / "results/alisim_master.csv"
OUT = HERE / "results/figures"

ALGOS = ["hll", "fastkmv", "probmh", "bindash", "setsketch"]
ALGO_LABEL = {
    "hll": "HLL",
    "fastkmv": "MinHash",
    "probmh": "ProbMH",
    "bindash": "BinDash",
    "setsketch": SETSKETCH_PANEL_LABEL,
}
ALGO_COLOR = {
    "hll": "#7570B3",
    "setsketch": "#4A4578",
    "fastkmv": "#E69F00",
    "probmh": "#D07A2D",
    "bindash": "#A8562E",
}

L_ORDER = ["1k", "100k", "1M", "2M", "4M", "8M"]
L_DISPLAY = {s: s for s in L_ORDER}
DEFAULT_K = 21

FS_BASE = 13
FS_PANEL_TITLE = 13
FS_AXIS = 14
FS_XLABEL = 16
FS_YLABEL = 14
FS_TICK = 12
FS_TICK_X_BOTTOM = 13
FS_LEGEND = 18
FS_SUP = 16

READ_COLS = [
    "mode", "Lstr", "BL", "algo", "kmer", "size", "bits",
    "ss_a", "ss_base", "bd_bbits", "bd_s64", "drlevel", "MAE_mean", "exactJ_mean",
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": FS_BASE,
    "axes.labelsize": FS_AXIS,
    "axes.titlesize": FS_PANEL_TITLE,
    "xtick.labelsize": FS_TICK,
    "ytick.labelsize": FS_TICK,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def k_sweep_mask(df: pd.DataFrame, algo: str) -> pd.Series:
    m = df["algo"] == algo
    if algo in ("fastkmv", "probmh"):
        return m & (df["size"] == 1024)
    if algo == "hll":
        return m & (df["bits"] == 13)
    if algo == "setsketch":
        return m & setsketch_param_mask(df, bits=SETSKETCH_BITS_DEFAULT)
    if algo == "bindash":
        return m & (df["bd_bbits"] == 16) & (df["bd_s64"] == 16)
    return m


def load_k_sweep(master_csv: Path | None = None) -> pd.DataFrame:
    df = pd.read_csv(master_csv or MASTER, usecols=READ_COLS)
    df = df[(df["mode"] == "fixed") & (df["algo"].isin(ALGOS))]
    parts = [df[k_sweep_mask(df, a)] for a in ALGOS]
    return pd.concat(parts, ignore_index=True)


def j_labels_for_bls(bls: list[float], master_csv: Path | None = None) -> dict[float, float]:
    """Pooled exact J per BL (fixed mode, all algos) — independent of param sweep."""
    df = pd.read_csv(master_csv or MASTER, usecols=["mode", "BL", "exactJ_mean"])
    df = df[df["mode"] == "fixed"]
    j_all = j_mean_by_bl(df)
    return {bl: j_all[bl] for bl in bls if bl in j_all}


def _bl_key(bl: float) -> float:
    """Stable dict key for AliSim BL (avoid float mismatch in lookup)."""
    return round(float(bl), 6)


def build_lookup(df: pd.DataFrame) -> dict[tuple[str, float, str, int], float]:
    g = df.groupby(["Lstr", "BL", "algo", "kmer"], as_index=False)["MAE_mean"].mean()
    return {
        (r.Lstr, _bl_key(r.BL), r.algo, int(r.kmer)): float(r.MAE_mean)
        for r in g.itertuples(index=False)
    }


def setsketch_default_lookup(master_csv: Path | None = None) -> dict[tuple[str, float], float]:
    df = pd.read_csv(master_csv or MASTER, usecols=READ_COLS)
    m = (
        (df["mode"] == "fixed")
        & (df["algo"] == "setsketch")
        & setsketch_param_mask(df, kmer=SETSKETCH_K_DEFAULT, bits=SETSKETCH_BITS_DEFAULT)
    )
    g = df[m].groupby(["Lstr", "BL"], as_index=False)["MAE_mean"].mean()
    return {(r.Lstr, _bl_key(r.BL)): float(r.MAE_mean) for r in g.itertuples(index=False)}


def plot_grid(
    lookup: dict,
    bls: list[float],
    col_title_fn,
    ss_default: dict[tuple[str, float], float],
    png_only: bool,
    stem: str,
) -> None:
    nL, nB = len(L_ORDER), len(bls)
    fig, axes = plt.subplots(
        nL, nB, figsize=(nB * 3.6, nL * 2.9), sharex=True, sharey=True, squeeze=False,
    )

    k_vals = list(range(11, 32, 2))
    y_floor, y_ceil = 5e-4, 1.0

    for i, lstr in enumerate(L_ORDER):
        for j, bl in enumerate(bls):
            ax = axes[i, j]
            bl_key = _bl_key(bl)
            for algo in ALGOS:
                xs, ys = [], []
                for k in k_vals:
                    key = (lstr, bl_key, algo, k)
                    if key not in lookup:
                        continue
                    xs.append(k)
                    ys.append(max(lookup[key], y_floor))
                if len(xs) == 1:
                    ax.scatter(xs, ys, color=ALGO_COLOR[algo], s=28, zorder=3,
                               label=ALGO_LABEL[algo])
                elif len(xs) >= 2:
                    ax.plot(xs, ys, "-o", ms=5.5, lw=1.8,
                            color=ALGO_COLOR[algo], label=ALGO_LABEL[algo])
                    # Open circle only when k=21 is the sole point (no curve yet)
                    if algo == "setsketch" and SETSKETCH_K_DEFAULT not in xs:
                        dk = (lstr, bl_key)
                        if dk in ss_default:
                            ax.scatter(
                                [SETSKETCH_K_DEFAULT], [max(ss_default[dk], y_floor)],
                                s=65, facecolors="none", edgecolors=ALGO_COLOR[algo],
                                linewidths=2.0, zorder=5,
                            )

            ax.axvline(DEFAULT_K, color="#aaaaaa", ls="--", lw=0.8, zorder=0)
            ax.set_yscale("log")
            ax.set_ylim(y_floor, y_ceil)
            ax.grid(True, which="both", color="#eeeeee", lw=0.5)
            ax.set_title(col_title_fn(lstr, bl), pad=5, fontsize=FS_PANEL_TITLE)
            if i == nL - 1:
                ax.set_xlabel("k", fontsize=FS_XLABEL, labelpad=6)
                ax.tick_params(axis="x", labelsize=FS_TICK_X_BOTTOM)
            if j == 0:
                ax.set_ylabel("MAE", fontsize=FS_YLABEL)

    handles = [
        plt.Line2D([0], [0], color=ALGO_COLOR[a], lw=2.2, marker="o", ms=7, label=ALGO_LABEL[a])
        for a in ALGOS
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=FS_LEGEND, bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(left=0.07, right=0.99, top=0.94, bottom=0.10, hspace=0.34, wspace=0.20)

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / stem
    fig.savefig(f"{out}.png", dpi=300)
    print(f"wrote {out}.png", flush=True)
    if not png_only:
        fig.savefig(f"{out}.pdf")
        print(f"wrote {out}.pdf", flush=True)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="", help="Path to alisim_master.csv")
    ap.add_argument("--png-only", action="store_true")
    args = ap.parse_args()

    master_csv = Path(args.master) if args.master else None
    df = load_k_sweep(master_csv)
    bls = sensitivity_mash_bls()
    lookup = build_lookup(df)
    j_map = j_labels_for_bls(bls, master_csv)
    cols = ", ".join(
        f"M={m:.2f}→BL={bl} ($J$={format_j_tick(j_map.get(bl, float('nan')))})"
        for m, bl in SENSITIVITY_MASH_COLUMNS
    )
    ss_def = setsketch_default_lookup(master_csv)
    ss_rows = df[(df["algo"] == "setsketch")]
    k_ss = sorted(ss_rows["kmer"].unique())
    print(
        f"loaded {len(df)} rows; SetSketch a=5 k={k_ss}; columns: {cols}",
        flush=True,
    )

    def title_j(lstr: str, bl: float) -> str:
        jv = j_map.get(bl, float("nan"))
        if np.isfinite(jv):
            return f"$L$={L_DISPLAY[lstr]}, $J$={format_j_tick(jv)}"
        return f"$L$={L_DISPLAY[lstr]}"

    def title_d(lstr: str, bl: float) -> str:
        return f"$L$={L_DISPLAY[lstr]}, $d$={bl:g}"

    plot_grid(lookup, bls, title_j, ss_def, args.png_only, "k_sensitivity_LxJ_grid")
    plot_grid(lookup, bls, title_d, ss_def, args.png_only, "k_sensitivity_LxD_grid")


if __name__ == "__main__":
    main()
