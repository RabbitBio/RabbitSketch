#!/usr/bin/env python3
"""
Sketch-size sensitivity at k=21 (x-axis = number of registers, not memory).

Per-algorithm register count:
  MinHash / ProbMH / FastKMV : size (K)
  HLL / SetSketch            : 2^bits
  BinDash                    : bd_s64 × 64 bins

Columns = fixed Mash levels 0.01, 0.05, 0.10, 0.20, 0.40 (nearest AliSim BL);
          panel titles show pooled exact J at that BL (not Mash).
Rows    = sequence length L

Output: results/figures/size_sensitivity_LxJ_grid.{png,pdf}
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
DEFAULT_K = 21
# Default register counts (library defaults at k=21).
DEFAULT_REGS = {
    "fastkmv": 1024,
    "probmh": 1024,
    "hll": 2 ** 13,
    "setsketch": 2 ** 13,
    "bindash": 16 * 64,  # s64=16 → 1024 bins
}

L_ORDER = ["1k", "100k", "1M", "2M", "4M", "8M"]
L_DISPLAY = {s: s for s in L_ORDER}

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
    "ss_a", "ss_base", "bd_bbits", "bd_s64", "MAE_mean", "exactJ_mean",
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": FS_BASE,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def size_sweep_mask(df: pd.DataFrame, algo: str) -> pd.Series:
    m = (df["algo"] == algo) & (df["kmer"] == DEFAULT_K)
    if algo in ("fastkmv", "probmh", "hll"):
        return m
    if algo == "setsketch":
        return m & setsketch_param_mask(df, kmer=DEFAULT_K)
    if algo == "bindash":
        return m & (df["bd_bbits"] == 16)
    return m


def add_n_regs(df: pd.DataFrame) -> pd.Series:
    """Sketch size as register/bin count (comparable x-axis across algorithms)."""
    nregs = pd.Series(np.nan, index=df.index, dtype=float)
    mh = df["algo"].isin(["fastkmv", "probmh"])
    nregs.loc[mh] = df.loc[mh, "size"].astype(float)
    hb = df["algo"].isin(["hll", "setsketch"])
    nregs.loc[hb] = 2.0 ** df.loc[hb, "bits"]
    bd = df["algo"] == "bindash"
    nregs.loc[bd] = df.loc[bd, "bd_s64"].astype(float) * 64.0
    return nregs


def load_size_sweep(master_csv: Path | None = None) -> pd.DataFrame:
    df = pd.read_csv(master_csv or MASTER, usecols=READ_COLS)
    df = df[(df["mode"] == "fixed") & (df["algo"].isin(ALGOS))]
    parts = [df[size_sweep_mask(df, a)] for a in ALGOS]
    out = pd.concat(parts, ignore_index=True)
    out["n_regs"] = add_n_regs(out)
    return out.dropna(subset=["n_regs"])


def j_labels_for_bls(bls: list[float], master_csv: Path | None = None) -> dict[float, float]:
    df = pd.read_csv(master_csv or MASTER, usecols=["mode", "BL", "exactJ_mean"])
    df = df[df["mode"] == "fixed"]
    j_all = j_mean_by_bl(df)
    return {bl: j_all[bl] for bl in bls if bl in j_all}


def build_lookup(df: pd.DataFrame) -> dict[tuple[str, float, str, float], float]:
    g = df.groupby(["Lstr", "BL", "algo", "n_regs"], as_index=False)["MAE_mean"].mean()
    return {
        (r.Lstr, float(r.BL), r.algo, float(r.n_regs)): float(r.MAE_mean)
        for r in g.itertuples(index=False)
    }


def setsketch_default_lookup(master_csv: Path | None = None) -> dict[tuple[str, float], float]:
    """MLE default (a=5, bits=13, k=21) for open-circle marker."""
    df = pd.read_csv(master_csv or MASTER, usecols=READ_COLS)
    m = (
        (df["mode"] == "fixed")
        & (df["algo"] == "setsketch")
        & setsketch_param_mask(df, kmer=DEFAULT_K, bits=SETSKETCH_BITS_DEFAULT)
    )
    sub = df[m]
    g = sub.groupby(["Lstr", "BL"], as_index=False)["MAE_mean"].mean()
    return {(r.Lstr, float(r.BL)): float(r.MAE_mean) for r in g.itertuples(index=False)}


def reg_ticks() -> tuple[list[float], list[str]]:
    ticks = [256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    labels = ["256", "512", "1k", "2k", "4k", "8k", "16k", "32k", "64k"]
    return ticks, labels


def plot_grid(
    lookup: dict,
    bls: list[float],
    j_map: dict[float, float],
    ss_default: dict[tuple[str, float], float],
    png_only: bool,
) -> None:
    nL, nB = len(L_ORDER), len(bls)
    fig, axes = plt.subplots(
        nL, nB, figsize=(nB * 3.6, nL * 2.9), sharex=True, sharey=True, squeeze=False,
    )
    xticks, xlabels = reg_ticks()
    y_floor, y_ceil = 5e-4, 1.0

    for i, lstr in enumerate(L_ORDER):
        for j, bl in enumerate(bls):
            ax = axes[i, j]
            for algo in ALGOS:
                regs = sorted(
                    {k[3] for k in lookup if k[0] == lstr and k[1] == bl and k[2] == algo}
                )
                pts = [(m, lookup[(lstr, bl, algo, m)]) for m in regs if (lstr, bl, algo, m) in lookup]
                if len(pts) < 2:
                    continue
                xs = [p[0] for p in pts]
                ys = [max(p[1], y_floor) for p in pts]
                ax.plot(xs, ys, "-o", ms=5.5, lw=1.8, color=ALGO_COLOR[algo], label=ALGO_LABEL[algo])
                if algo == "setsketch":
                    key = (lstr, bl)
                    if key in ss_default:
                        ax.scatter(
                            [DEFAULT_REGS["setsketch"]], [max(ss_default[key], y_floor)],
                            s=65, facecolors="none", edgecolors=ALGO_COLOR[algo],
                            linewidths=2.0, zorder=5,
                        )
                else:
                    dregs = DEFAULT_REGS[algo]
                    idx = min(range(len(xs)), key=lambda t: abs(xs[t] - dregs))
                    ax.scatter([xs[idx]], [ys[idx]], s=65, facecolors="none",
                               edgecolors=ALGO_COLOR[algo], linewidths=2.0, zorder=5)

            ax.axvline(8192, color="#cccccc", ls="--", lw=0.8, zorder=0)
            ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            ax.set_xlim(200, 70000)
            ax.set_ylim(y_floor, y_ceil)
            ax.set_xticks(xticks)
            ax.grid(True, which="both", color="#eeeeee", lw=0.5)
            jv = j_map.get(bl, float("nan"))
            ax.set_title(
                f"$L$={L_DISPLAY[lstr]}, $J$={format_j_tick(jv)}"
                if np.isfinite(jv) else f"$L$={L_DISPLAY[lstr]}",
                pad=5, fontsize=FS_PANEL_TITLE,
            )
            if i == nL - 1:
                ax.set_xticklabels(xlabels, fontsize=FS_TICK_X_BOTTOM)
                ax.set_xlabel("Sketch size", fontsize=FS_XLABEL, labelpad=6)
            if j == 0:
                ax.set_ylabel("MAE", fontsize=FS_YLABEL)

    handles = [
        plt.Line2D([0], [0], color=ALGO_COLOR[a], lw=2.2, marker="o", ms=7, label=ALGO_LABEL[a])
        for a in ALGOS
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=FS_LEGEND, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(
        "Sketch-size sensitivity at $k=21$ (x-axis = register count; open circles = default size)",
        fontsize=FS_SUP, y=0.995,
    )
    fig.subplots_adjust(left=0.07, right=0.99, top=0.94, bottom=0.10, hspace=0.34, wspace=0.20)

    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / "size_sensitivity_LxJ_grid"
    fig.savefig(f"{stem}.png", dpi=300)
    print(f"wrote {stem}.png", flush=True)
    if not png_only:
        fig.savefig(f"{stem}.pdf")
        print(f"wrote {stem}.pdf", flush=True)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="", help="Path to alisim_master.csv")
    ap.add_argument("--png-only", action="store_true")
    args = ap.parse_args()

    master_csv = Path(args.master) if args.master else None
    df = load_size_sweep(master_csv)
    bls = sensitivity_mash_bls()
    lookup = build_lookup(df)
    j_map = j_labels_for_bls(bls, master_csv)
    cols = ", ".join(
        f"M={m:.2f}→BL={bl} ($J$={format_j_tick(j_map.get(bl, float('nan')))})"
        for m, bl in SENSITIVITY_MASH_COLUMNS
    )
    print(f"loaded {len(df)} size-sweep rows; columns: {cols}", flush=True)
    ss_def = setsketch_default_lookup(master_csv)
    plot_grid(lookup, bls, j_map, ss_def, args.png_only)


if __name__ == "__main__":
    main()
