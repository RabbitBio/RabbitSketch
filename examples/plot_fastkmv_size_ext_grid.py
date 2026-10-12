#!/usr/bin/env python3
"""
FastKMV-only sketch-size sensitivity (JC69, k=21).

Merges:
  - Pre-8k sizes (64 … 4096) from alisim_sweep.JC69.pairs1000.selected6.csv
  - 8k … 8M extension from alisim_sweep.JC69.pairs1000.fastkmv_size_ext.csv
At size=8192 the extension rerun takes precedence.

Output: results/figures/fastkmv_size_full_LxJ_uptoL2M.{png,pdf}
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

_MPL_CACHE = os.environ.get("MPLCONFIGDIR", "/tmp/mpl-cache-rabbitsketch")
os.environ["MPLCONFIGDIR"] = _MPL_CACHE
Path(_MPL_CACHE).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_common import SENSITIVITY_MASH_COLUMNS, format_j_tick

HERE = Path(__file__).resolve().parent
EXT_CSV = HERE / "results/alisim_sweep.JC69.pairs1000.fastkmv_size_ext.csv"
BASE_CSV = HERE / "results/alisim_sweep.JC69.pairs1000.selected6.csv"
OUT = HERE / "results/figures"

L_ORDER = ["L1k", "L10k", "L50k", "L100k", "L500k", "L1M", "L2M"]
L_DISPLAY = {s: s[1:] if s.startswith("L") else s for s in L_ORDER}
COLOR = "#E69F00"
LABEL = "MinHash (FastKMV)"
DEFAULT_LIB = 1024
DEFAULT_EQ = 8192

READ_COLS = [
    "dir", "model", "branch_len", "algo", "kmer", "size", "bits",
    "ss_a", "ss_base", "bd_bbits", "bd_s64", "drlevel",
    "n", "mean_exact_J", "mean_est_J", "MAE", "RMSE",
]

FS_BASE = 13
FS_PANEL_TITLE = 13
FS_XLABEL = 16
FS_YLABEL = 14
FS_TICK_X_BOTTOM = 11
FS_SUP = 15

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": FS_BASE,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def parse_lstr(path: str) -> str | None:
    m = re.search(r"fixed_(L[^/]+)/", str(path))
    return m.group(1) if m else None


def scenario_key(df: pd.DataFrame) -> pd.Series:
    return df["dir"].astype(str) + "|" + df["branch_len"].astype(str)


def load_base_pre8k(base_csv: Path, ext_keys: set[str]) -> pd.DataFrame:
    df = pd.read_csv(base_csv, usecols=READ_COLS, on_bad_lines="skip", low_memory=False)
    for col in ("branch_len", "kmer", "size", "MAE", "mean_exact_J"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["dir", "branch_len", "kmer", "size", "MAE", "algo"])
    df = df[
        (df["algo"] == "fastkmv")
        & (df["kmer"] == 21)
        & (df["size"] < DEFAULT_EQ)
        & scenario_key(df).isin(ext_keys)
    ].copy()
    return df


def load_merged(base_csv: Path, ext_csv: Path, max_l: str = "L2M") -> pd.DataFrame:
    ext = pd.read_csv(ext_csv)
    ext_keys = set(scenario_key(ext))
    pre = load_base_pre8k(base_csv, ext_keys)
    merged = pd.concat([pre, ext], ignore_index=True)
    merged = merged.drop_duplicates(subset=["dir", "branch_len", "size"], keep="last")
    merged["Lstr"] = merged["dir"].map(parse_lstr)
    merged = merged[merged["Lstr"].notna()]
    if max_l in L_ORDER:
        merged = merged[merged["Lstr"].isin(L_ORDER[: L_ORDER.index(max_l) + 1])]
    return merged.sort_values(["Lstr", "branch_len", "size"])


def nearest_bl(target: float, available: list[float]) -> float:
    return min(available, key=lambda x: abs(x - target))


def column_bls(df: pd.DataFrame) -> list[tuple[float, float]]:
    avail = sorted(df["branch_len"].unique())
    return [(target, nearest_bl(target, avail)) for _, target in SENSITIVITY_MASH_COLUMNS]


def j_map_for_bls(df: pd.DataFrame, data_bls: list[float]) -> dict[float, float]:
    g = df.groupby("branch_len")["mean_exact_J"].mean()
    return {float(bl): float(g.loc[bl]) for bl in data_bls if bl in g.index}


def build_lookup(df: pd.DataFrame) -> dict[tuple[str, float, float], float]:
    g = df.groupby(["Lstr", "branch_len", "size"], as_index=False)["MAE"].mean()
    return {
        (r.Lstr, float(r.branch_len), float(r.size)): float(r.MAE)
        for r in g.itertuples(index=False)
    }


def reg_ticks() -> tuple[list[float], list[str]]:
    ticks = [
        64, 128, 256, 512, 1024, 2048, 4096, 8192,
        16384, 32768, 65536, 131072, 262144, 524288,
        1048576, 2097152, 4194304, 8388608,
    ]
    labels = [
        "64", "128", "256", "512", "1k", "2k", "4k", "8k",
        "16k", "32k", "64k", "128k", "256k", "512k",
        "1M", "2M", "4M", "8M",
    ]
    return ticks, labels


def plot_grid(
    df: pd.DataFrame,
    lookup: dict[tuple[str, float, float], float],
    col_bls: list[tuple[float, float]],
    j_by_data_bl: dict[float, float],
    stem: str,
    png_only: bool,
) -> None:
    nL, nB = len(L_ORDER), len(col_bls)
    fig, axes = plt.subplots(
        nL, nB, figsize=(nB * 3.6, nL * 2.8), sharex=True, sharey=True, squeeze=False,
    )
    xticks, xlabels = reg_ticks()
    y_floor, y_ceil = 5e-4, 1.0

    for i, lstr in enumerate(L_ORDER):
        for j, (_, data_bl) in enumerate(col_bls):
            ax = axes[i, j]
            pts = sorted(
                (sz, lookup[(lstr, data_bl, sz)])
                for sz in {k[2] for k in lookup if k[0] == lstr and k[1] == data_bl}
                if (lstr, data_bl, sz) in lookup
            )
            if len(pts) >= 2:
                xs = [p[0] for p in pts]
                ys = [max(p[1], y_floor) for p in pts]
                ax.plot(xs, ys, "-o", ms=4.5, lw=2.0, color=COLOR, label=LABEL)
                for mark, ms in ((DEFAULT_LIB, 55), (DEFAULT_EQ, 70)):
                    mv = next((v for s, v in pts if s == mark), None)
                    if mv is not None:
                        ax.scatter(
                            [mark], [max(mv, y_floor)], s=ms,
                            facecolors="none", edgecolors=COLOR, linewidths=1.8, zorder=5,
                        )

            ax.axvline(DEFAULT_EQ, color="#cccccc", ls="--", lw=0.8, zorder=0)
            ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            ax.set_xlim(50, 10_000_000)
            ax.set_ylim(y_floor, y_ceil)
            ax.set_xticks(xticks)
            ax.grid(True, which="both", color="#eeeeee", lw=0.5)
            jv = df.loc[(df["Lstr"] == lstr) & (df["branch_len"] == data_bl), "mean_exact_J"].mean()
            if not np.isfinite(jv):
                jv = j_by_data_bl.get(data_bl, float("nan"))
            ax.set_title(
                f"$L$={L_DISPLAY[lstr]}, $J$={format_j_tick(jv)}" if np.isfinite(jv) else f"$L$={L_DISPLAY[lstr]}",
                pad=5, fontsize=FS_PANEL_TITLE,
            )
            if i == nL - 1:
                ax.set_xticklabels(xlabels, fontsize=FS_TICK_X_BOTTOM, rotation=45, ha="right")
                ax.set_xlabel("Sketch size", fontsize=FS_XLABEL, labelpad=8)
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel("MAE", fontsize=FS_YLABEL)

    fig.suptitle(
        "FastKMV sketch-size sensitivity (JC69, $k=21$; 64→8M; "
        "small circle = 1k default, large = 8k)",
        fontsize=FS_SUP, y=0.995,
    )
    fig.subplots_adjust(left=0.08, right=0.99, top=0.94, bottom=0.14, hspace=0.35, wspace=0.22)

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
    ap.add_argument("--ext-csv", default=str(EXT_CSV))
    ap.add_argument("--base-csv", default=str(BASE_CSV))
    ap.add_argument("--max-l", default="L2M")
    ap.add_argument("--stem", default="fastkmv_size_full_LxJ_uptoL2M")
    ap.add_argument("--png-only", action="store_true")
    args = ap.parse_args()

    df = load_merged(Path(args.base_csv), Path(args.ext_csv), args.max_l)
    if df.empty:
        raise SystemExit("no merged data")

    col_bls = column_bls(df)
    j_by_data_bl = j_map_for_bls(df, [db for _, db in col_bls])
    lookup = build_lookup(df)

    n_sc = df.groupby(["Lstr", "branch_len"]).ngroups
    print(
        f"merged {len(df)} rows, {n_sc} scenarios, L≤{args.max_l}, "
        f"sizes {int(df['size'].min())}..{int(df['size'].max())} ({df['size'].nunique()} points)",
        flush=True,
    )
    for target, data in col_bls:
        print(
            f"  col BL target={target} → data={data} ($J$={format_j_tick(j_by_data_bl.get(data, float('nan')))})",
            flush=True,
        )

    plot_grid(df, lookup, col_bls, j_by_data_bl, args.stem, args.png_only)


if __name__ == "__main__":
    main()
