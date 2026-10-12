"""Shared constants for publication figures.

SetSketch: always joint-MLE (Ertl 2021); default a=5, k=21, bits=13, base=2.
AliSim sweeps run the same k/bits grids at a ∈ {0.5, 1, 5, 10, 15, 20}.
Legacy inclusion–exclusion sketch rows are not used in master CSV or figures.
"""

from __future__ import annotations

import math

import matplotlib.colors as mcolors

SETSKETCH_A_DEFAULT = 5.0
SETSKETCH_BASE_DEFAULT = 2.0
SETSKETCH_BITS_DEFAULT = 13
SETSKETCH_K_DEFAULT = 21

SETSKETCH_PANEL_LABEL = "SetSketch"  # joint-MLE estimator

# Shared landscape figure styling (default + oracle heatmaps)
LANDSCAPE_ALGO_COLORS: dict[str, str] = {
    "hll": "#4D6271",
    "setsketch": "#7D96A2",
    "fastkmv": "#D0A47E",
    "probmh": "#B9825A",
    "bindash": "#8E5F4A",
    "kssd": "#5B8C83",
    "minhash": "#D0A47E",  # alias for legacy rows
}

LANDSCAPE_MAE_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "landscape_mae",
    ["#F7F6F4", "#ECE9E6", "#DDD6CF", "#C8B6A5", "#9BAFBA", "#6F8794", "#3F5563"],
    N=256,
)

LANDSCAPE_FIG_RC = {
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 13,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#333333",
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
}


def landscape_mae_norm(values) -> mcolors.LogNorm:
    import numpy as np

    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    vmin = max(5e-4, float(np.min(vals))) if len(vals) else 5e-4
    vmax = float(np.max(vals)) if len(vals) else 1.0
    return mcolors.LogNorm(vmin=vmin, vmax=vmax)


# Landscape x-axis: design Mash distance → AliSim BL (nearest empirical M at L10k, 11 models)
LANDSCAPE_MASH_COLUMNS: list[tuple[float, float]] = [
    (0.001, 0.0005),
    (0.005, 0.0025),
    (0.01, 0.005),
    (0.05, 0.025),
    (0.075, 0.038),
    (0.10, 0.053),
    (0.15, 0.080),
    (0.20, 0.116),
    (0.25, 0.150),
    (0.30, 0.200),
    (0.35, 0.280),
    (0.40, 0.31),
]

# Varlen landscape (3 models; lengths drawn from [L/4, L×4] per sequence)
VARLEN_LANDSCAPE_COLUMNS: list[tuple[float, float]] = [
    (0.01, 0.005),
    (0.075, 0.038),
    (0.15, 0.080),
    (0.20, 0.116),
    (0.25, 0.150),
    (0.35, 0.310),
]

# k/size sensitivity grids: pick columns by design Mash, label panels with exact J
SENSITIVITY_MASH_COLUMNS: list[tuple[float, float]] = [
    (0.01, 0.005),
    (0.05, 0.025),
    (0.10, 0.053),
    (0.20, 0.116),
    (0.40, 0.31),
]
# Legacy aliases
MASH_THRESHOLDS = [m for m, _ in SENSITIVITY_MASH_COLUMNS]
MASH_BL_ORDER = [bl for _, bl in SENSITIVITY_MASH_COLUMNS]


def sensitivity_mash_bls() -> list[float]:
    return [bl for _, bl in SENSITIVITY_MASH_COLUMNS]


def format_j_tick(j: float) -> str:
    """Panel/column label from pooled exact Jaccard similarity."""
    j = float(j)
    if j >= 0.1:
        return f"{j:.2f}"
    if j >= 0.01:
        return f"{j:.3f}".rstrip("0").rstrip(".")
    return f"{j:.3f}"


def j_mean_by_bl(df) -> dict[float, float]:
    """Mean exactJ_mean per AliSim BL (11-model average in master CSV)."""
    g = df.groupby("BL")["exactJ_mean"].mean()
    return {float(bl): float(v) for bl, v in g.items()}


def varlen_landscape_bls() -> list[float]:
    return [bl for _, bl in VARLEN_LANDSCAPE_COLUMNS]


def varlen_landscape_xtick_labels() -> list[str]:
    return [format_mash_tick(m) for m, _ in VARLEN_LANDSCAPE_COLUMNS]


def landscape_mash_bls() -> list[float]:
    return [bl for _, bl in LANDSCAPE_MASH_COLUMNS]


def format_mash_tick(m: float) -> str:
    """X tick = Mash distance only (no AliSim d tags)."""
    if m < 0.01:
        return f"{m:.3f}"
    if m < 0.1:
        return f"{m:.3f}".rstrip("0").rstrip(".")
    return f"{m:.2f}"


def landscape_mash_xtick_labels() -> list[str]:
    return [format_mash_tick(m) for m, _ in LANDSCAPE_MASH_COLUMNS]


def mash_distance(j: float, k: int = 21) -> float:
    j = max(float(j), 1e-300)
    return -math.log(2.0 * j / (1.0 + j)) / k


def mash_xtick_labels(bls: list[float], ej_by_bl, k: int = 21) -> list[str]:
    """Empirical Mash from exact J (used by sensitivity grids, not landscape heatmaps)."""
    labels = []
    for bl in bls:
        key = min(ej_by_bl.index, key=lambda x: abs(float(x) - bl))
        ej = float(ej_by_bl.loc[key]) if key in ej_by_bl.index else float(ej_by_bl[key])
        labels.append(format_mash_tick(mash_distance(ej, k)))
    return labels


def setsketch_param_mask(df, *, kmer: int | None = None, bits: int | None = None):
    """Boolean mask for SetSketch MLE rows at default a/base (optional k, bits)."""
    m = (
        (df["algo"] == "setsketch")
        & (df["ss_a"] == SETSKETCH_A_DEFAULT)
        & (df["ss_base"] == SETSKETCH_BASE_DEFAULT)
    )
    if kmer is not None:
        m = m & (df["kmer"] == kmer)
    if bits is not None:
        m = m & (df["bits"] == bits)
    return m
