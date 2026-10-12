#!/usr/bin/env python3
"""Aggregate results/alisim_sweep.csv → alisim_master.csv + alisim_defaults.csv.

SetSketch = joint-MLE only (sweep algo ``setsketch_mle`` → ``setsketch``).
Inclusion–exclusion ``setsketch`` rows in the sweep are dropped.
Default params: k=21, bits=13, a=5, base=2.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SWEEP = HERE / "results" / "alisim_sweep.csv"
MASTER = HERE / "results" / "alisim_master.csv"
DEFAULTS = HERE / "results" / "alisim_defaults.csv"


def parse_L(d: str) -> tuple[str, int, str]:
    parts = d.strip("/").split("/")
    md_L = parts[-2]
    mode, Ls = md_L.split("_", 1)
    s = Ls.lstrip("L")
    if s.endswith("M"):
        L = int(float(s[:-1]) * 1_000_000)
    elif s.endswith("k"):
        L = int(float(s[:-1]) * 1_000)
    else:
        L = int(s)
    return mode, L, Ls


def normalize_sweep_row(row: dict) -> dict | None:
    """SetSketch ≡ joint-MLE: rename setsketch_mle; drop legacy inclusion–exclusion setsketch."""
    algo = row["algo"]
    if algo == "setsketch_mle":
        row = dict(row)
        row["algo"] = "setsketch"
        return row
    if algo == "setsketch":
        return None
    return row


def is_default(algo: str, kmer: int, size: int, bits: int, ss_a: float, ss_base: float,
               bd_bbits: int, bd_s64: int, drlevel: int) -> bool:
    if algo == "minhash":
        return kmer == 21 and size == 1024
    if algo == "fastkmv":
        return kmer == 21 and size == 1024
    if algo == "hll":
        return kmer == 21 and bits == 13
    if algo == "setsketch":
        return kmer == 21 and bits == 13 and ss_a == 5.0 and ss_base == 2.0
    if algo == "bindash":
        # s64=16 → 16×64=1024 bins, ~2 KB (library default; not s64=64 ~8 KB)
        return kmer == 21 and bd_bbits == 16 and bd_s64 == 16
    if algo == "kssd":
        return drlevel == 3
    if algo == "probmh":
        return kmer == 21 and size == 1024
    return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--sweep",
        default="",
        help="Input sweep CSV (default: results/alisim_sweep.csv)",
    )
    ap.add_argument(
        "--models",
        default="",
        help="Comma-separated model names (e.g. JC69,HKY_ts2). Default: all models in sweep.",
    )
    ap.add_argument(
        "--master",
        default="",
        help="Output master CSV path (default: results/alisim_master.csv)",
    )
    ap.add_argument(
        "--defaults",
        default="",
        help="Output defaults CSV path (default: results/alisim_defaults.csv)",
    )
    args = ap.parse_args()

    model_filter: set[str] | None = None
    if args.models.strip():
        model_filter = {m.strip() for m in args.models.split(",") if m.strip()}

    master_path = Path(args.master) if args.master else MASTER
    defaults_path = Path(args.defaults) if args.defaults else DEFAULTS
    sweep_path = Path(args.sweep) if args.sweep else SWEEP

    if not sweep_path.is_file():
        raise SystemExit(f"missing {sweep_path}")

    groups: dict[tuple, list[dict]] = defaultdict(list)
    with sweep_path.open() as f:
        for row in csv.DictReader(f):
            if model_filter is not None and row["model"] not in model_filter:
                continue
            row = normalize_sweep_row(row)
            if row is None:
                continue
            try:
                mode, L, Lstr = parse_L(row["dir"])
                groups[(
                    mode,
                    L,
                    Lstr,
                    float(row["branch_len"]),
                    row["algo"],
                    int(row["kmer"]),
                    int(row["size"]),
                    int(row["bits"]),
                    float(row["ss_a"]),
                    float(row["ss_base"]),
                    int(row["bd_bbits"]),
                    int(row["bd_s64"]),
                    int(row["drlevel"]),
                )].append({
                    "model": row["model"],
                    "n": int(row["n"]),
                    "exactJ": float(row["mean_exact_J"]),
                    "estJ": float(row["mean_est_J"]),
                    "mae": float(row["MAE"]),
                    "rmse": float(row["RMSE"]) if row.get("RMSE") not in ("", None) else None,
                })
            except Exception:
                continue

    header = [
        "mode", "L", "Lstr", "BL", "algo", "kmer", "size", "bits",
        "ss_a", "ss_base", "bd_bbits", "bd_s64", "drlevel", "is_default",
        "n_models", "n_pairs_mean", "exactJ_mean", "estJ_mean",
        "MAE_mean", "MAE_std", "RMSE_mean", "bias_mean",
    ]
    out_rows = []
    for key, items in sorted(groups.items()):
        (mode, L, Lstr, BL, algo, kmer, size, bits,
         ss_a, ss_base, bd_bbits, bd_s64, drlevel) = key
        maes = [x["mae"] for x in items]
        rmses = [x["rmse"] for x in items if x["rmse"] is not None]
        ej = [x["exactJ"] for x in items]
        est = [x["estJ"] for x in items]
        bias = [x["estJ"] - x["exactJ"] for x in items]
        n_pairs = statistics.mean(x["n"] for x in items)
        rmse_mean = statistics.mean(rmses) if rmses else ""
        lstr_out = Lstr.lstrip("L") if str(Lstr).startswith("L") else Lstr
        out_rows.append({
            "mode": mode,
            "L": L,
            "Lstr": lstr_out,
            "BL": BL,
            "algo": algo,
            "kmer": kmer,
            "size": size,
            "bits": bits,
            "ss_a": ss_a,
            "ss_base": ss_base,
            "bd_bbits": bd_bbits,
            "bd_s64": bd_s64,
            "drlevel": drlevel,
            "is_default": int(is_default(algo, kmer, size, bits, ss_a, ss_base,
                                          bd_bbits, bd_s64, drlevel)),
            "n_models": len(items),
            "n_pairs_mean": n_pairs,
            "exactJ_mean": statistics.mean(ej),
            "estJ_mean": statistics.mean(est),
            "MAE_mean": statistics.mean(maes),
            "MAE_std": statistics.pstdev(maes) if len(maes) > 1 else 0.0,
            "RMSE_mean": rmse_mean,
            "bias_mean": statistics.mean(bias),
        })

    master_path.parent.mkdir(parents=True, exist_ok=True)
    with master_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(out_rows)

    def_rows = [r for r in out_rows if r["is_default"]]
    with defaults_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(def_rows)

    if model_filter:
        print(f"models: {', '.join(sorted(model_filter))}")
    print(f"wrote {master_path}  ({len(out_rows)} rows)")
    print(f"wrote {defaults_path}  ({len(def_rows)} default rows)")
    bls = sorted({r["BL"] for r in def_rows if r["mode"] == "fixed"})
    print(f"fixed default BL levels: {len(bls)}  ({bls[:3]} ... {bls[-3:]})")


if __name__ == "__main__":
    main()
