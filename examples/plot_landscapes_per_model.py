#!/usr/bin/env python3
"""Per evolutionary model: default + oracle accuracy landscapes (PNG only).

The pooled figures average over 11 models; this script emits one pair per model.

Outputs:
  results/figures/per_model/<MODEL>/accuracy_landscape_default.png
  results/figures/per_model/<MODEL>/accuracy_landscape_oracle.png
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SWEEP = HERE / "results" / "alisim_sweep.csv"
PER_MODEL_MASTER = HERE / "results" / "per_model"
FIG_ROOT = HERE / "results" / "figures" / "per_model"

MODEL_ORDER = [
    "JC69",
    "HKY_ts2",
    "HKY_ts5",
    "GTR_ATrich",
    "GTR_G02",
    "GTR_G05",
    "GTR_G20",
    "GTR_highGC",
    "GTR_G_indel1",
    "GTR_G_indel2",
    "GTR_G_indel3",
]


def models_in_sweep() -> list[str]:
    with SWEEP.open() as f:
        names = sorted({row["model"] for row in csv.DictReader(f)})
    order = {m: i for i, m in enumerate(MODEL_ORDER)}
    return sorted(names, key=lambda m: order.get(m, 999))


def build_model_master(model: str, rebuild: bool) -> Path:
    PER_MODEL_MASTER.mkdir(parents=True, exist_ok=True)
    master = PER_MODEL_MASTER / f"alisim_master_{model}.csv"
    if master.is_file() and not rebuild:
        return master
    subprocess.run(
        [
            sys.executable,
            str(HERE / "build_master_csv.py"),
            "--models",
            model,
            "--master",
            str(master),
            "--defaults",
            str(PER_MODEL_MASTER / f"alisim_defaults_{model}.csv"),
        ],
        check=True,
        cwd=HERE,
    )
    return master


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--models",
        default="",
        help="Comma-separated subset (default: all 11 models in sweep)",
    )
    ap.add_argument("--rebuild-master", action="store_true")
    ap.add_argument("--no-kssd-oracle", action="store_true",
                    help="Oracle winner/panels exclude KSSD (like oracle_no_kssd)")
    args = ap.parse_args()

    if args.models.strip():
        models = [m.strip() for m in args.models.split(",") if m.strip()]
    else:
        models = models_in_sweep()

    from plot_accuracy_landscape import plot_landscape
    from plot_accuracy_landscape_oracle import ALGOS_ALL, plot_oracle

    algos = [a for a in ALGOS_ALL if not (args.no_kssd_oracle and a == "kssd")]

    print(f"Generating per-model landscapes for {len(models)} models (PNG only)")
    for model in models:
        print(f"\n=== {model} ===")
        master = build_model_master(model, rebuild=args.rebuild_master)
        out_dir = FIG_ROOT / model
        plot_landscape(
            "MAE_mean",
            "accuracy_landscape_default",
            bindash_mem="2k",
            subtitle=f"{model}  |  Default params  |  SetSketch=MLE a=5",
            master_csv=master,
            out_dir=out_dir,
            png_only=True,
        )
        plot_oracle(
            "MAE_mean",
            "accuracy_landscape_oracle",
            algos=algos,
            master_csv=master,
            subtitle=f"{model}  |  Oracle (best params)  |  SetSketch=MLE",
            out_dir=out_dir,
            png_only=True,
            save_configs=False,
        )

    print(f"\nDone. Figures under {FIG_ROOT}/")


if __name__ == "__main__":
    main()
