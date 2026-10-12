#!/usr/bin/env python3
"""Remove rows for a given algo from sweep CSV(s). Writes backup if missing."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


def strip_file(path: Path, algo: str, backup: bool) -> tuple[int, int]:
    if not path.exists():
        return 0, 0
    bak = path.with_suffix(path.suffix + f".pre_{algo}_strip.bak")
    if backup and not bak.exists():
        shutil.copy2(path, bak)
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if not fieldnames:
        return 0, 0
    algo_col = "algo" if "algo" in fieldnames else None
    if not algo_col:
        raise SystemExit(f"{path}: no algo column")
    kept = [r for r in rows if r.get(algo_col) != algo]
    removed = len(rows) - len(kept)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(kept)
    return len(rows), removed


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--algo", required=True)
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("csv", nargs="+", type=Path)
    args = ap.parse_args()
    for p in args.csv:
        before, removed = strip_file(p, args.algo, backup=not args.no_backup)
        print(f"{p}: {before} rows, removed {removed} ({args.algo})")


if __name__ == "__main__":
    main()
