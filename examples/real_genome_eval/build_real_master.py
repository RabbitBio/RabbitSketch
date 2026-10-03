#!/usr/bin/env python3
"""
Build a long-format CSV combining:
  pairs.tsv            (200 pair definitions with tier/taxonomy)
  fastani_truth.tsv    (ground truth)
  sketch_results.csv   (per-pair × config sketch outputs)

Output:
  real_genome_master.csv   — one row per (pair × config), with FastANI ANI joined
"""

import csv
from pathlib import Path

HERE = Path(__file__).parent

pairs = {}
with open(HERE / "pairs.tsv") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for r in rdr:
        pairs[r["pair_id"]] = r

truth = {}
with open(HERE / "fastani_truth.tsv") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for r in rdr:
        v = float(r["fastani"])
        truth[r["pair_id"]] = {
            "ani":   None if v < 0 else v / 100.0,
            "frac":  float(r["frac_aligned"]),
            "total": int(float(r["total_frags"])),
        }

out_path = HERE / "real_genome_master.csv"
cols = ["pair_id", "tier", "accA", "accB", "orgA", "orgB",
        "fastani_ani", "fastani_frac_aligned", "fastani_total_frags",
        "layer", "tag", "params_str", "mash_dist", "sketch_ani",
        "abs_error", "sec"]

n = 0
with open(out_path, "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=cols)
    w.writeheader()
    with open(HERE / "sketch_results.csv") as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            pid = r["pair_id"]
            if pid not in pairs: continue
            t = truth.get(pid, {"ani": None, "frac": 0.0, "total": 0})
            ani_truth = t["ani"]
            ani_sketch = float(r["sketch_ani"])
            err = abs(ani_sketch - ani_truth) if ani_truth is not None else ""
            w.writerow({
                "pair_id":          pid,
                "tier":             pairs[pid]["tier"],
                "accA":             pairs[pid]["accA"],
                "accB":             pairs[pid]["accB"],
                "orgA":             pairs[pid]["orgA"],
                "orgB":             pairs[pid]["orgB"],
                "fastani_ani":      f"{ani_truth:.6f}" if ani_truth is not None else "",
                "fastani_frac_aligned": f"{t['frac']:.4f}",
                "fastani_total_frags":  t["total"],
                "layer":            r["layer"],
                "tag":              r["algo"],
                "params_str":       r["params"],
                "mash_dist":        f"{float(r['mash_dist']):.6f}",
                "sketch_ani":       f"{ani_sketch:.6f}",
                "abs_error":        f"{err:.6f}" if err != "" else "",
                "sec":              f"{float(r['sec']):.3f}",
            })
            n += 1

print(f"Wrote {out_path}  ({n} rows)")
