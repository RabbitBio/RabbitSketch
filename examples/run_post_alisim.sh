#!/usr/bin/env bash
# Post-processing pipeline: run AFTER results/alisim_sweep.csv is final.
#
#   1) Build master CSVs   (alisim_master, alisim_defaults)
#   2) Build real-genome master CSV (idempotent, requires sketch_results.csv)
#   3) Run textual analysis  → results/analysis_report.txt
#   4) Generate 12 figures   → results/figures/figNN_*.png
#
# Usage:  bash run_post_alisim.sh
set -euo pipefail
cd "$(dirname "$0")"

CSV=results/alisim_sweep.csv
if [[ ! -s "$CSV" ]]; then
    echo "Error: $CSV missing or empty.  Run run_alisim_eval.sh first."
    exit 1
fi

echo "════════════════════════════════════════════════════════════════════════"
echo "  [1/4] Building master CSVs"
echo "════════════════════════════════════════════════════════════════════════"
python3 build_master_csv.py

echo
echo "════════════════════════════════════════════════════════════════════════"
echo "  [2/4] Building real-genome master CSV"
echo "════════════════════════════════════════════════════════════════════════"
if [[ -s real_genome_eval/sketch_results.csv ]]; then
    python3 real_genome_eval/build_real_master.py
else
    echo "  skipped: real_genome_eval/sketch_results.csv missing"
fi

echo
echo "════════════════════════════════════════════════════════════════════════"
echo "  [3/4] Running textual analysis"
echo "════════════════════════════════════════════════════════════════════════"
python3 analyze_alisim_results.py > results/analysis_report.txt
echo "  wrote results/analysis_report.txt  ($(wc -l < results/analysis_report.txt) lines)"

echo
echo "════════════════════════════════════════════════════════════════════════"
echo "  [4/4] Generating 12 figures"
echo "════════════════════════════════════════════════════════════════════════"
python3 make_all_figures.py

echo
echo "════════════════════════════════════════════════════════════════════════"
echo "  Done."
echo "  Outputs:"
echo "    results/alisim_master.csv        — all parameter sweeps (aggregated)"
echo "    results/alisim_defaults.csv      — default-param subset"
echo "    results/analysis_report.txt      — human-readable analysis"
echo "    results/figures/figNN_*.png      — 12 figures"
echo "    real_genome_eval/real_genome_master.csv  — real bacterial pairs joined"
echo "════════════════════════════════════════════════════════════════════════"
