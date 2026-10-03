#!/usr/bin/env bash
# One-button driver for the real-bacterial-genome benchmark.
#
# Pipeline:
#   1) select_real_pairs.py   →  pairs.tsv          (~200 pairs)
#   2) run_fastani_truth.sh   →  fastani_truth.tsv  (ground truth ANI)
#   3) run_sketch_real.sh     →  sketch_results.csv (7 algos × default + spot-check)
#   4) analyze_real_results.py → printed report
#
# Total wall-clock estimate (16-way parallelism):
#   pair selection      <  1 s
#   FastANI             ~ 10-20 min  (FastANI is the slowest step)
#   sketching           ~  5-10 min
#   analysis            <  5 s
# Disk: pair-level temp files only; no persistent intermediate sketches.
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${THREADS:-16}"

echo "================================================================="
echo "  Real bacterial-genome benchmark   (THREADS=$THREADS)"
echo "  $(date)"
echo "================================================================="

# 1) pair selection (idempotent)
if [[ ! -s pairs.tsv ]]; then
    echo "[1/4] Selecting pairs…"
    python3 select_real_pairs.py
else
    echo "[1/4] pairs.tsv already exists ($(($(wc -l < pairs.tsv) - 1)) pairs)"
fi

# 2) FastANI truth (idempotent if user keeps the output file)
if [[ ! -s fastani_truth.tsv ]]; then
    echo ""
    echo "[2/4] FastANI ground truth (this is the slowest stage)…"
    THREADS="$THREADS" bash run_fastani_truth.sh
else
    echo "[2/4] fastani_truth.tsv already exists  ($(($(wc -l < fastani_truth.tsv) - 1)) rows)"
fi

# 3) Sketch runs (default + spot-check)
if [[ ! -s sketch_results.csv ]]; then
    echo ""
    echo "[3/4] Sketch benchmarks (Layer-1 defaults + Layer-2 spot-checks)…"
    PARALLEL="$THREADS" bash run_sketch_real.sh
else
    echo "[3/4] sketch_results.csv already exists ($(($(wc -l < sketch_results.csv) - 1)) rows)"
fi

# 4) Analysis
echo ""
echo "[4/4] Analysis report:"
echo ""
python3 analyze_real_results.py

echo ""
echo "================================================================="
echo "  Real-genome benchmark done at $(date)"
echo "  Outputs:"
echo "    pairs.tsv             — pair list with tier"
echo "    fastani_truth.tsv     — FastANI ANI ground truth"
echo "    sketch_results.csv    — per-(pair × config) sketch ANI"
echo "================================================================="
