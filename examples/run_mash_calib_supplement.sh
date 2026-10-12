#!/usr/bin/env bash
# Full algorithm sweep at 4 Mash-calibrated branch lengths (k=21 targets):
#   d=0.038  →  M≈0.075
#   d=0.080  →  M≈0.150
#   d=0.116  →  M≈0.200
#   d=0.31   →  M≈0.400  (new; existing d=0.28/0.38 are M≈0.35/0.39)
#
# Fixed: 11 models × 9 lengths × 4 BL = 396 combos
# Varlen: 3 models × 9 lengths × 4 BL = 108 combos
# Appends to results/alisim_sweep.csv (all algos incl. setsketch + setsketch_mle from sweep exe)
#
# Usage:
#   bash run_mash_calib_supplement.sh [THREADS] [PAIRS_MAX]
#   nohup bash run_mash_calib_supplement.sh 64 1000 > results/mash_calib_supplement.log 2>&1 &
set -euo pipefail
cd "$(dirname "$0")"

for _iqbin in "${HOME}/anaconda3/envs/cpp20/bin" \
               "/home/user_home/tongz/anaconda3/envs/cpp20/bin"; do
    [[ -d "$_iqbin" ]] && export PATH="$_iqbin:$PATH" && break
done

THREADS="${1:-64}"
PAIRS_MAX="${2:-1000}"
LOG="${LOG:-results/mash_calib_supplement.log}"

# shellcheck source=alisim_bl_config.sh
source "$(dirname "$0")/alisim_bl_config.sh"

echo "══════════════════════════════════════════════════════════════════"
echo " Mash calibration supplement (all algorithms)"
echo " branch_len: ${MASH_CALIB_BRANCH_LENS[*]}"
echo " targets M:  ${MASH_CALIB_TARGET_M[*]}"
echo " THREADS=$THREADS  PAIRS_MAX=$PAIRS_MAX"
echo " Log: $LOG"
echo "══════════════════════════════════════════════════════════════════"

export MASH_CALIB_SUPPLEMENT=1
mkdir -p results
bash run_alisim_eval.sh "$THREADS" "$PAIRS_MAX"
