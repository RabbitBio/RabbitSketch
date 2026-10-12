#!/usr/bin/env bash
# Wait for SetSketch MLE backfill (792 fixed + 81 varlen), then run Mash-calib
# supplement (all algos at M≈0.075, 0.150, 0.200).
#
# Usage:
#   nohup bash run_after_mle_then_mash_calib.sh 64 1000 \
#     > results/after_mle_mash_calib.log 2>&1 &
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${1:-64}"
PAIRS_MAX="${2:-1000}"
LOG="${LOG:-results/after_mle_mash_calib.log}"
MARKER_DIR="/home/bigssd/zt/ref"
FIXED_NEED=792
VARLEN_NEED=81

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

log "Waiting for MLE backfill: ${FIXED_NEED} fixed + ${VARLEN_NEED} varlen (.mle_backfill_done)"

while true; do
    n=$(find "$MARKER_DIR" -name '.mle_backfill_done' 2>/dev/null | wc -l)
    prog=$(grep '^\[' results/setsketch_mle_backfill.log 2>/dev/null | tail -1 || true)
    log "  markers=$n / $((FIXED_NEED + VARLEN_NEED))  ${prog:-no log yet}"
    if (( n >= FIXED_NEED + VARLEN_NEED )); then
        log "MLE backfill complete."
        break
    fi
    if ! pgrep -f "run_setsketch_mle_backfill.sh" >/dev/null 2>&1; then
        if (( n < FIXED_NEED + VARLEN_NEED )); then
            log "WARNING: backfill process gone but markers=$n — continuing anyway."
        fi
        break
    fi
    sleep 120
done

log "Starting Mash calibration supplement (all algorithms)..."
bash run_mash_calib_supplement.sh "$THREADS" "$PAIRS_MAX" 2>&1 | tee -a "$LOG"

log "Done. Next: python3 build_master_csv.py && replot figures."
