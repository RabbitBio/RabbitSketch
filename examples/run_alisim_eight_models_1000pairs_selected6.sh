#!/usr/bin/env bash
# Run 6 algorithms (fastkmv, hll, setsketch_mle, bindash, probmh, kssd) on the
# eleven AliSim substitution models, 1000 pairs per scenario.
#
# Models (8 substitution + 3 indel):
#   JC69 HKY_ts2 HKY_ts5 GTR_ATrich GTR_G02 GTR_G05 GTR_G20 GTR_highGC
#   GTR_G_indel1 GTR_G_indel2 GTR_G_indel3
#
# Usage:
#   bash run_alisim_eight_models_1000pairs_selected6.sh [THREADS] [REFRESH] [N_PARALLEL]
#   THREADS      total CPU threads available (default 64)
#   REFRESH=1    remove per-model done markers and truncate CSVs before run
#   REFRESH=0    skip scenarios that already have a done marker (default)
#   N_PARALLEL   how many models to run concurrently (default 3;
#                  threads-per-model = THREADS / N_PARALLEL = 64/3 ≈ 21)
#
# Environment overrides:
#   SKIP_MODELS=JC69              comma-separated list of models to skip
#   ONLY_MODELS=GTR_G_indel1,...  run only these models (comma-separated)
#
# Example (full run, background):
#   nohup bash run_alisim_eight_models_1000pairs_selected6.sh 64 1 3 \
#     > results/all_models_pairs1000.log 2>&1 &
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${1:-64}"
REFRESH="${2:-0}"
N_PARALLEL="${3:-3}"
THREADS_PER_MODEL=$(( THREADS / N_PARALLEL ))
SKIP_MODELS="${SKIP_MODELS:-}"
ONLY_MODELS="${ONLY_MODELS:-}"

ALL_MODELS=(JC69 HKY_ts2 HKY_ts5 GTR_ATrich GTR_G02 GTR_G05 GTR_G20 GTR_highGC
            GTR_G_indel1 GTR_G_indel2 GTR_G_indel3)

if [[ -n "$ONLY_MODELS" ]]; then
    IFS=',' read -ra MODELS <<< "$ONLY_MODELS"
else
    MODELS=("${ALL_MODELS[@]}")
fi

should_skip_model() {
    local m="$1"
    [[ -z "$SKIP_MODELS" ]] && return 1
    local s
    for s in ${SKIP_MODELS//,/ }; do
        [[ "$s" == "$m" ]] && return 0
    done
    return 1
}

echo "=== AliSim sweep: ${#MODELS[@]} model(s), N_PARALLEL=$N_PARALLEL, THREADS_PER_MODEL=$THREADS_PER_MODEL ==="

# ---------------------------------------------------------------------------
# Phase 1: REFRESH cleanup (serial, so log is clean)
# ---------------------------------------------------------------------------
if [[ "$REFRESH" == "1" ]]; then
    echo ""
    echo "--- REFRESH: clearing markers and CSVs ---"
    for model in "${MODELS[@]}"; do
        should_skip_model "$model" && continue
        OUT_CSV="results/alisim_sweep.${model}.pairs1000.selected6.csv"
        MARKER_GLOB=".pairs1000_selected6_${model}_done"
        rm -f "$OUT_CSV"
        find /home/bigssd/zt/ref/"$model" -name "$MARKER_GLOB" -delete 2>/dev/null || true
        echo "  cleared $model"
    done
    echo "--- REFRESH done ---"
fi

# ---------------------------------------------------------------------------
# Phase 2: parallel model execution
#   At most N_PARALLEL models running simultaneously.
#   Each model writes its own log: results/<model>_pairs1000_selected6.log
#   'wait -n' (bash 4.3+) waits for any one background job to finish.
# ---------------------------------------------------------------------------
active_pids=()
active_models=()

for model in "${MODELS[@]}"; do
    if should_skip_model "$model"; then
        echo "=== Skip model $model (SKIP_MODELS) ==="
        continue
    fi

    # Wait until a slot opens up
    while (( ${#active_pids[@]} >= N_PARALLEL )); do
        wait -n 2>/dev/null || true      # wait for any one job to finish
        # Prune finished pids
        new_pids=(); new_models=()
        for idx in "${!active_pids[@]}"; do
            if kill -0 "${active_pids[$idx]}" 2>/dev/null; then
                new_pids+=("${active_pids[$idx]}")
                new_models+=("${active_models[$idx]}")
            else
                echo "=== Finished model: ${active_models[$idx]} ==="
            fi
        done
        active_pids=("${new_pids[@]+"${new_pids[@]}"}")
        active_models=("${new_models[@]+"${new_models[@]}"}")
    done

    LOG="results/${model}_pairs1000_selected6.log"
    mkdir -p results
    echo ""
    echo "========== Starting model: $model (log → $LOG) =========="
    bash run_alisim_one_model_1000pairs_selected_algos.sh "$model" "$THREADS_PER_MODEL" 0 \
        > "$LOG" 2>&1 &
    active_pids+=($!)
    active_models+=("$model")
done

# Wait for all remaining jobs
echo ""
echo "--- Waiting for remaining models: ${active_models[*]+"${active_models[*]}"} ---"
for pid in "${active_pids[@]+"${active_pids[@]}"}"; do
    wait "$pid" || echo "[WARN] model job pid=$pid exited non-zero"
done

echo ""
echo "All models finished at $(date -Iseconds)"
