#!/usr/bin/env bash
# Re-run selected algorithms on ONE model with fixed 1000 pairs per scenario.
#
# Algorithms (fixed): fastkmv,hll,setsketch_mle,bindash,probmh,kssd
# Usage:
#   bash run_alisim_one_model_1000pairs_selected_algos.sh [MODEL] [THREADS] [MAX_COMBOS]
# Example:
#   bash run_alisim_one_model_1000pairs_selected_algos.sh JC69 64 0
#
# Optimizations vs. v1:
#   1. _split_and_meta replaced by scripts/split_pairs_fast.py (16-thread I/O pool,
#      no per-pair subprocess spawning – saves ~150 s on every L8M scenario).
#   2. gen_alisim_batch_varlen precomputes ALL pair lengths with a single python3
#      call instead of 1000 individual calls (saves ~100 s per varlen scenario).
# Temp iqtree output goes under /home/bigssd/zt/tmp (NOT /tmp on root — L8M peaks ~32 GB).
set -euo pipefail
cd "$(dirname "$0")"

MODEL_FILTER="${1:-JC69}"
THREADS="${2:-64}"
MAX_COMBOS="${3:-0}"
FORCE_PAIRS=1000
ALGOS="fastkmv,hll,setsketch_mle,bindash,probmh,kssd"

OUTBASE="/home/bigssd/zt/ref"
TMPBASE="/home/bigssd/zt/tmp"
RESDIR="results"
DATASET_CSV="$RESDIR/alisim_dataset_summary.csv"
SWEEP_EXE="./exe_eval_alisim_sweep"
OUT_CSV="$RESDIR/alisim_sweep.${MODEL_FILTER}.pairs1000.selected6.csv"
MARKER_NAME=".pairs1000_selected6_${MODEL_FILTER}_done"

SPLIT_PY="$(dirname "$0")/scripts/split_pairs_fast.py"
[[ -f "$SPLIT_PY" ]] || { echo "[ERROR] $SPLIT_PY not found"; exit 1; }

for _iqbin in "${HOME}/anaconda3/envs/cpp20/bin" \
               "/home/user_home/tongz/anaconda3/envs/cpp20/bin"; do
    [[ -d "$_iqbin" ]] && export PATH="$_iqbin:$PATH" && break
done

declare -A MODEL_STR MODEL_INDEL
MODEL_STR["JC69"]="JC"
MODEL_STR["HKY_ts2"]="HKY{2.0}"
MODEL_STR["HKY_ts5"]="HKY{5.0}"
MODEL_STR["GTR_ATrich"]="GTR{1.5,4.0,1.2,0.8,3.5,1.0}+F{0.3,0.2,0.2,0.3}"
MODEL_STR["GTR_G02"]="GTR+G{0.2}"
MODEL_STR["GTR_G05"]="GTR+G{0.5}"
MODEL_STR["GTR_G20"]="GTR+G{2.0}"
MODEL_STR["GTR_highGC"]="GTR+G{1.0}+F{0.15,0.35,0.35,0.15}"
MODEL_STR["GTR_G_indel1"]="GTR+G{1.0}"
MODEL_STR["GTR_G_indel2"]="GTR+G{1.0}"
MODEL_STR["GTR_G_indel3"]="GTR+G{1.0}"
MODEL_INDEL["GTR_G_indel1"]="0.005,0.005"
MODEL_INDEL["GTR_G_indel2"]="0.01,0.01"
MODEL_INDEL["GTR_G_indel3"]="0.03,0.03"

mkdir -p "$RESDIR" "$OUTBASE" "$TMPBASE"
[[ -x "$SWEEP_EXE" ]] || { echo "[ERROR] $SWEEP_EXE missing"; exit 1; }

# Remove stale iqtree scratch dirs from interrupted runs (older than 1 day).
find "$TMPBASE" -maxdepth 1 -type d -name 'alisim_1k6_*' -mtime +1 -exec rm -rf {} + 2>/dev/null || true
find "$TMPBASE" -maxdepth 1 -type d -name 'alisim_1k6_var_*' -mtime +1 -exec rm -rf {} + 2>/dev/null || true

IQTREE=""
for cand in iqtree2 iqtree3 iqtree; do
    if p="$(command -v "$cand" 2>/dev/null)" && [[ -n "$p" && -x "$p" ]]; then IQTREE="$p"; break; fi
done
[[ -n "$IQTREE" && -x "$IQTREE" ]] || { echo "[ERROR] iqtree not found"; exit 1; }

# ---------------------------------------------------------------------------
# gen_alisim_batch_fixed
#   Runs n_pairs iqtree simulations (up to THREADS in parallel), then calls
#   the fast Python splitter in one shot.
# ---------------------------------------------------------------------------
gen_alisim_batch_fixed() {
    local model_name="$1" iqt_model="$2" indel_rates="$3"
    local len="$4" bl="$5" outdir="$6" n_pairs="$7"
    local tmpd="${TMPBASE}/alisim_1k6_${$}_${model_name}_${len}_${bl//./p}"
    mkdir -p "$tmpd" "$outdir"
    trap 'rm -rf "$tmpd"' RETURN
    printf "(A:%.6f,B:%.6f);\n" "$bl" "$bl" > "$tmpd/tree2.nwk"
    local pids=() slot=0
    for ((i=0; i<n_pairs; i++)); do
        local cmd=("$IQTREE" --alisim "${tmpd}/p${i}_sim" -m "$iqt_model" -t "$tmpd/tree2.nwk"
            --length "$len" --out-format fasta --seed $((42+i)) -nt 1 -quiet)
        [[ -n "$indel_rates" ]] && cmd+=(--indel "$indel_rates")
        "${cmd[@]}" >/dev/null 2>&1 &
        pids+=($!)
        slot=$((slot + 1))
        if (( slot >= THREADS )); then
            wait "${pids[0]}" || true; pids=("${pids[@]:1}"); slot=$((slot - 1))
        fi
    done
    for p in "${pids[@]}"; do wait "$p" || true; done
    python3 "$SPLIT_PY" "$tmpd" "$outdir" "$bl" "$n_pairs"
    rm -rf "$tmpd"
    trap - RETURN
}

# ---------------------------------------------------------------------------
# gen_alisim_batch_varlen
#   Same as fixed, but each pair gets a log-uniform random length.
#   Optimization: ALL lengths are precomputed with ONE python3 call (instead of
#   one call per pair), saving ~100 s of interpreter startup overhead per scenario.
# ---------------------------------------------------------------------------
gen_alisim_batch_varlen() {
    local model_name="$1" iqt_model="$2"
    local nominal_len="$3" bl="$4" outdir="$5" n_pairs="$6"
    local tmpd="${TMPBASE}/alisim_1k6_var_${$}_${model_name}_${nominal_len}_${bl//./p}"
    mkdir -p "$tmpd" "$outdir"
    trap 'rm -rf "$tmpd"' RETURN
    printf "(A:%.6f,B:%.6f);\n" "$bl" "$bl" > "$tmpd/tree2.nwk"

    # Precompute all pair lengths in one Python call
    IFS=' ' read -ra VARLEN_SIZES <<< "$(python3 - <<PYEOF
import random, math
results = []
nom = $nominal_len
for i in range($n_pairs):
    random.seed(i * 7919 + 31337)
    lo = max(200, nom // 4)
    hi = min(nom * 4, 8_000_000)
    if hi <= lo:
        results.append(lo)
    else:
        results.append(int(math.exp(random.uniform(math.log(lo), math.log(hi)))))
print(' '.join(map(str, results)))
PYEOF
)"

    local pids=() slot=0
    for ((i=0; i<n_pairs; i++)); do
        local pair_len="${VARLEN_SIZES[$i]}"
        local cmd=("$IQTREE" --alisim "${tmpd}/p${i}_sim" -m "$iqt_model" -t "$tmpd/tree2.nwk"
            --length "$pair_len" --out-format fasta --seed $((42+i)) -nt 1 -quiet)
        "${cmd[@]}" >/dev/null 2>&1 &
        pids+=($!)
        slot=$((slot + 1))
        if (( slot >= THREADS )); then
            wait "${pids[0]}" || true; pids=("${pids[@]:1}"); slot=$((slot - 1))
        fi
    done
    for p in "${pids[@]}"; do wait "$p" || true; done
    python3 "$SPLIT_PY" "$tmpd" "$outdir" "$bl" "$n_pairs"
    rm -rf "$tmpd"
    trap - RETURN
}

run_selected_algos() {
    local outdir="$1" model_name="$2"
    "$SWEEP_EXE" --dir "$outdir" --model "$model_name" --threads "$THREADS" --algos "$ALGOS" \
      | awk -v csv="$OUT_CSV" '
        /^ALISIM_SWEEP/ {
            dir=""; model=""; bl=""; algo=""; kmer=""; size=""; bits=""
            ss_a=""; ss_base=""; bd_bbits=""; bd_s64=""; drlevel=""
            n=""; mej=""; mest=""; mae=""; rmse=""
            for(i=2;i<=NF;i++){
                split($i,kv,"="); k=kv[1]; v=kv[2]
                if(k=="dir")          dir=v
                if(k=="model")        model=v
                if(k=="branch_len")   bl=v
                if(k=="algo")         algo=v
                if(k=="kmer")         kmer=v
                if(k=="size")         size=v
                if(k=="bits")         bits=v
                if(k=="ss_a")         ss_a=v
                if(k=="ss_base")      ss_base=v
                if(k=="bd_bbits")     bd_bbits=v
                if(k=="bd_s64")       bd_s64=v
                if(k=="drlevel")      drlevel=v
                if(k=="n")            n=v
                if(k=="mean_exact_J") mej=v
                if(k=="mean_est_J")   mest=v
                if(k=="MAE")          mae=v
                if(k=="RMSE")         rmse=v
            }
            printf "%s,%s,%.6f,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
                dir,model,bl,algo,kmer,size,bits,ss_a,ss_base,
                bd_bbits,bd_s64,drlevel,n,mej,mest,mae,rmse >> csv
        }'
}

if [[ ! -f "$OUT_CSV" ]]; then
    echo "dir,model,branch_len,algo,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,n,mean_exact_J,mean_est_J,MAE,RMSE" > "$OUT_CSV"
fi

n_total=$(awk -F, -v m="$MODEL_FILTER" 'NR>1 && $3==m {c++} END{print c+0}' "$DATASET_CSV")
echo "One-model selected6 run: MODEL=$MODEL_FILTER FORCE_PAIRS=$FORCE_PAIRS combos=$n_total"
echo "ALGOS=$ALGOS"
echo "CSV=$OUT_CSV"

n=0
while IFS=, read -r dataset len_mode model nominal_len branch_len len_tag bl_tag _pairs _rest; do
    [[ "$dataset" == "dataset" ]] && continue
    [[ "$model" != "$MODEL_FILTER" ]] && continue

    outdir="$OUTBASE/${model}/${len_mode}_${len_tag}/${bl_tag}"
    marker="$outdir/$MARKER_NAME"
    if [[ -f "$marker" ]]; then
        echo "  SKIP $dataset"
        continue
    fi

    n=$((n + 1))
    if (( MAX_COMBOS > 0 && n > MAX_COMBOS )); then
        echo "MAX_COMBOS=$MAX_COMBOS reached"
        break
    fi

    n_pairs="$FORCE_PAIRS"
    echo "[$n/$n_total] RUN $dataset  n=$n_pairs"
    mkdir -p "$outdir"
    iqt="${MODEL_STR[$model]:-}"
    [[ -n "$iqt" ]] || { echo "  [WARN] unknown model $model"; continue; }

    cleanup_scenario_data() {
        rm -f "$outdir"/pair_*.fa "$outdir"/pairs_meta.tsv 2>/dev/null || true
    }
    trap cleanup_scenario_data RETURN

    if [[ "$len_mode" == "fixed" ]]; then
        gen_alisim_batch_fixed "$model" "$iqt" "${MODEL_INDEL[$model]:-}" \
            "$nominal_len" "$branch_len" "$outdir" "$n_pairs"
    else
        gen_alisim_batch_varlen "$model" "$iqt" "$nominal_len" "$branch_len" "$outdir" "$n_pairs"
    fi

    run_selected_algos "$outdir" "$model"
    cleanup_scenario_data
    trap - RETURN
    touch "$marker"
done < "$DATASET_CSV"

echo "Done at $(date -Iseconds)"
