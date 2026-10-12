#!/usr/bin/env bash
# FastKMV size extension on JC69: k=21, sketch sizes 8k→8M (doubling).
#
# Sizes (registers): 8192, 16384, 32768, 65536, 131072, 262144,
#                    524288, 1048576, 2097152, 4194304, 8388608
#
# Regenerates 1000 pairs per scenario, evaluates FastKMV only, appends CSV.
# Temp data under /home/bigssd/zt/tmp (not /tmp).
#
# Usage:
#   bash run_jc69_fastkmv_size_ext_backfill.sh [THREADS] [MAX_COMBOS]
#   nohup bash run_jc69_fastkmv_size_ext_backfill.sh 64 0 \
#     > results/jc69_fastkmv_size_ext.log 2>&1 &
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${1:-64}"
MAX_COMBOS="${2:-0}"
MODEL_FILTER="JC69"
FORCE_PAIRS=1000

OUTBASE="/home/bigssd/zt/ref"
TMPBASE="/home/bigssd/zt/tmp"
RESDIR="results"
DATASET_CSV="$RESDIR/alisim_dataset_summary.csv"
FKMV_EXE="./exe_eval_fastkmv_size_ext"
OUT_CSV="$RESDIR/alisim_sweep.JC69.pairs1000.fastkmv_size_ext.csv"
MARKER_NAME=".fastkmv_size_ext_JC69_done"
SPLIT_PY="$(dirname "$0")/scripts/split_pairs_fast.py"

for _iqbin in "${HOME}/anaconda3/envs/cpp20/bin" \
               "/home/user_home/tongz/anaconda3/envs/cpp20/bin"; do
    [[ -d "$_iqbin" ]] && export PATH="$_iqbin:$PATH" && break
done

declare -A MODEL_STR
MODEL_STR["JC69"]="JC"

mkdir -p "$RESDIR" "$OUTBASE" "$TMPBASE"
make -s exe_eval_fastkmv_size_ext
[[ -x "$FKMV_EXE" ]] || { echo "[ERROR] $FKMV_EXE missing"; exit 1; }
[[ -f "$SPLIT_PY" ]] || { echo "[ERROR] $SPLIT_PY not found"; exit 1; }

IQTREE=""
for cand in iqtree2 iqtree3 iqtree; do
    if p="$(command -v "$cand" 2>/dev/null)" && [[ -n "$p" && -x "$p" ]]; then IQTREE="$p"; break; fi
done
[[ -n "$IQTREE" && -x "$IQTREE" ]] || { echo "[ERROR] iqtree not found"; exit 1; }

gen_alisim_batch_fixed() {
    local model_name="$1" iqt_model="$2" len="$3" bl="$4" outdir="$5" n_pairs="$6"
    local tmpd="${TMPBASE}/alisim_fkmv_ext_${$}_${model_name}_${len}_${bl//./p}"
    mkdir -p "$tmpd" "$outdir"
    trap 'rm -rf "$tmpd"' RETURN
    printf "(A:%.6f,B:%.6f);\n" "$bl" "$bl" > "$tmpd/tree2.nwk"
    local pids=() slot=0
    for ((i=0; i<n_pairs; i++)); do
        local cmd=("$IQTREE" --alisim "${tmpd}/p${i}_sim" -m "$iqt_model" -t "$tmpd/tree2.nwk"
            --length "$len" --out-format fasta --seed $((42+i)) -nt 1 -quiet)
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

gen_alisim_batch_varlen() {
    local model_name="$1" iqt_model="$2" nominal_len="$3" bl="$4" outdir="$5" n_pairs="$6"
    local tmpd="${TMPBASE}/alisim_fkmv_ext_var_${$}_${model_name}_${nominal_len}_${bl//./p}"
    mkdir -p "$tmpd" "$outdir"
    trap 'rm -rf "$tmpd"' RETURN
    printf "(A:%.6f,B:%.6f);\n" "$bl" "$bl" > "$tmpd/tree2.nwk"
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

run_fkmv_size_ext() {
    local outdir="$1" model_name="$2" bl="$3"
    "$FKMV_EXE" --dir "$outdir" --model "$model_name" --threads "$THREADS" --branch-len "$bl" \
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

[[ -f "$OUT_CSV" ]] || echo "dir,model,branch_len,algo,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,n,mean_exact_J,mean_est_J,MAE,RMSE" > "$OUT_CSV"

n_total=$(awk -F, -v m="$MODEL_FILTER" 'NR>1 && $3==m {c++} END{print c+0}' "$DATASET_CSV")
echo "JC69 FastKMV size ext: k=21, sizes 8k..8M (x2), pairs=$FORCE_PAIRS"
echo "THREADS=$THREADS  scenarios(csv rows)=$n_total  CSV=$OUT_CSV"

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

    echo "[$n] RUN $dataset  n=$FORCE_PAIRS"
    mkdir -p "$outdir"
    iqt="${MODEL_STR[$model]:-JC}"

    cleanup_scenario_data() {
        rm -f "$outdir"/pair_*.fa "$outdir"/pairs_meta.tsv 2>/dev/null || true
    }
    trap cleanup_scenario_data RETURN

    if [[ "$len_mode" == "fixed" ]]; then
        gen_alisim_batch_fixed "$model" "$iqt" "$nominal_len" "$branch_len" "$outdir" "$FORCE_PAIRS"
    else
        gen_alisim_batch_varlen "$model" "$iqt" "$nominal_len" "$branch_len" "$outdir" "$FORCE_PAIRS"
    fi

    run_fkmv_size_ext "$outdir" "$model" "$branch_len"
    cleanup_scenario_data
    trap - RETURN
    touch "$marker"
done < "$DATASET_CSV"

echo "Done at $(date -Iseconds)"
