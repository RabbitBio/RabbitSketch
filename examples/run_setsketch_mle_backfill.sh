#!/usr/bin/env bash
# Backfill setsketch_mle rows for the original 8 branch lengths (fixed + varlen).
# FASTA are regenerated (deleted after main eval); appends to results/alisim_sweep.csv.
#
# Convention: SetSketch in all figures = joint-MLE (algo setsketch_mle → setsketch in master).
#
# Usage:
#   bash run_setsketch_mle_backfill.sh [THREADS] [PAIRS_MAX]
#   nohup bash run_setsketch_mle_backfill.sh 64 1000 > results/setsketch_mle_backfill.log 2>&1 &
#
# After completion:
#   python3 build_master_csv.py
#   MPLCONFIGDIR=/tmp/mpl-cache-rabbitsketch python3 plot_accuracy_landscape.py
#   MPLCONFIGDIR=/tmp/mpl-cache-rabbitsketch python3 plot_accuracy_landscape_oracle.py
set -euo pipefail
cd "$(dirname "$0")"
# shellcheck source=alisim_bl_config.sh
source "$(dirname "$0")/alisim_bl_config.sh"

for _iqbin in "${HOME}/anaconda3/envs/cpp20/bin" \
               "/home/user_home/tongz/anaconda3/envs/cpp20/bin"; do
    [[ -d "$_iqbin" ]] && export PATH="$_iqbin:$PATH" && break
done

THREADS="${1:-64}"
PAIRS_MAX="${2:-1000}"
LOG="${LOG:-results/setsketch_mle_backfill.log}"

OUTBASE="/home/bigssd/zt/ref"
RESDIR="results"
SWEEP_CSV="$RESDIR/alisim_sweep.csv"
SWEEP_EXE="./exe_eval_alisim_sweep"
MARKER_NAME=".mle_backfill_done"

mkdir -p "$RESDIR" "$OUTBASE"
[[ -x "$SWEEP_EXE" ]] || { echo "[ERROR] $SWEEP_EXE missing"; exit 1; }

IQTREE=""
for cand in iqtree2 iqtree3 iqtree; do
    if p="$(command -v "$cand" 2>/dev/null)" && [[ -n "$p" && -x "$p" ]]; then IQTREE="$p"; break; fi
done
[[ -n "$IQTREE" && -x "$IQTREE" ]] || { echo "[ERROR] iqtree not found"; exit 1; }
echo "Using IQ-TREE: $IQTREE"

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

MODEL_ORDER=(JC69 HKY_ts2 HKY_ts5 GTR_ATrich GTR_G02 GTR_G05 GTR_G20 GTR_highGC
             GTR_G_indel1 GTR_G_indel2 GTR_G_indel3)
VARLEN_MODELS=(JC69 GTR_G05 GTR_ATrich)
LENGTHS=(1000 10000 50000 100000 500000 1000000 2000000 4000000 8000000)
declare -A LEN_TAG
LEN_TAG[1000]="L1k"
LEN_TAG[10000]="L10k"
LEN_TAG[50000]="L50k"
LEN_TAG[100000]="L100k"
LEN_TAG[500000]="L500k"
LEN_TAG[1000000]="L1M"
LEN_TAG[2000000]="L2M"
LEN_TAG[4000000]="L4M"
LEN_TAG[8000000]="L8M"

get_pairs_for_len() {
    local len="$1"
    if   (( len <=   1000 )); then echo $(( PAIRS_MAX * 2 ))
    elif (( len <=  10000 )); then echo "$PAIRS_MAX"
    elif (( len <=  50000 )); then echo $(( PAIRS_MAX / 2  ))
    elif (( len <= 100000 )); then echo $(( PAIRS_MAX / 4  ))
    elif (( len <= 500000 )); then echo $(( PAIRS_MAX / 5  ))
    elif (( len <= 1000000)); then echo $(( PAIRS_MAX / 10 ))
    elif (( len <= 2000000)); then echo $(( PAIRS_MAX / 20 ))
    elif (( len <= 4000000)); then echo $(( PAIRS_MAX / 20 ))
    else                           echo $(( PAIRS_MAX / 50 ))
    fi
}

_split_and_meta() {
    local tmpd="$1" outdir="$2" bl="$3" n_pairs="$4"
    local meta="$outdir/pairs_meta.tsv"
    printf "pair_id\tbranch_len\tlen_A\tlen_B\n" > "$meta"
    for ((i=0; i<n_pairs; i++)); do
        local fa="${tmpd}/p${i}_sim.fa"
        [[ -f "$fa" ]] || continue
        awk -v oa="$outdir/pair_${i}_A.fa" -v ob="$outdir/pair_${i}_B.fa" '
            /^>/ { rec++; next }
            rec==1 { seqA=seqA $0 }
            rec==2 { seqB=seqB $0 }
            END { printf ">A\n%s\n", seqA > oa; printf ">B\n%s\n", seqB > ob }
        ' "$fa"
        local la lb
        la=$(grep -v "^>" "$outdir/pair_${i}_A.fa" | tr -d '\n' | wc -c)
        lb=$(grep -v "^>" "$outdir/pair_${i}_B.fa" | tr -d '\n' | wc -c)
        printf "%d\t%.6f\t%d\t%d\n" "$i" "$bl" "$la" "$lb" >> "$meta"
    done
}

gen_alisim_batch_fixed() {
    local model_name="$1" iqt_model="$2" indel_rates="$3"
    local len="$4" bl="$5" outdir="$6" n_pairs="$7"
    local tmpd="/tmp/alisim_mlebf_$$_${model_name}_${len}_${bl//./p}"
    mkdir -p "$tmpd" "$outdir"
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
    _split_and_meta "$tmpd" "$outdir" "$bl" "$n_pairs"
    rm -rf "$tmpd"
}

gen_alisim_batch_varlen() {
    local model_name="$1" iqt_model="$2"
    local nominal_len="$3" bl="$4" outdir="$5" n_pairs="$6"
    local tmpd="/tmp/alisim_mlebf_vl_$$_${model_name}_${nominal_len}_${bl//./p}"
    mkdir -p "$tmpd" "$outdir"
    printf "(A:%.6f,B:%.6f);\n" "$bl" "$bl" > "$tmpd/tree2.nwk"
    local pids=() slot=0
    for ((i=0; i<n_pairs; i++)); do
        local pair_len
        pair_len=$(python3 -c "
import random, math
random.seed($i * 7919 + 31337)
lo = max(200, $nominal_len // 4)
hi = min($nominal_len * 4, 8_000_000)
print(lo if hi <= lo else int(math.exp(random.uniform(math.log(lo), math.log(hi)))))")
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
    _split_and_meta "$tmpd" "$outdir" "$bl" "$n_pairs"
    rm -rf "$tmpd"
}

run_mle_sweep() {
    local outdir="$1" model_name="$2"
    "$SWEEP_EXE" --dir "$outdir" --model "$model_name" --threads "$THREADS" \
                 --algos setsketch_mle \
      | awk -v csv="$SWEEP_CSV" '
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
            printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
                   dir,model,bl,algo,kmer,size,bits,ss_a,ss_base,
                   bd_bbits,bd_s64,drlevel,n,mej,mest,mae,rmse >> csv
        }'
}

process_combo() {
    local len_mode="$1" model_name="$2" iqt_model="$3" indel_rates="$4"
    local nominal_len="$5" bl="$6"
    local lt="${LEN_TAG[$nominal_len]}" blt="${BL_TAG[$bl]}"
    local outdir="$OUTBASE/${model_name}/${len_mode}_${lt}/${blt}"
    local marker="$outdir/$MARKER_NAME"

    if [[ -f "$marker" ]]; then
        echo "  SKIP ${model_name}/${len_mode}_${lt}/${blt} (MLE backfill done)"
        return 0
    fi

    local n_pairs
    n_pairs=$(get_pairs_for_len "$nominal_len")
    echo "  RUN  ${model_name}/${len_mode}_${lt}/${blt}  n=$n_pairs"
    if [[ "$len_mode" == "fixed" ]]; then
        gen_alisim_batch_fixed "$model_name" "$iqt_model" "$indel_rates" \
            "$nominal_len" "$bl" "$outdir" "$n_pairs"
    else
        gen_alisim_batch_varlen "$model_name" "$iqt_model" \
            "$nominal_len" "$bl" "$outdir" "$n_pairs"
    fi
    run_mle_sweep "$outdir" "$model_name"
    rm -f "$outdir"/pair_*.fa
    touch "$marker"
}

[[ -f "$SWEEP_CSV" ]] || echo "dir,model,branch_len,algo,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,n,mean_exact_J,mean_est_J,MAE,RMSE" > "$SWEEP_CSV"

n_fixed=$(( ${#MODEL_ORDER[@]} * ${#LENGTHS[@]} * ${#ORIGINAL_BRANCH_LENS[@]} ))
n_varlen=$(( ${#VARLEN_MODELS[@]} * ${#LENGTHS[@]} * ${#VARLEN_BLS[@]} ))
echo "SetSketch MLE backfill: original BLs fixed=$n_fixed varlen=$n_varlen"
echo "THREADS=$THREADS  PAIRS_MAX=$PAIRS_MAX  log=$LOG"

n=0
for model in "${MODEL_ORDER[@]}"; do
    for L in "${LENGTHS[@]}"; do
        for bl in "${ORIGINAL_BRANCH_LENS[@]}"; do
            n=$((n + 1))
            echo "[$n/$n_fixed] fixed"
            process_combo fixed "$model" "${MODEL_STR[$model]}" "${MODEL_INDEL[$model]:-}" \
                "$L" "$bl"
        done
    done
done

for model in "${VARLEN_MODELS[@]}"; do
    for L in "${LENGTHS[@]}"; do
        for bl in "${VARLEN_BLS[@]}"; do
            echo "varlen $model L=$L bl=$bl"
            process_combo varlen "$model" "${MODEL_STR[$model]}" "" "$L" "$bl"
        done
    done
done

echo "Done. Rebuild: python3 build_master_csv.py && replot landscapes."
