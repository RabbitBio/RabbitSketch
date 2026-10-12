#!/usr/bin/env bash
# Legacy backfill: SetSketch MLE k sweep at a=5 (exe_eval_alisim_sweep now runs
# the same k/bits grids at a=5 and a=20 by default; use this only for old CSVs).
# Backfill SetSketch MLE k-mer sweep at default a=5, base=2, bits=13.
# Targets k-sensitivity grid: 11 models × 6 lengths × 5 Mash columns (BLs).
# k = 11,13,15,17,19,21,23,25,27,29,31 (matches plot_k_sensitivity_grid.py).
#
# Usage:
#   bash run_setsketch_mle_k_backfill.sh [THREADS] [PAIRS_MAX]
#   nohup bash run_setsketch_mle_k_backfill.sh 64 1000 > results/setsketch_mle_k_backfill.log 2>&1 &
#
# After completion (also run automatically at end of this script):
#   python3 build_master_csv.py
#   MPLCONFIGDIR=/tmp/mpl-cache-rabbitsketch python3 plot_k_sensitivity_grid.py
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
LOG="${LOG:-results/setsketch_mle_k_backfill.log}"

OUTBASE="/home/bigssd/zt/ref"
RESDIR="results"
SWEEP_CSV="$RESDIR/alisim_sweep.csv"
K_EXE="./exe_eval_ss_mle_k_supplement"
MARKER_NAME=".mle_k_a5_done"

SENSITIVITY_LENGTHS=(1000 100000 1000000 2000000 4000000 8000000)
SENSITIVITY_BLS=(0.005 0.025 0.053 0.116 0.31)

mkdir -p "$RESDIR" "$OUTBASE"
make -s exe_eval_ss_mle_k_supplement
[[ -x "$K_EXE" ]] || { echo "[ERROR] $K_EXE missing"; exit 1; }

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

declare -A LEN_TAG
LEN_TAG[1000]="L1k"
LEN_TAG[100000]="L100k"
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
    local tmpd="/tmp/alisim_kbf_$$_${model_name}_${len}_${bl//./p}"
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

append_sweep_lines() {
    awk -v csv="$SWEEP_CSV" '
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

run_k_sweep() {
    local outdir="$1" model_name="$2" bl="$3"
    "$K_EXE" --dir "$outdir" --model "$model_name" --threads "$THREADS" --branch-len "$bl" \
      | append_sweep_lines
}

process_combo() {
    local model_name="$1" iqt_model="$2" indel_rates="$3"
    local nominal_len="$4" bl="$5"
    local lt="${LEN_TAG[$nominal_len]}" blt="${BL_TAG[$bl]}"
    local outdir="$OUTBASE/${model_name}/fixed_${lt}/${blt}"
    local marker="$outdir/$MARKER_NAME"

    if [[ -f "$marker" ]]; then
        echo "  SKIP ${model_name}/fixed_${lt}/${blt} (k a=5 done)"
        return 0
    fi

    local n_pairs
    n_pairs=$(get_pairs_for_len "$nominal_len")
    echo "  RUN  ${model_name}/fixed_${lt}/${blt}  n=$n_pairs"
    if [[ ! -f "$outdir/pair_0_A.fa" ]]; then
        gen_alisim_batch_fixed "$model_name" "$iqt_model" "${MODEL_INDEL[$model_name]:-}" \
            "$nominal_len" "$bl" "$outdir" "$n_pairs"
    fi
    run_k_sweep "$outdir" "$model_name" "$bl"
    rm -f "$outdir"/pair_*.fa
    touch "$marker"
}

[[ -f "$SWEEP_CSV" ]] || echo "dir,model,branch_len,algo,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,n,mean_exact_J,mean_est_J,MAE,RMSE" > "$SWEEP_CSV"

n_total=$(( ${#MODEL_ORDER[@]} * ${#SENSITIVITY_LENGTHS[@]} * ${#SENSITIVITY_BLS[@]} ))
echo "SetSketch MLE k@a=5 backfill: $n_total combos (k-sensitivity grid)"
echo "THREADS=$THREADS  PAIRS_MAX=$PAIRS_MAX  log=$LOG"

n=0
for model in "${MODEL_ORDER[@]}"; do
    for L in "${SENSITIVITY_LENGTHS[@]}"; do
        for bl in "${SENSITIVITY_BLS[@]}"; do
            n=$((n + 1))
            echo "[$n/$n_total]"
            process_combo "$model" "${MODEL_STR[$model]}" "${MODEL_INDEL[$model]:-}" "$L" "$bl"
        done
    done
done

echo "Done. Rebuilding master + k sensitivity plot..."
python3 build_master_csv.py
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mpl-cache-rabbitsketch}"
python3 plot_k_sensitivity_grid.py
echo "All finished."
