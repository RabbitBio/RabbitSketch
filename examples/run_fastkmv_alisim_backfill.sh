#!/usr/bin/env bash
# Re-run FastKMV on the full AliSim grid (9 lengths × models × branch lengths).
# Regenerates FASTA per combo (deleted after original sweep), evaluates 18
# FastKMV parameter points, appends to alisim_sweep.csv.
#
# Usage:
#   bash run_fastkmv_alisim_backfill.sh [THREADS] [PAIRS_MAX] [MAX_COMBOS]
#   MAX_COMBOS=0 → all rows in alisim_dataset_summary.csv (default)
#
#   nohup bash run_fastkmv_alisim_backfill.sh 64 1000 > results/fastkmv_alisim_backfill.log 2>&1 &
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
MAX_COMBOS="${3:-0}"

OUTBASE="/home/bigssd/zt/ref"
RESDIR="results"
SWEEP_CSV="$RESDIR/alisim_sweep.csv"
DATASET_CSV="$RESDIR/alisim_dataset_summary.csv"
FKMV_EXE="./exe_eval_fastkmv_alisim_supplement"
MARKER_NAME=".fastkmv_lex_redesign_done"
LOG="${LOG:-$RESDIR/fastkmv_alisim_backfill.log}"

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

mkdir -p "$RESDIR" "$OUTBASE"
make -s exe_eval_fastkmv_alisim_supplement
[[ -x "$FKMV_EXE" ]] || { echo "[ERROR] $FKMV_EXE missing"; exit 1; }

IQTREE=""
for cand in iqtree2 iqtree3 iqtree; do
    if p="$(command -v "$cand" 2>/dev/null)" && [[ -n "$p" && -x "$p" ]]; then IQTREE="$p"; break; fi
done
[[ -n "$IQTREE" && -x "$IQTREE" ]] || { echo "[ERROR] iqtree not found"; exit 1; }
echo "Using IQ-TREE: $IQTREE"

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
    local tmpd="/tmp/alisim_fkmv_$$_${model_name}_${len}_${bl//./p}"
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
    local tmpd="/tmp/alisim_fkmv_var_$$_${model_name}_${nominal_len}_${bl//./p}"
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
if hi <= lo: print(lo)
else: print(int(math.exp(random.uniform(math.log(lo), math.log(hi)))))")
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

run_fkmv_sweep() {
    local outdir="$1" model_name="$2" bl="$3"
    "$FKMV_EXE" --dir "$outdir" --model "$model_name" --threads "$THREADS" --branch-len "$bl" \
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
            printf "%s,%s,%.6f,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
                dir,model,bl,algo,kmer,size,bits,ss_a,ss_base,
                bd_bbits,bd_s64,drlevel,n,mej,mest,mae,rmse >> csv
        }'
}

[[ -f "$SWEEP_CSV" ]] || echo "dir,model,branch_len,algo,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,n,mean_exact_J,mean_est_J,MAE,RMSE" > "$SWEEP_CSV"

if [[ ! -f "$RESDIR/.fastkmv_strip_done" ]]; then
    echo "Stripping old fastkmv rows from $SWEEP_CSV ..."
    python3 scripts/strip_algo_rows.py --algo fastkmv "$SWEEP_CSV"
    touch "$RESDIR/.fastkmv_strip_done"
fi

n_total=$(tail -n +2 "$DATASET_CSV" | wc -l)
echo "FastKMV AliSim backfill: up to $n_total combos from $DATASET_CSV"
echo "THREADS=$THREADS  PAIRS_MAX=$PAIRS_MAX  log=$LOG"

n=0
while IFS=, read -r dataset len_mode model nominal_len branch_len len_tag bl_tag pairs _rest; do
    [[ "$dataset" == "dataset" ]] && continue
    [[ -z "$dataset" ]] && continue

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

    n_pairs="$pairs"
    if [[ -z "$n_pairs" || "$n_pairs" == "0" ]]; then
        n_pairs=$(get_pairs_for_len "$nominal_len")
    fi

    echo "[$n] RUN $dataset  n=$n_pairs"
    mkdir -p "$outdir"
    iqt="${MODEL_STR[$model]:-}"
    [[ -n "$iqt" ]] || { echo "  [WARN] unknown model $model"; continue; }

    if [[ "$len_mode" == "fixed" ]]; then
        gen_alisim_batch_fixed "$model" "$iqt" "${MODEL_INDEL[$model]:-}" \
            "$nominal_len" "$branch_len" "$outdir" "$n_pairs"
    else
        gen_alisim_batch_varlen "$model" "$iqt" "$nominal_len" "$branch_len" "$outdir" "$n_pairs"
    fi

    run_fkmv_sweep "$outdir" "$model" "$branch_len"
    rm -f "$outdir"/pair_*.fa "$outdir"/pairs_meta.tsv
    touch "$marker"
done < "$DATASET_CSV"

echo "AliSim FastKMV backfill pass done at $(date -Iseconds)"
echo "Next: python3 build_master_csv.py && replot figures"
