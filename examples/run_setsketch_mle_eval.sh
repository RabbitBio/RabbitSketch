#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_setsketch_mle_eval.sh
#
# Targeted re-run for SetSketch with both inclusion-exclusion and joint-MLE
# (Ertl 2021) estimators. Produces head-to-head data for two representative
# evolutionary models, all 9 fixed lengths × 8 branch lengths × full SetSketch
# parameter grid.
#
# Output:  results/alisim_setsketch_mle.csv   (same schema as alisim_sweep.csv)
#
# Usage:   bash run_setsketch_mle_eval.sh [THREADS] [PAIRS_MAX]
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"
# shellcheck source=alisim_bl_config.sh
source "$(dirname "$0")/alisim_bl_config.sh"

THREADS="${1:-32}"
PAIRS_MAX="${2:-300}"

# Adaptive pair count (smaller than full alisim run – we only need stable means).
get_pairs_for_len() {
    local len="$1"
    if   (( len <=   1000 )); then echo $(( PAIRS_MAX * 2 ))
    elif (( len <=  10000 )); then echo "$PAIRS_MAX"
    elif (( len <=  50000 )); then echo $(( PAIRS_MAX / 2  ))
    elif (( len <= 100000 )); then echo $(( PAIRS_MAX / 3  ))
    elif (( len <= 500000 )); then echo $(( PAIRS_MAX / 5  ))
    elif (( len <= 1000000)); then echo $(( PAIRS_MAX / 10 ))
    elif (( len <= 2000000)); then echo $(( PAIRS_MAX / 20 ))
    elif (( len <= 4000000)); then echo $(( PAIRS_MAX / 30 ))
    else                           echo $(( PAIRS_MAX / 50 ))
    fi
}

OUTBASE="/home/bigssd/zt/ref_setsketch_mle"
RESDIR="results"
SWEEP_CSV="$RESDIR/alisim_setsketch_mle.csv"
SWEEP_EXE="./exe_eval_alisim_sweep"

mkdir -p "$RESDIR" "$OUTBASE"
[[ -x "$SWEEP_EXE" ]] || { echo "[ERROR] $SWEEP_EXE not found"; exit 1; }

# locate iqtree
IQTREE=""
for cand in iqtree2 iqtree3 iqtree; do
    if p="$(command -v "$cand" 2>/dev/null)" && [[ -n "$p" && -x "$p" ]]; then
        IQTREE="$p"
        break
    fi
done
if [[ -z "$IQTREE" ]]; then
    for d in "$(conda info --base 2>/dev/null)/envs/cpp20/bin" \
             "${HOME}/anaconda3/envs/cpp20/bin" \
             "/home/user_home/tongz/anaconda3/envs/cpp20/bin"; do
        [[ -z "$d" || "$d" == "/envs/cpp20/bin" || ! -d "$d" ]] && continue
        for cand in iqtree2 iqtree3 iqtree; do
            if [[ -x "$d/$cand" ]]; then IQTREE="$d/$cand"; break 2; fi
        done
    done
fi
[[ -z "$IQTREE" || ! -x "$IQTREE" ]] && { echo "[ERROR] iqtree not found."; exit 1; }
echo "Using IQ-TREE: $IQTREE"

# ── two representative models ────────────────────────────────────────────────
declare -A MODEL_STR
MODEL_STR["JC69"]="JC"
MODEL_STR["GTR_G05"]="GTR+G{0.5}"
MODEL_ORDER=(JC69 GTR_G05)

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

# BRANCH_LENS, BL_TAG from alisim_bl_config.sh

# CSV header (same schema as alisim_sweep.csv)
if [[ ! -f "$SWEEP_CSV" ]]; then
    echo "dir,model,branch_len,algo,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,n,mean_exact_J,mean_est_J,MAE,RMSE" \
        > "$SWEEP_CSV"
fi

# ── alisim batch generation (fixed-length only) ──────────────────────────────
gen_alisim_batch_fixed() {
    local model_name="$1" iqt_model="$2" len="$3" bl="$4" outdir="$5" n_pairs="$6"
    local tmpd="/tmp/alisim_ssmle_$$_${model_name}_${len}_${bl//./p}"
    mkdir -p "$tmpd" "$outdir"
    printf "(A:%.6f,B:%.6f);\n" "$bl" "$bl" > "$tmpd/tree2.nwk"
    local pids=() slot=0
    for ((i=0; i<n_pairs; i++)); do
        "$IQTREE" --alisim "${tmpd}/p${i}_sim" -m "$iqt_model" -t "$tmpd/tree2.nwk" \
                  --length "$len" --out-format fasta --seed $((42+i)) -nt 1 -quiet \
                  >/dev/null 2>&1 &
        pids+=($!)
        slot=$((slot + 1))
        if (( slot >= THREADS )); then
            wait "${pids[0]}" || true; pids=("${pids[@]:1}"); slot=$((slot - 1))
        fi
    done
    for p in "${pids[@]}"; do wait "$p" || true; done

    local meta="$outdir/pairs_meta.tsv"
    printf "pair_id\tbranch_len\tlen_A\tlen_B\n" > "$meta"
    for ((i=0; i<n_pairs; i++)); do
        local fa="${tmpd}/p${i}_sim.fa"
        [[ -f "$fa" ]] || continue
        awk -v oa="$outdir/pair_${i}_A.fa" -v ob="$outdir/pair_${i}_B.fa" '
            /^>/ { rec++; next }
            rec==1 { seqA=seqA $0 }
            rec==2 { seqB=seqB $0 }
            END { printf ">A\n%s\n", seqA > oa
                  printf ">B\n%s\n", seqB > ob }
        ' "$fa"
        local la lb
        la=$(grep -v "^>" "$outdir/pair_${i}_A.fa" | tr -d '\n' | wc -c)
        lb=$(grep -v "^>" "$outdir/pair_${i}_B.fa" | tr -d '\n' | wc -c)
        printf "%d\t%.6f\t%d\t%d\n" "$i" "$bl" "$la" "$lb" >> "$meta"
    done
    rm -rf "$tmpd"
}

# ── parameter sweep (only SetSketch + SetSketch_MLE) ─────────────────────────
run_setsketch_sweep() {
    local outdir="$1" model_name="$2"
    "$SWEEP_EXE" --dir "$outdir" --model "$model_name" --threads "$THREADS" \
                 --algos setsketch,setsketch_mle \
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

# ── main loop ────────────────────────────────────────────────────────────────
total_combos=$(( ${#MODEL_ORDER[@]} * ${#LENGTHS[@]} * ${#BRANCH_LENS[@]} ))
echo "Plan: ${#MODEL_ORDER[@]} models × ${#LENGTHS[@]} lengths × ${#BRANCH_LENS[@]} BLs = $total_combos combos"
echo "Threads: $THREADS  PAIRS_MAX: $PAIRS_MAX"
echo "Output: $SWEEP_CSV"
echo "─────────────────────────────────────────────"

n_done=0
for model in "${MODEL_ORDER[@]}"; do
    iqt_model="${MODEL_STR[$model]}"
    for L in "${LENGTHS[@]}"; do
        lt="${LEN_TAG[$L]}"
        for bl in "${BRANCH_LENS[@]}"; do
            blt="${BL_TAG[$bl]}"
            outdir="$OUTBASE/$model/fixed_$lt/$blt"
            done_marker="$outdir/.done"
            n_pairs="$(get_pairs_for_len "$L")"
            n_done=$((n_done + 1))

            if [[ -f "$done_marker" ]]; then
                echo "[$n_done/$total_combos] $model $lt $blt  SKIP (already done)"
                continue
            fi

            echo "[$n_done/$total_combos] $model $lt $blt  n=$n_pairs"
            printf '  gen...'
            t0=$SECONDS
            gen_alisim_batch_fixed "$model" "$iqt_model" "$L" "$bl" "$outdir" "$n_pairs"
            printf 'done(%ds)  sweep...' "$((SECONDS - t0))"
            t0=$SECONDS
            run_setsketch_sweep "$outdir" "$model"
            printf 'done(%ds)  clean...' "$((SECONDS - t0))"
            rm -f "$outdir"/pair_*.fa
            touch "$done_marker"
            echo "done"
        done
    done
done

echo
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║  Done: SetSketch MLE eval finished                              ║"
echo "║  End:  $(date '+%Y-%m-%d %H:%M:%S')                              ║"
echo "║  CSV:  $SWEEP_CSV                              ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
