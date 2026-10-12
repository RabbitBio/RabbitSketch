#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_alisim_eval.sh
#
# Two evaluation modes per combo (model × length × branch_len):
#
#   FIXED  – Both sequences are exactly L bp (standard Jaccard scenario).
#             Covers all 11 evolutionary models.
#
#   VARLEN – Within each batch each pair's length is drawn log-uniformly from
#             [L/4, L×4], simulating heterogeneous-length inputs (all sequences
#             for A and B can differ in length across pairs).
#             Run on 3 representative substitution models × 3 BL values.
#
# Lengths: 1k 10k 50k 100k 500k 1M 2M 4M 8M
#
# Space management: FASTA files are deleted after each combo; .done markers
# prevent re-running completed combos on resume.
#
# Output:
#   results/alisim_sweep.csv         – per-(combo,algo,params) accuracy rows
#   results/alisim_dataset_summary.csv
#
# Usage:
#   bash run_alisim_eval.sh [THREADS] [PAIRS_MAX]
#   PAIRS_MAX is the pair count for shortest sequences; scaled down for longer.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"
# shellcheck source=alisim_bl_config.sh
source "$(dirname "$0")/alisim_bl_config.sh"

# Supplement mode: only new Mash-distance branch lengths (see run_mash_dist_supplement.sh)
if [[ "${MASH_SUPPLEMENT_ONLY:-0}" == "1" ]]; then
    BRANCH_LENS=("${SUPPLEMENT_BRANCH_LENS[@]}")
    VARLEN_BLS=("${VARLEN_BLS_SUPPLEMENT[@]}")
fi
# Mash calibration: M≈0.075, 0.150, 0.200 (all algos; see run_mash_calib_supplement.sh)
if [[ "${MASH_CALIB_SUPPLEMENT:-0}" == "1" ]]; then
    BRANCH_LENS=("${MASH_CALIB_BRANCH_LENS[@]}")
    VARLEN_BLS=("${MASH_CALIB_BRANCH_LENS[@]}")
fi

THREADS="${1:-16}"
PAIRS_MAX="${2:-1000}"

# Adaptive pair count: scale down for large L so wall-time stays manageable.
get_pairs_for_len() {
    local len="$1"
    if   (( len <=   1000 )); then echo $(( PAIRS_MAX * 2 ))   # 2000
    elif (( len <=  10000 )); then echo "$PAIRS_MAX"            # 1000
    elif (( len <=  50000 )); then echo $(( PAIRS_MAX / 2  ))  # 500
    elif (( len <= 100000 )); then echo $(( PAIRS_MAX / 4  ))  # 250
    elif (( len <= 500000 )); then echo $(( PAIRS_MAX / 5  ))  # 200
    elif (( len <= 1000000)); then echo $(( PAIRS_MAX / 10 ))  # 100
    elif (( len <= 2000000)); then echo $(( PAIRS_MAX / 20 ))  # 50
    elif (( len <= 4000000)); then echo $(( PAIRS_MAX / 20 ))  # 50
    else                           echo $(( PAIRS_MAX / 50 ))  # 20
    fi
}

OUTBASE="/home/bigssd/zt/ref"
RESDIR="results"
SWEEP_CSV="$RESDIR/alisim_sweep.csv"
DATASET_CSV="$RESDIR/alisim_dataset_summary.csv"
SWEEP_EXE="./exe_eval_alisim_sweep"

mkdir -p "$RESDIR" "$OUTBASE"

if [[ ! -x "$SWEEP_EXE" ]]; then
    echo "[ERROR] $SWEEP_EXE not found. Run: make exe_eval_alisim_sweep"; exit 1
fi

# ── locate iqtree ─────────────────────────────────────────────────────────────
IQTREE=""
for cand in iqtree2 iqtree3 iqtree; do
    if p="$(command -v "$cand" 2>/dev/null)" && [[ -n "$p" && -x "$p" ]]; then
        IQTREE="$p"
        break
    fi
done
if [[ -z "$IQTREE" ]]; then
    CONDA_BIN=""
    if command -v conda &>/dev/null; then
        CONDA_BIN="$(conda info --base 2>/dev/null)/envs/cpp20/bin"
    fi
    for d in "$CONDA_BIN" \
             "${HOME}/anaconda3/envs/cpp20/bin" \
             "/home/user_home/tongz/anaconda3/envs/cpp20/bin"; do
        [[ -z "$d" || ! -d "$d" ]] && continue
        for cand in iqtree2 iqtree3 iqtree; do
            if [[ -x "$d/$cand" ]]; then IQTREE="$d/$cand"; break 2; fi
        done
    done
fi
if [[ -z "$IQTREE" ]] || [[ ! -x "$IQTREE" ]]; then
    echo "[ERROR] iqtree not found."; exit 1
fi
echo "Using IQ-TREE: $IQTREE"

# ── evolutionary models ───────────────────────────────────────────────────────
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

MODEL_INDEL["JC69"]=""
MODEL_INDEL["HKY_ts2"]=""
MODEL_INDEL["HKY_ts5"]=""
MODEL_INDEL["GTR_ATrich"]=""
MODEL_INDEL["GTR_G02"]=""
MODEL_INDEL["GTR_G05"]=""
MODEL_INDEL["GTR_G20"]=""
MODEL_INDEL["GTR_highGC"]=""
MODEL_INDEL["GTR_G_indel1"]="0.005,0.005"
MODEL_INDEL["GTR_G_indel2"]="0.01,0.01"
MODEL_INDEL["GTR_G_indel3"]="0.03,0.03"

# Deterministic order
MODEL_ORDER=(JC69 HKY_ts2 HKY_ts5 GTR_ATrich GTR_G02 GTR_G05 GTR_G20 GTR_highGC
             GTR_G_indel1 GTR_G_indel2 GTR_G_indel3)

# 3 representative models for variable-length experiment
VARLEN_MODELS=(JC69 GTR_G05 GTR_ATrich)

# ── lengths ───────────────────────────────────────────────────────────────────
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

# BRANCH_LENS, BL_TAG, VARLEN_BLS from alisim_bl_config.sh (may be overridden above)

# ── CSV headers ───────────────────────────────────────────────────────────────
if [[ ! -f "$SWEEP_CSV" ]]; then
    echo "dir,model,branch_len,algo,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,n,mean_exact_J,mean_est_J,MAE,RMSE" \
        > "$SWEEP_CSV"
fi
if [[ ! -f "$DATASET_CSV" ]]; then
    echo "dataset,len_mode,model,nominal_len,branch_len,len_tag,bl_tag,pairs,files,total_bp,generated_at,deleted_at" \
        > "$DATASET_CSV"
fi

# ── generate one FIXED-length alisim batch ────────────────────────────────────
# Each pair is exactly $len bp.
gen_alisim_batch_fixed() {
    local model_name="$1" iqt_model="$2" indel_rates="$3"
    local len="$4" bl="$5" outdir="$6" n_pairs="$7"
    local tmpd="/tmp/alisim_$$_${model_name}_${len}_${bl//./p}"
    mkdir -p "$tmpd" "$outdir"
    printf "(A:%.6f,B:%.6f);\n" "$bl" "$bl" > "$tmpd/tree2.nwk"

    local pids=() slot=0
    for ((i=0; i<n_pairs; i++)); do
        local cmd=("$IQTREE" --alisim "${tmpd}/p${i}_sim"
            -m "$iqt_model" -t "$tmpd/tree2.nwk"
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

# ── generate one VARIABLE-length alisim batch ─────────────────────────────────
# Each pair i has its own length sampled log-uniformly from
#   [max(200, len/4), min(len*4, 8M)]
# (capped at 8M to keep single-call time bounded; minimum 200 bp).
# A and B within the same pair share that pair's length; lengths differ
# across pairs in the batch.
gen_alisim_batch_varlen() {
    local model_name="$1" iqt_model="$2"
    local nominal_len="$3" bl="$4" outdir="$5" n_pairs="$6"
    local tmpd="/tmp/alisim_$$_varlen_${model_name}_${nominal_len}_${bl//./p}"
    mkdir -p "$tmpd" "$outdir"
    printf "(A:%.6f,B:%.6f);\n" "$bl" "$bl" > "$tmpd/tree2.nwk"

    local pids=() slot=0
    for ((i=0; i<n_pairs; i++)); do
        # log-uniform sample in [max(200, L/4), min(L*4, 8M)]
        local pair_len
        pair_len=$(python3 -c "
import random, math
random.seed($i * 7919 + 31337)
lo = max(200, $nominal_len // 4)
hi = min($nominal_len * 4, 8_000_000)
if hi <= lo: print(lo)
else: print(int(math.exp(random.uniform(math.log(lo), math.log(hi)))))")
        local cmd=("$IQTREE" --alisim "${tmpd}/p${i}_sim"
            -m "$iqt_model" -t "$tmpd/tree2.nwk"
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

# ── helper: split multi-fasta → pair_N_A/B.fa + pairs_meta.tsv ───────────────
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
            END { printf ">A\n%s\n", seqA > oa
                  printf ">B\n%s\n", seqB > ob }
        ' "$fa"
        local la lb
        la=$(grep -v "^>" "$outdir/pair_${i}_A.fa" | tr -d '\n' | wc -c)
        lb=$(grep -v "^>" "$outdir/pair_${i}_B.fa" | tr -d '\n' | wc -c)
        printf "%d\t%.6f\t%d\t%d\n" "$i" "$bl" "$la" "$lb" >> "$meta"
    done
}

# ── record dataset metadata ───────────────────────────────────────────────────
record_dataset_summary() {
    local outdir="$1" len_mode="$2" model_name="$3"
    local nominal_len="$4" bl="$5" generated_at="$6" n_pairs="$7"
    local lt="${LEN_TAG[$nominal_len]}" blt="${BL_TAG[$bl]}"
    local dataset="${model_name}/${len_mode}_${lt}/${blt}"
    local files total_bp
    files=$(ls "$outdir"/pair_*.fa 2>/dev/null | wc -l)
    total_bp=$(awk 'BEGIN{s=0} /^>/ {next} {s+=length($0)} END{print s+0}' \
               "$outdir"/pair_*.fa 2>/dev/null || echo 0)
    printf '%s,%s,%s,%d,%.6f,%s,%s,%d,%d,%d,%s,\n' \
        "$dataset" "$len_mode" "$model_name" "$nominal_len" "$bl" \
        "$lt" "$blt" "$n_pairs" "$files" "$total_bp" \
        "$generated_at" >> "$DATASET_CSV"
}

mark_dataset_deleted() {
    local dataset="$1" deleted_at="$2"
    python3 -c "
import csv
path='$DATASET_CSV'; dataset='$dataset'; deleted_at='$deleted_at'
rows=[]
with open(path,newline='') as f: rows=list(csv.reader(f))
hdr=rows[0]; di=hdr.index('dataset'); del_i=hdr.index('deleted_at')
for row in rows[1:]:
    if row[di]==dataset and (len(row)<=del_i or row[del_i]==''):
        while len(row)<=del_i: row.append('')
        row[del_i]=deleted_at
with open(path,'w',newline='') as f: csv.writer(f).writerows(rows)
" 2>/dev/null || true
}

# ── run full parameter sweep on one combo directory ───────────────────────────
run_full_param_sweep() {
    local outdir="$1" model_name="$2"
    "$SWEEP_EXE" --dir "$outdir" --model "$model_name" --threads "$THREADS" \
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
                dir,model,bl,algo,kmer,size,bits,
                ss_a,ss_base,bd_bbits,bd_s64,drlevel,
                n,mej,mest,mae,rmse >> csv
        }' 2>/dev/null
}

# ── per-combo driver ──────────────────────────────────────────────────────────
# Sets global _COMBO_STATUS="skip"|"done"
process_combo() {
    local len_mode="$1" model_name="$2" iqt_model="$3" indel_rates="$4"
    local nominal_len="$5" bl="$6"
    local lt="${LEN_TAG[$nominal_len]}" blt="${BL_TAG[$bl]}"
    local outdir="$OUTBASE/${model_name}/${len_mode}_${lt}/${blt}"
    local dataset="${model_name}/${len_mode}_${lt}/${blt}"

    if [[ -f "$outdir/.done" ]]; then _COMBO_STATUS="skip"; return 0; fi

    local n_pairs generated_at
    n_pairs=$(get_pairs_for_len "$nominal_len")
    printf "  %-16s %-10s %-8s [n=%4d]  gen..." \
           "$model_name" "${len_mode}_${lt}" "$blt" "$n_pairs"
    generated_at="$(date -Iseconds)"

    if [[ "$len_mode" == "fixed" ]]; then
        gen_alisim_batch_fixed "$model_name" "$iqt_model" "$indel_rates" \
                               "$nominal_len" "$bl" "$outdir" "$n_pairs"
    else
        gen_alisim_batch_varlen "$model_name" "$iqt_model" \
                                "$nominal_len" "$bl" "$outdir" "$n_pairs"
    fi

    record_dataset_summary "$outdir" "$len_mode" "$model_name" \
                           "$nominal_len" "$bl" "$generated_at" "$n_pairs"
    printf "sweep..."
    run_full_param_sweep "$outdir" "$model_name"
    printf "clean..."
    rm -f "$outdir"/pair_*.fa
    mark_dataset_deleted "$dataset" "$(date -Iseconds)"
    touch "$outdir/.done"
    printf "done\n"
    _COMBO_STATUS="done"
}

# ── main ──────────────────────────────────────────────────────────────────────
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║  AliSim full-parameter evaluation                               ║"
echo "╠══════════════════════════════════════════════════════════════════╣"
printf "║  THREADS=%-4s  PAIRS_MAX=%-6s  Start: %-22s ║\n" \
       "$THREADS" "$PAIRS_MAX" "$(date '+%Y-%m-%d %H:%M:%S')"
echo "║  Lengths: 1k 10k 50k 100k 500k 1M 2M 4M 8M                     ║"
printf "║  Fixed-len : 11 models × 9 lengths × %2d BLs                  ║\n" "${#BRANCH_LENS[@]}"
printf "║  Varlen    :  3 models × 9 lengths × %2d BLs                  ║\n" "${#VARLEN_BLS[@]}"
if [[ "${MASH_SUPPLEMENT_ONLY:-0}" == "1" ]]; then
    echo "║  (Mash-distance supplement — new BLs only)                     ║"
fi
echo "║  Varlen range: [max(200,L/4), min(L*4, 8M)]                    ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""

DONE=0; SKIP=0
_COMBO_STATUS=""

# ── Phase 1: Fixed-length (all models) ────────────────────────────────────────
echo "▶ Phase 1: Fixed-length combos"
for model_name in "${MODEL_ORDER[@]}"; do
    iqt_model="${MODEL_STR[$model_name]}"
    indel_rates="${MODEL_INDEL[$model_name]}"
    for len in "${LENGTHS[@]}"; do
        for bl in "${BRANCH_LENS[@]}"; do
            _COMBO_STATUS=""
            process_combo "fixed" "$model_name" "$iqt_model" "$indel_rates" \
                          "$len" "$bl"
            [[ "$_COMBO_STATUS" == "skip" ]] && SKIP=$((SKIP+1)) || DONE=$((DONE+1))
        done
    done
done

# ── Phase 2: Variable-length (representative models) ──────────────────────────
echo ""
echo "▶ Phase 2: Variable-length combos"
for model_name in "${VARLEN_MODELS[@]}"; do
    iqt_model="${MODEL_STR[$model_name]}"
    for len in "${LENGTHS[@]}"; do
        for bl in "${VARLEN_BLS[@]}"; do
            _COMBO_STATUS=""
            process_combo "varlen" "$model_name" "${iqt_model}" "" \
                          "$len" "$bl"
            [[ "$_COMBO_STATUS" == "skip" ]] && SKIP=$((SKIP+1)) || DONE=$((DONE+1))
        done
    done
done

echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
printf "║  Done: %-4d combos processed, %-4d skipped                      ║\n" \
       "$DONE" "$SKIP"
printf "║  End:  %-43s   ║\n" "$(date '+%Y-%m-%d %H:%M:%S')"
printf "║  CSV:  %-43s   ║\n" "$SWEEP_CSV"
echo "╚══════════════════════════════════════════════════════════════════╝"
