#!/usr/bin/env bash
# Backfill SetSketch MLE grid on JC69 (1000 pairs): k/bits sweeps at
# a ∈ {0.5, 1, 5, 10, 15, 20} plus b-sweep at a=20 (requires new exe_eval_alisim_sweep).
#
# Strips existing setsketch_mle rows from the JC69 CSV, re-runs only scenarios
# that are missing the full grid, appends fresh rows.
#
# Usage:
#   bash run_jc69_setsketch_mle_grid_backfill.sh [THREADS]
#   FORCE_ALL=1  re-run every JC69 scenario (not only incomplete)
#
# Example:
#   nohup bash run_jc69_setsketch_mle_grid_backfill.sh 64 \
#     > results/jc69_setsketch_mle_grid_backfill.log 2>&1 &
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${1:-64}"
MODEL="JC69"
FORCE_PAIRS=1000
FORCE_ALL="${FORCE_ALL:-0}"
ALGOS="setsketch_mle"

OUTBASE="/home/bigssd/zt/ref"
RESDIR="results"
DATASET_CSV="$RESDIR/alisim_dataset_summary.csv"
SWEEP_EXE="./exe_eval_alisim_sweep"
OUT_CSV="$RESDIR/alisim_sweep.JC69.pairs1000.selected6.csv"
MARKER_NAME=".pairs1000_selected6_JC69_done"
SS_MARKER=".pairs1000_ss_mle_grid_JC69_done"
LOG="${LOG:-$RESDIR/jc69_setsketch_mle_grid_backfill.log}"

# Expected unique (ss_a, kmer, bits, ss_base) combos per scenario (ss_base=2 only).
EXPECTED_COMBOS=122

for _iqbin in "${HOME}/anaconda3/envs/cpp20/bin" \
               "/home/user_home/tongz/anaconda3/envs/cpp20/bin"; do
    [[ -d "$_iqbin" ]] && export PATH="$_iqbin:$PATH" && break
done

declare -A MODEL_STR MODEL_INDEL
MODEL_STR["JC69"]="JC"
MODEL_INDEL["JC69"]=""

mkdir -p "$RESDIR" "$OUTBASE"
[[ -x "$SWEEP_EXE" ]] || { echo "[ERROR] $SWEEP_EXE missing; run: make exe_eval_alisim_sweep"; exit 1; }

IQTREE=""
for cand in iqtree2 iqtree3 iqtree; do
    if p="$(command -v "$cand" 2>/dev/null)" && [[ -n "$p" && -x "$p" ]]; then IQTREE="$p"; break; fi
done
[[ -n "$IQTREE" && -x "$IQTREE" ]] || { echo "[ERROR] iqtree not found"; exit 1; }

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
    local tmpd="/tmp/alisim_jc69ss_$$_${model_name}_${len}_${bl//./p}"
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
    local tmpd="/tmp/alisim_jc69ss_var_$$_${model_name}_${nominal_len}_${bl//./p}"
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

append_sweep_lines() {
    awk -v csv="$OUT_CSV" '
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

strip_ss_rows_for_dir() {
    local outdir="$1"
    [[ -f "$OUT_CSV" ]] || return 0
    python3 - "$outdir" "$OUT_CSV" <<'PY'
import sys, csv, shutil
from pathlib import Path
outdir, path = sys.argv[1], Path(sys.argv[2])
rows = list(csv.DictReader(path.open()))
kept = [r for r in rows
        if not (r.get("dir") == outdir and r.get("algo") == "setsketch_mle")]
if len(kept) == len(rows):
    sys.exit(0)
bak = path.with_suffix(path.suffix + ".pre_ss_strip.bak")
if not bak.exists():
    shutil.copy2(path, bak)
with path.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(kept)
print(f"stripped {len(rows)-len(kept)} setsketch_mle rows for {outdir}")
PY
}

need_backfill_dir() {
    local outdir="$1"
    if [[ "$FORCE_ALL" == "1" ]]; then return 0; fi
    if [[ -f "$outdir/$SS_MARKER" ]]; then return 1; fi
    if [[ ! -f "$OUT_CSV" ]]; then return 0; fi
    python3 - "$outdir" "$OUT_CSV" "$EXPECTED_COMBOS" <<'PY'
import sys, csv
outdir, csv_path, need = sys.argv[1], sys.argv[2], int(sys.argv[3])
keys = set()
with open(csv_path) as f:
    for row in csv.DictReader(f):
        if row.get("dir") != outdir or row.get("algo") != "setsketch_mle":
            continue
        if float(row.get("ss_base", 2)) != 2.0:
            continue
        keys.add((float(row["ss_a"]), int(row["kmer"]), int(row["bits"])))
sys.exit(0 if len(keys) < need else 1)
PY
}

exec >>"$LOG" 2>&1
echo "=== JC69 SetSketch MLE grid backfill $(date -Iseconds) threads=$THREADS ==="

if [[ ! -f "$OUT_CSV" ]]; then
    echo "dir,model,branch_len,algo,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,n,mean_exact_J,mean_est_J,MAE,RMSE" > "$OUT_CSV"
fi

n_run=0 n_skip=0
while IFS=, read -r dataset len_mode model nominal_len branch_len len_tag bl_tag _pairs _rest; do
    [[ "$dataset" == "dataset" ]] && continue
    [[ "$model" != "$MODEL" ]] && continue

    outdir="$OUTBASE/${model}/${len_mode}_${len_tag}/${bl_tag}"
    if ! need_backfill_dir "$outdir"; then
        n_skip=$((n_skip + 1))
        continue
    fi

    n_run=$((n_run + 1))
    echo "[$n_run] RUN $dataset  n=$FORCE_PAIRS (setsketch_mle grid)"
    rm -f "$outdir/$SS_MARKER"
    strip_ss_rows_for_dir "$outdir"

    iqt="${MODEL_STR[$model]}"
    if [[ "$len_mode" == "fixed" ]]; then
        gen_alisim_batch_fixed "$model" "$iqt" "" \
            "$nominal_len" "$branch_len" "$outdir" "$FORCE_PAIRS"
    else
        gen_alisim_batch_varlen "$model" "$iqt" "$nominal_len" "$branch_len" "$outdir" "$FORCE_PAIRS"
    fi

    "$SWEEP_EXE" --dir "$outdir" --model "$model" --threads "$THREADS" --algos "$ALGOS" \
        | append_sweep_lines

    rm -f "$outdir"/pair_*.fa "$outdir"/pairs_meta.tsv
    touch "$outdir/$SS_MARKER"
done < "$DATASET_CSV"

echo "Done: ran=$n_run skipped=$n_skip at $(date -Iseconds)"
