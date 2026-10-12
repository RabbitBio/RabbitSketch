#!/usr/bin/env bash
# Backfill JC69/fixed_L1k/d0001 (BL=0.0005) @ 1000 pairs, append to selected6 CSV.
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${1:-64}"
MODEL="JC69"
NOMINAL_LEN=1000
BRANCH_LEN=0.0005
FORCE_PAIRS=1000
ALGOS="fastkmv,hll,setsketch_mle,bindash,probmh,kssd"

OUTBASE="/home/bigssd/zt/ref"
OUTDIR="$OUTBASE/JC69/fixed_L1k/d0001"
OUT_CSV="results/alisim_sweep.JC69.pairs1000.selected6.csv"
SWEEP_EXE="./exe_eval_alisim_sweep"
MARKER="$OUTDIR/.pairs1000_selected6_JC69_done"

for _iqbin in "${HOME}/anaconda3/envs/cpp20/bin" \
               "/home/user_home/tongz/anaconda3/envs/cpp20/bin"; do
    [[ -d "$_iqbin" ]] && export PATH="$_iqbin:$PATH" && break
done

IQTREE=""
for cand in iqtree2 iqtree3 iqtree; do
    if p="$(command -v "$cand" 2>/dev/null)" && [[ -n "$p" && -x "$p" ]]; then IQTREE="$p"; break; fi
done
[[ -x "$SWEEP_EXE" ]] || { echo "[ERROR] $SWEEP_EXE missing"; exit 1; }
[[ -n "$IQTREE" ]] || { echo "[ERROR] iqtree not found"; exit 1; }

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

echo "Supplement: $MODEL/fixed_L1k/d0001  n=$FORCE_PAIRS  BL=$BRANCH_LEN"
rm -f "$MARKER"

tmpd="/tmp/alisim_supp_l1k_d0001_$$"
mkdir -p "$tmpd" "$OUTDIR"
printf "(A:%.6f,B:%.6f);\n" "$BRANCH_LEN" "$BRANCH_LEN" > "$tmpd/tree2.nwk"

pids=() slot=0
for ((i=0; i<FORCE_PAIRS; i++)); do
    iqtree3 --alisim "${tmpd}/p${i}_sim" -m "JC" -t "$tmpd/tree2.nwk" \
        --length "$NOMINAL_LEN" --out-format fasta --seed $((42+i)) -nt 1 -quiet \
        >/dev/null 2>&1 &
    pids+=($!)
    slot=$((slot + 1))
    if (( slot >= THREADS )); then
        wait "${pids[0]}" || true; pids=("${pids[@]:1}"); slot=$((slot - 1))
    fi
done
for p in "${pids[@]}"; do wait "$p" || true; done
_split_and_meta "$tmpd" "$OUTDIR" "$BRANCH_LEN" "$FORCE_PAIRS"
rm -rf "$tmpd"

"$SWEEP_EXE" --dir "$OUTDIR" --model "$MODEL" --threads "$THREADS" --algos "$ALGOS" \
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

rm -f "$OUTDIR"/pair_*.fa "$OUTDIR"/pairs_meta.tsv
touch "$MARKER"
echo "Appended to $OUT_CSV"
