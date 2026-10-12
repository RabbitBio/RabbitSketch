#!/usr/bin/env bash
# Re-run synthetic parameter sweeps for FastKMV only (Sections A–D).
#
# Usage:
#   bash run_fastkmv_param_sweep.sh [THREADS] [PAIRS] [LEN]
#   nohup bash run_fastkmv_param_sweep.sh 32 300 1000000 > results/fastkmv_param_sweep.log 2>&1 &
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${1:-32}"
PAIRS="${2:-300}"
LEN="${3:-1000000}"
EXE="./exe_eval_param_sweep"
OUTDIR="results"

mkdir -p "$OUTDIR"
make -s exe_eval_param_sweep
[[ -x "$EXE" ]] || { echo "[ERROR] $EXE missing"; exit 1; }

parse_and_append() {
    local dim="$1" pval="$2" csv_file="$3" mode="${4:-jaccard}"
    while IFS= read -r line; do
        [[ "$line" =~ ^SWEEP\  ]] || continue
        local _algo _kmer _size _bits _ss_a _ss_base _bd_bbits _bd_s64 _drlevel _len _len_b _rate _n _mae _rmse
        _algo=$(    echo "$line" | grep -oP 'algo=\K\S+')
        _kmer=$(    echo "$line" | grep -oP 'kmer=\K[0-9]+')
        _size=$(    echo "$line" | grep -oP 'size=\K[0-9]+')
        _bits=$(    echo "$line" | grep -oP 'bits=\K[0-9]+')
        _ss_a=$(    echo "$line" | grep -oP 'ss_a=\K[0-9.]+')
        _ss_base=$( echo "$line" | grep -oP 'ss_base=\K[0-9.]+')
        _bd_bbits=$(echo "$line" | grep -oP 'bd_bbits=\K[0-9]+')
        _bd_s64=$(  echo "$line" | grep -oP 'bd_s64=\K[0-9]+')
        _drlevel=$( echo "$line" | grep -oP 'drlevel=\K[0-9]+')
        _len=$(     echo "$line" | grep -oP '\blen=\K[0-9]+')
        _len_b=$(   echo "$line" | grep -oP 'len_b=\K[0-9]+')
        _rate=$(    echo "$line" | grep -oP 'rate=\K[0-9.eE+-]+')
        _n=$(       echo "$line" | grep -oP '\bn=\K[0-9]+')
        _mae=$(     echo "$line" | grep -oP 'MAE=\K[0-9.]+')
        _rmse=$(    echo "$line" | grep -oP 'RMSE=\K[0-9.]+')
        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "$mode" "$_algo" "$dim" "$pval" \
            "$_kmer" "$_size" "$_bits" "$_ss_a" "$_ss_base" \
            "$_bd_bbits" "$_bd_s64" "$_drlevel" \
            "$_len" "${_len_b:-$_len}" \
            "$_rate" "$_mae" "$_rmse" >> "$csv_file"
    done
}

run_point() {
    local csv_file="$1" dim="$2" pval="$3"
    shift 3
    local out
    out=$("$EXE" --pairs "$PAIRS" --len "$LEN" --threads "$THREADS" "$@" 2>&1 >/dev/null) || true
    parse_and_append "$dim" "$pval" "$csv_file" "jaccard" <<< "$out"
}

run_point_len() {
    local csv_file="$1" dim="$2" pval="$3" custom_len="$4"
    shift 4
    local out
    out=$("$EXE" --pairs "$PAIRS" --len "$custom_len" --threads "$THREADS" "$@" 2>&1 >/dev/null) || true
    parse_and_append "$dim" "$pval" "$csv_file" "jaccard" <<< "$out"
}

CSV_A="$OUTDIR/param_sweep.csv"
CSV_B="$OUTDIR/crosslen_sweep.csv"
CSV_C="$OUTDIR/equal_memory.csv"
CSV_D="$OUTDIR/lowcomp_sweep.csv"
HDR="mode,algo,sweep_dim,param_value,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,len,len_b,rate,MAE,RMSE"

for f in "$CSV_A" "$CSV_B" "$CSV_C" "$CSV_D"; do
    if [[ -f "$f" ]]; then
        python3 scripts/strip_algo_rows.py --algo fastkmv "$f"
    else
        echo "$HDR" > "$f"
    fi
done

echo "=== FastKMV param sweep: THREADS=$THREADS PAIRS=$PAIRS LEN=$LEN ==="
echo "Start: $(date -Iseconds)"

echo ""
echo "── Section A: size sweep ──"
for s in 64 96 128 192 256 384 512 768 1024 1536 2048 3072 4096 6144 8192; do
    echo -n "  size=$s  "
    run_point "$CSV_A" fastkmv "size" --algo fastkmv --size "$s" --kmer 21
    echo "done"
done

echo "── Section A: kmer sweep ──"
for k in 11 13 15 17 19 21 23 25 27 29 31 33; do
    echo -n "  kmer=$k  "
    run_point "$CSV_A" fastkmv "kmer" --algo fastkmv --size 1024 --kmer "$k"
    echo "done"
done

echo ""
echo "── Section B: cross-length size sweep ──"
CROSS_LENGTHS="1000 5000 20000 100000 500000 1000000 4000000"
for xlen in $CROSS_LENGTHS; do
    tag="L${xlen}"
    echo -n "  $tag  "
    for s in 64 128 256 512 1024 2048 4096 8192; do
        run_point_len "$CSV_B" "size@${tag}" "$s" "$xlen" --algo fastkmv --size "$s" --kmer 21
    done
    echo "done"
done

echo ""
echo "── Section C: equal-memory tiers ──"
declare -A TIER_N TIER_BITS
TIER_N=(  [128]=128   [256]=256   [512]=512   [1024]=1024  [2048]=2048  [4096]=4096  [8192]=8192 )
for tier in 128 256 512 1024 2048 4096 8192; do
    N="${TIER_N[$tier]}"
    echo -n "  tier N=$N  "
    run_point_len "$CSV_C" "mem_tier" "$N" "$LEN" --algo fastkmv --size "$N" --kmer 21
    echo "done"
done

echo ""
echo "── Section D: low-complexity ──"
for frac in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do
    echo -n "  lowcomp_frac=$frac  "
    out=$("$EXE" --pairs "$PAIRS" --len "$LEN" --threads "$THREADS" \
        --mode lowcomp --lowcomp-frac "$frac" --algo fastkmv --size 1024 --kmer 21 \
        2>&1 >/dev/null) || true
    parse_and_append "lowcomp_frac" "$frac" "$CSV_D" "lowcomp" <<< "$out"
    echo "done"
done

echo ""
echo "Done: $(date -Iseconds)"
