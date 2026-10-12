#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_containment_sweep.sh
#
# Systematic evaluation of containment estimation accuracy across:
#   1. Length-ratio sweep:   fixed len_A=1M, vary len_B from 1k to 1M
#   2. Sketch-size × ratio:  for each algo, vary size AND length ratio
#   3. Mutation rate × ratio: how containment accuracy degrades under mutations
#
# All 7 algorithms evaluated (KSSD uses Jaccard as proxy since no native
# containment method).
#
# Output:
#   results/containment_ratio.csv    (length-ratio sweep, default params)
#   results/containment_size.csv     (sketch-size × ratio for each algo)
#   results/containment_rate.csv     (mutation rate × ratio)
#
# Usage:
#   bash run_containment_sweep.sh [THREADS] [PAIRS]
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${1:-32}"
PAIRS="${2:-300}"
LEN_A=1000000
EXE="./exe_eval_param_sweep"
OUTDIR="results"

mkdir -p "$OUTDIR"

if [[ ! -x "$EXE" ]]; then
    echo "[ERROR] $EXE not found – run 'make exe_eval_param_sweep' first"
    exit 1
fi

CSV_RATIO="$OUTDIR/containment_ratio.csv"
CSV_SIZE="$OUTDIR/containment_size.csv"
CSV_RATE="$OUTDIR/containment_rate.csv"

HDR="mode,algo,sweep_dim,param_value,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,len_a,len_b,rate,MAE,RMSE"
for f in "$CSV_RATIO" "$CSV_SIZE" "$CSV_RATE"; do
    : > "$f"; echo "$HDR" >> "$f"
done

# ─── parser ──────────────────────────────────────────────────────────────────
parse_containment() {
    local dim="$1" pval="$2" csv_file="$3"
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
        _mae=$(     echo "$line" | grep -oP 'MAE=\K[0-9.]+')
        _rmse=$(    echo "$line" | grep -oP 'RMSE=\K[0-9.]+')
        printf 'containment,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "$_algo" "$dim" "$pval" \
            "$_kmer" "$_size" "$_bits" "$_ss_a" "$_ss_base" \
            "$_bd_bbits" "$_bd_s64" "$_drlevel" \
            "$_len" "${_len_b:-$_len}" \
            "$_rate" "$_mae" "$_rmse" >> "$csv_file"
    done
}

# Run one containment point for a given algo + len_b + extra args
run_cont() {
    local csv_file="$1" dim="$2" pval="$3" algo="$4" len_b="$5"
    shift 5
    local out
    out=$("$EXE" --mode containment --algo "$algo" \
              --len "$LEN_A" --len-b "$len_b" \
              --pairs "$PAIRS" --threads "$THREADS" "$@" 2>&1 >/dev/null) || true
    parse_containment "$dim" "$pval" "$csv_file" <<< "$out"
}

# ═══════════════════════════════════════════════════════════════════════════════
# Part 1: Length-ratio sweep  (default params, all 7 algos)
# Vary len_B / len_A from 0.001 to 1.0 (B is a subregion of A)
# ═══════════════════════════════════════════════════════════════════════════════
echo "=== Part 1: Length-ratio sweep (default params) ==="

LENGTH_BS=(1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000)

for len_b in "${LENGTH_BS[@]}"; do
    ratio=$(echo "scale=4; $len_b / $LEN_A" | bc)
    echo -n "  len_B=$len_b (ratio=$ratio)  "

    run_cont "$CSV_RATIO" "len_b" "$len_b" minhash    "$len_b" --size 1024  --kmer 21
    run_cont "$CSV_RATIO" "len_b" "$len_b" fastkmv    "$len_b" --size 1024  --kmer 21
    run_cont "$CSV_RATIO" "len_b" "$len_b" hll        "$len_b" --bits 13    --kmer 32
    run_cont "$CSV_RATIO" "len_b" "$len_b" setsketch  "$len_b" --bits 13    --ss-a 5.0 --ss-base 2.0 --kmer 32
    run_cont "$CSV_RATIO" "len_b" "$len_b" bindash    "$len_b" --bd-s64 16  --bd-bbits 16 --kmer 21
    run_cont "$CSV_RATIO" "len_b" "$len_b" kssd       "$len_b" --drlevel 3
    run_cont "$CSV_RATIO" "len_b" "$len_b" probmh     "$len_b" --size 1024  --kmer 21

    echo "done"
done

echo "Part 1 done: $CSV_RATIO"

# ═══════════════════════════════════════════════════════════════════════════════
# Part 2: Sketch-size × length-ratio  (key algos: MinHash, FastKMV, PMH)
# Shows whether larger sketches recover containment accuracy for short B
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "=== Part 2: Sketch-size × length-ratio ==="

SIZES=(128 256 512 1024 2048 4096)
# Focus on the challenging short-B regime
SHORT_LEN_BS=(1000 5000 20000 100000 500000)

for len_b in "${SHORT_LEN_BS[@]}"; do
    echo "  len_B=$len_b:"
    for s in "${SIZES[@]}"; do
        echo -n "    size=$s  "
        run_cont "$CSV_SIZE" "size@lenB${len_b}" "$s" minhash   "$len_b" --size "$s" --kmer 21
        run_cont "$CSV_SIZE" "size@lenB${len_b}" "$s" fastkmv   "$len_b" --size "$s" --kmer 21
        run_cont "$CSV_SIZE" "size@lenB${len_b}" "$s" probmh    "$len_b" --size "$s" --kmer 21
        echo "done"
    done

    # HLL/SetSketch bits sweep
    for bits in 8 10 12 13 14 16; do
        echo -n "    bits=$bits  "
        run_cont "$CSV_SIZE" "bits@lenB${len_b}" "$bits" hll       "$len_b" --bits "$bits" --kmer 32
        run_cont "$CSV_SIZE" "bits@lenB${len_b}" "$bits" setsketch "$len_b" --bits "$bits" --ss-a 5.0 --ss-base 2.0 --kmer 32
        echo "done"
    done
done

echo "Part 2 done: $CSV_SIZE"

# ═══════════════════════════════════════════════════════════════════════════════
# Part 3: Mutation rate × length-ratio
# At fixed default params, 3×3 grid: len_B ∈ {10k,100k,1M} × rate ∈ {0,0.01,0.05,0.10}
# Shows how mutation + short-B degrades containment jointly
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "=== Part 3: Mutation rate × length-ratio ==="

CONT_LEN_BS=(10000 100000 500000 1000000)

for len_b in "${CONT_LEN_BS[@]}"; do
    echo -n "  len_B=$len_b  "
    for algo in minhash fastkmv hll setsketch bindash kssd probmh; do
        case "$algo" in
            hll|setsketch) args=(--bits 13 --kmer 32) ;;
            minhash|fastkmv|probmh) args=(--size 1024 --kmer 21) ;;
            bindash) args=(--bd-s64 16 --bd-bbits 16 --kmer 21) ;;
            kssd)    args=(--drlevel 3) ;;
        esac
        run_cont "$CSV_RATE" "rate×lenB${len_b}" "$len_b" "$algo" "$len_b" "${args[@]}"
    done
    echo "done"
done

echo "Part 3 done: $CSV_RATE"

echo ""
echo "=== Containment sweep complete ==="
echo "Part 1 (ratio sweep):      $CSV_RATIO"
echo "Part 2 (size × ratio):     $CSV_SIZE"
echo "Part 3 (rate × ratio):     $CSV_RATE"
