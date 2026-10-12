#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_param_sweep_all.sh  –  comprehensive evaluation suite
#
# Section A: Parameter sensitivity (default length, Jaccard mode)
#   – for each algorithm: sweep every key parameter over a fine grid
# Section B: Parameter × Length cross (primary size param, all lengths)
#   – reveals how the optimal parameter changes with sequence length
# Section C: Equal-memory fair comparison
#   – all algorithms compared at the same memory budget (per-sketch bytes)
# Section D: Low-complexity contamination
#   – all algorithms at varying lowcomp fraction, default params
#
# Output:
#   results/param_sweep.csv         (Section A)
#   results/crosslen_sweep.csv      (Section B)
#   results/equal_memory.csv        (Section C)
#   results/lowcomp_sweep.csv       (Section D)
#
# Usage:
#   bash run_param_sweep_all.sh [THREADS] [PAIRS] [LEN]
#   THREADS : OpenMP threads per job       (default: 32)
#   PAIRS   : pairs per mutation rate      (default: 300)
#   LEN     : default sequence length (A)  (default: 1000000)
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${1:-32}"
PAIRS="${2:-300}"
LEN="${3:-1000000}"
EXE="./exe_eval_param_sweep"
OUTDIR="results"

mkdir -p "$OUTDIR"

if [[ ! -x "$EXE" ]]; then
    echo "[ERROR] $EXE not found – run 'make exe_eval_param_sweep' first"
    exit 1
fi

# ─── generic parser helper ────────────────────────────────────────────────────
# parse_sweep_line: extracts key=value fields from one SWEEP line → CSV row
# args: dim pval csv_file
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

# ─── run_point: run exe + capture + append ────────────────────────────────────
# args: csv_file dim pval [--mode M] --algo A [algo-args...]
run_point() {
    local csv_file="$1" dim="$2" pval="$3"
    shift 3
    local mode_flag="jaccard"
    # peek at args for --mode
    for arg in "$@"; do
        [[ "$arg" == "jaccard" || "$arg" == "containment" || "$arg" == "lowcomp" ]] && mode_flag="$arg"
    done
    local out
    out=$("$EXE" --pairs "$PAIRS" --len "$LEN" --threads "$THREADS" "$@" 2>&1 >/dev/null) || true
    parse_and_append "$dim" "$pval" "$csv_file" "$mode_flag" <<< "$out"
}

# Same but with a custom length
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

HDR_MAIN="mode,algo,sweep_dim,param_value,kmer,size,bits,ss_a,ss_base,bd_bbits,bd_s64,drlevel,len,len_b,rate,MAE,RMSE"
for f in "$CSV_A" "$CSV_B" "$CSV_C" "$CSV_D"; do
    : > "$f"
    echo "$HDR_MAIN" >> "$f"
done

echo "=== run_param_sweep_all: THREADS=$THREADS PAIRS=$PAIRS LEN=$LEN ==="
echo "Start: $(date -Iseconds)"

# ═══════════════════════════════════════════════════════════════════════════════
# MinHash  –  sweep: sketch size s, k-mer size k
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "── MinHash: sketch size sweep ──"
for s in 64 96 128 192 256 384 512 768 1024 1536 2048 3072 4096 6144 8192; do
    echo -n "  size=$s  "
    run_point "$CSV_A" minhash "size" --algo minhash "$s" --size "$s" --kmer 21
    echo "done"
done

echo "── MinHash: k-mer size sweep ──"
for k in 11 13 15 17 19 21 23 25 27 29 31 33; do
    echo -n "  kmer=$k  "
    run_point "$CSV_A" minhash "kmer" --algo minhash "$k" --size 1024 --kmer "$k"
    echo "done"
done

# ═══════════════════════════════════════════════════════════════════════════════
# FastKMV  –  sweep: sketch size k, k-mer size
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "── FastKMV: sketch size sweep ──"
for s in 64 96 128 192 256 384 512 768 1024 1536 2048 3072 4096 6144 8192; do
    echo -n "  size=$s  "
    run_point "$CSV_A" fastkmv "size" --algo fastkmv "$s" --size "$s" --kmer 21
    echo "done"
done

echo "── FastKMV: k-mer size sweep ──"
for k in 11 13 15 17 19 21 23 25 27 29 31 33; do
    echo -n "  kmer=$k  "
    run_point "$CSV_A" fastkmv "kmer" --algo fastkmv "$k" --size 1024 --kmer "$k"
    echo "done"
done

# ═══════════════════════════════════════════════════════════════════════════════
# HLL  –  sweep: bits (np), k-mer size
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "── HLL: bits (np) sweep ──"
for b in 8 9 10 11 12 13 14 15 16; do
    echo -n "  bits=$b  "
    run_point "$CSV_A" hll "bits" --algo hll "$b" --bits "$b" --kmer 32
    echo "done"
done

echo "── HLL: k-mer size sweep ──"
for k in 15 17 19 21 23 25 27 29 31 33 35 37; do
    echo -n "  kmer=$k  "
    run_point "$CSV_A" hll "kmer" --algo hll "$k" --bits 13 --kmer "$k"
    echo "done"
done

# ═══════════════════════════════════════════════════════════════════════════════
# SetSketch  –  sweep: bits, a, base, k-mer
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "── SetSketch: bits sweep ──"
for b in 8 9 10 11 12 13 14 15 16; do
    echo -n "  bits=$b  "
    run_point "$CSV_A" setsketch "bits" --algo setsketch "$b" --bits "$b" --ss-a 5.0 --ss-base 2.0 --kmer 32
    echo "done"
done

echo "── SetSketch: a parameter sweep ──"
for a in 1 2 3 5 7 10 15 20 30 50 80 120; do
    echo -n "  ss_a=$a  "
    run_point "$CSV_A" setsketch "ss_a" --algo setsketch "$a" --bits 13 --ss-a "$a" --ss-base 2.0 --kmer 32
    echo "done"
done

echo "── SetSketch: base parameter sweep ──"
for base in 1.2 1.5 1.8 2.0 2.5 3.0 4.0 5.0 8.0; do
    echo -n "  ss_base=$base  "
    run_point "$CSV_A" setsketch "ss_base" --algo setsketch "$base" --bits 13 --ss-a 5.0 --ss-base "$base" --kmer 32
    echo "done"
done

echo "── SetSketch: k-mer size sweep ──"
for k in 15 17 19 21 23 25 27 29 31 33 35 37; do
    echo -n "  kmer=$k  "
    run_point "$CSV_A" setsketch "kmer" --algo setsketch "$k" --bits 13 --ss-a 5.0 --ss-base 2.0 --kmer "$k"
    echo "done"
done

# ═══════════════════════════════════════════════════════════════════════════════
# BinDash  –  sweep: bbits, sketch64, k-mer
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "── BinDash: bbits sweep ──"
for bb in 2 4 6 8 10 12 14 16 20 24 32; do
    echo -n "  bd_bbits=$bb  "
    run_point "$CSV_A" bindash "bd_bbits" --algo bindash "$bb" --bd-bbits "$bb" --bd-s64 16 --kmer 21
    echo "done"
done

echo "── BinDash: sketch64 (bins) sweep ──"
for s64 in 2 4 6 8 10 12 16 20 24 32; do
    echo -n "  bd_s64=$s64  "
    run_point "$CSV_A" bindash "bd_s64" --algo bindash "$s64" --bd-bbits 16 --bd-s64 "$s64" --kmer 21
    echo "done"
done

echo "── BinDash: k-mer size sweep ──"
for k in 11 13 15 17 19 21 23 25 27 29 31 33; do
    echo -n "  kmer=$k  "
    run_point "$CSV_A" bindash "kmer" --algo bindash "$k" --bd-bbits 16 --bd-s64 16 --kmer "$k"
    echo "done"
done

# ═══════════════════════════════════════════════════════════════════════════════
# KSSD  –  sweep: drlevel (only user-controllable parameter)
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "── KSSD: drlevel sweep ──"
for d in 1 2 3 4 5 6; do
    echo -n "  drlevel=$d  "
    run_point "$CSV_A" kssd "drlevel" --algo kssd "$d" --drlevel "$d"
    echo "done"
done

# ═══════════════════════════════════════════════════════════════════════════════
# ProbMinHash  –  sweep: m (register count), k-mer
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "── ProbMinHash: register count (m) sweep ──"
for m in 64 96 128 192 256 384 512 768 1024 1536 2048 3072 4096; do
    echo -n "  m=$m  "
    run_point "$CSV_A" probmh "size" --algo probmh "$m" --size "$m" --kmer 21
    echo "done"
done

echo "── ProbMinHash: k-mer size sweep ──"
for k in 11 13 15 17 19 21 23 25 27 29 31 33; do
    echo -n "  kmer=$k  "
    run_point "$CSV_A" probmh "kmer" "$k" --algo probmh --size 1024 --kmer "$k"
    echo "done"
done

echo ""
echo "=== Section A done: $CSV_A ==="

# ═══════════════════════════════════════════════════════════════════════════════
# Section B – Parameter × Length cross sweep
# For each algorithm, sweep primary size parameter over ALL sequence lengths.
# This reveals how optimal parameter changes with sequence length.
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "╔══════════════════════════════════════════════════════╗"
echo "║  Section B: Parameter × Length cross sweep           ║"
echo "╚══════════════════════════════════════════════════════╝"

CROSS_LENGTHS="1000 5000 20000 100000 500000 1000000 4000000"

for xlen in $CROSS_LENGTHS; do
    tag="L${xlen}"
    echo ""
    echo "── Length = $xlen bp ──"

    echo -n "  MinHash size  "
    for s in 64 128 256 512 1024 2048 4096 8192; do
        run_point_len "$CSV_B" "size@$tag" "$s" "$xlen" --algo minhash --size "$s" --kmer 21
    done; echo "done"

    echo -n "  FastKMV size  "
    for s in 64 128 256 512 1024 2048 4096 8192; do
        run_point_len "$CSV_B" "size@$tag" "$s" "$xlen" --algo fastkmv --size "$s" --kmer 21
    done; echo "done"

    echo -n "  HLL bits  "
    for b in 8 9 10 11 12 13 14 15 16; do
        run_point_len "$CSV_B" "bits@$tag" "$b" "$xlen" --algo hll --bits "$b" --kmer 32
    done; echo "done"

    echo -n "  SetSketch bits  "
    for b in 8 9 10 11 12 13 14 15 16; do
        run_point_len "$CSV_B" "bits@$tag" "$b" "$xlen" --algo setsketch --bits "$b" --ss-a 5.0 --ss-base 2.0 --kmer 32
    done; echo "done"

    echo -n "  BinDash s64  "
    for s64 in 2 4 8 12 16 24 32; do
        run_point_len "$CSV_B" "bd_s64@$tag" "$s64" "$xlen" --algo bindash --bd-s64 "$s64" --bd-bbits 16 --kmer 21
    done; echo "done"

    echo -n "  KSSD drlevel  "
    for d in 1 2 3 4 5 6; do
        run_point_len "$CSV_B" "drlevel@$tag" "$d" "$xlen" --algo kssd --drlevel "$d"
    done; echo "done"

    echo -n "  PMH size  "
    for m in 64 128 256 512 1024 2048 4096; do
        run_point_len "$CSV_B" "size@$tag" "$m" "$xlen" --algo probmh --size "$m" --kmer 21
    done; echo "done"
done

echo ""
echo "=== Section B done: $CSV_B ==="

# ═══════════════════════════════════════════════════════════════════════════════
# Section C – Equal-memory fair comparison
# All algorithms compared at the same per-sketch memory tier.
# Tier defined by "equivalent register/hash count" N.
# Approximate sketch sizes:
#   MinHash(N):    N × 8 bytes
#   FastKMV(N):    N × 8 bytes
#   HLL(bits):     2^bits bytes  →  bits = log2(N)
#   SetSketch(bits): same as HLL
#   BinDash(s64, bbits=16):  s64 × 64 × 2 bytes  →  s64 = N/64 (min 2)
#   ProbMinHash(m): m × 8 bytes
#   KSSD:          use nearest natural drlevel (reference point only)
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "╔══════════════════════════════════════════════════════╗"
echo "║  Section C: Equal-memory fair comparison             ║"
echo "╚══════════════════════════════════════════════════════╝"

# Tiers: (tier_name  N  HLL_bits  bd_s64)
declare -A TIER_N TIER_BITS TIER_S64
TIER_N=(  [128]=128   [256]=256   [512]=512   [1024]=1024  [2048]=2048  [4096]=4096  [8192]=8192 )
TIER_BITS=([128]=7    [256]=8     [512]=9     [1024]=10    [2048]=11    [4096]=12    [8192]=13   )
TIER_S64=( [128]=2    [256]=4     [512]=8     [1024]=16    [2048]=32    [4096]=64    [8192]=128  )

for tier in 128 256 512 1024 2048 4096 8192; do
    N="${TIER_N[$tier]}"
    bits="${TIER_BITS[$tier]}"
    s64="${TIER_S64[$tier]}"
    mem_kb=$(( N * 8 / 1024 > 0 ? N * 8 / 1024 : 1 ))
    echo ""
    echo "── Tier N=$N (~${mem_kb} KB per sketch) ──"

    run_point_len "$CSV_C" "mem_tier" "$N" "$LEN" --algo minhash    --size "$N"   --kmer 21
    run_point_len "$CSV_C" "mem_tier" "$N" "$LEN" --algo fastkmv    --size "$N"   --kmer 21
    run_point_len "$CSV_C" "mem_tier" "$N" "$LEN" --algo hll        --bits "$bits" --kmer 32
    run_point_len "$CSV_C" "mem_tier" "$N" "$LEN" --algo setsketch  --bits "$bits" --ss-a 5.0 --ss-base 2.0 --kmer 32
    run_point_len "$CSV_C" "mem_tier" "$N" "$LEN" --algo bindash    --bd-s64 "$s64" --bd-bbits 16 --kmer 21
    run_point_len "$CSV_C" "mem_tier" "$N" "$LEN" --algo probmh     --size "$N"   --kmer 21
    echo "  done"
done
# KSSD at natural sizes (drlevel 2/3/4 for reference)
for d in 2 3 4; do
    run_point_len "$CSV_C" "mem_tier" "kssd_d${d}" "$LEN" --algo kssd --drlevel "$d"
done

echo ""
echo "=== Section C done: $CSV_C ==="

# ═══════════════════════════════════════════════════════════════════════════════
# Section D – Low-complexity contamination
# Sequences have a low-complexity biased prefix.
# Sweep lowcomp_frac from 0 (clean) to 0.9 (mostly low-complexity).
# Uses default params for all algorithms.
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "╔══════════════════════════════════════════════════════╗"
echo "║  Section D: Low-complexity contamination sweep       ║"
echo "╚══════════════════════════════════════════════════════╝"

for frac in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do
    echo -n "  lowcomp_frac=$frac  "
    for algo in minhash fastkmv hll setsketch bindash kssd probmh; do
        _lc_args=(--mode lowcomp --lowcomp-frac "$frac" --algo "$algo")
        case "$algo" in
            hll|setsketch) _lc_args+=(--bits 13 --kmer 32) ;;
            minhash|fastkmv|probmh) _lc_args+=(--size 1024 --kmer 21) ;;
            bindash)  _lc_args+=(--bd-s64 16 --bd-bbits 16 --kmer 21) ;;
            kssd)     _lc_args+=(--drlevel 3) ;;
        esac
        _lc_out=$("$EXE" --pairs "$PAIRS" --len "$LEN" --threads "$THREADS" "${_lc_args[@]}" 2>&1 >/dev/null) || true
        parse_and_append "lowcomp_frac" "$frac" "$CSV_D" "lowcomp" <<< "$_lc_out"
    done
    echo "done"
done

echo ""
echo "=== Section D done: $CSV_D ==="

echo ""
echo "╔══════════════════════════════════════════════════════╗"
echo "║  All sections complete                               ║"
echo "╚══════════════════════════════════════════════════════╝"
echo "Done: $(date -Iseconds)"
echo "Section A (param sweep):        $CSV_A"
echo "Section B (cross-length):       $CSV_B"
echo "Section C (equal-memory):       $CSV_C"
echo "Section D (low-complexity):     $CSV_D"
