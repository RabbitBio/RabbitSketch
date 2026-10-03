#!/usr/bin/env bash
# Run all 7 sketch algorithms on the real-genome pair list (pairs.tsv).
#
# Two layers:
#   L1 = default params (all 200 pairs, all 7 algos)
#   L2 = parameter spot-checks (50-pair subset, key params per algo)
#
# Output:
#   sketch_results.csv
#     columns: pair_id, layer, algo, params, mash_dist, sketch_ani, sec
#
# Sketch-derived ANI is reported as  1 - mash_distance  (the value rabbitsketch
# already returns); rows with dist=1.0 mean "too divergent" / sketch-empty
# intersection, kept for completeness.
set -euo pipefail
cd "$(dirname "$0")"

BIN="${BIN:-/home/user_home/tongz/Finish/RabbitSketch/examples/rabbitsketch}"
PAIRS_TSV="pairs.tsv"
OUT_CSV="sketch_results.csv"
TMPDIR_OUT="${TMPDIR_OUT:-/tmp/sketch_real_$$}"
PARALLEL="${PARALLEL:-16}"
THREADS_PER="${THREADS_PER:-1}"   # threads used by rabbitsketch per pair

[[ -f "$PAIRS_TSV" ]] || { echo "Need $PAIRS_TSV"; exit 1; }
[[ -x "$BIN" ]]      || { echo "rabbitsketch binary not executable: $BIN"; exit 1; }

mkdir -p "$TMPDIR_OUT"
printf "pair_id,layer,algo,params,mash_dist,sketch_ani,sec\n" > "$OUT_CSV"

# ── Configurations to test ───────────────────────────────────────────────────
# Each entry:  "tag|algo|extra-args"
# Layer L1 (defaults — single tag per algo):
L1=(
    "minhash_def|--minhash|-k 21 -s 1000"
    "fastkmv_def|--fastkmv|-k 21 -K 1024"
    "hll_def|--hll|-k 21 -b 13"
    "setsketch_def|--setsketch|-k 21 -b 13 -a 5 -B 2"
    "bindash_def|--bindash|-k 21 --bbits 16 --sketch64 16"
    "kssd_def|--kssd|--drlevel 3"
    "probmh_def|--probminhash|-k 21 -m 1024"
)

# Layer L2 (parameter spot-checks — informed by alisim findings):
#   k-mer sensitivity:  k=15  and  k=27  (vs default k=21)
#   SetSketch a:        a=0.5 and  a=1.5 (vs default a=5)  [short/medium fill]
#   KSSD drlevel:       drlevel=1 and 2  (vs default drlevel=3)
#   sketch size:        s=256 and s=4096 (vs default 1000/1024)
L2=(
    # k-mer sensitivity (one per algo using k=15 / k=27)
    "minhash_k15|--minhash|-k 15 -s 1000"
    "minhash_k27|--minhash|-k 27 -s 1000"
    "fastkmv_k15|--fastkmv|-k 15 -K 1024"
    "fastkmv_k27|--fastkmv|-k 27 -K 1024"
    "hll_k15|--hll|-k 15 -b 13"
    "hll_k27|--hll|-k 27 -b 13"
    "bindash_k15|--bindash|-k 15 --bbits 16 --sketch64 16"
    "bindash_k27|--bindash|-k 27 --bbits 16 --sketch64 16"
    "probmh_k15|--probminhash|-k 15 -m 1024"
    "probmh_k27|--probminhash|-k 27 -m 1024"

    # SetSketch a-sensitivity (key alisim finding)
    "setsketch_a0p5|--setsketch|-k 21 -b 13 -a 0.5 -B 2"
    "setsketch_a1p5|--setsketch|-k 21 -b 13 -a 1.5 -B 2"
    "setsketch_k15|--setsketch|-k 15 -b 13 -a 5 -B 2"
    "setsketch_k27|--setsketch|-k 27 -b 13 -a 5 -B 2"

    # KSSD drlevel sensitivity
    "kssd_dr1|--kssd|--drlevel 1"
    "kssd_dr2|--kssd|--drlevel 2"

    # Sketch size sensitivity
    "minhash_s256|--minhash|-k 21 -s 256"
    "minhash_s4096|--minhash|-k 21 -s 4096"
    "probmh_s256|--probminhash|-k 21 -m 256"
    "probmh_s4096|--probminhash|-k 21 -m 4096"
)

# ── How many L2 pairs to use ─────────────────────────────────────────────────
# Take first 50 pairs that have a successful FastANI value (if file present),
# otherwise just first 50 pair IDs.
L2_PAIR_FILE="$TMPDIR_OUT/l2_pair_ids.txt"
# Use ALL FastANI-aligned pairs (stratified across tiers) for parameter sweeps.
# This avoids the "first-N is all Tier A" sampling bias and gives clean
# parameter-vs-tier breakdowns.
if [[ -f fastani_truth.tsv ]]; then
    awk -F'\t' 'NR>1 && $4!="-1.0" {print $1}' fastani_truth.tsv > "$L2_PAIR_FILE"
else
    awk -F'\t' 'NR>1 {print $1}' "$PAIRS_TSV" > "$L2_PAIR_FILE"
fi
echo "L2 will use $(wc -l < "$L2_PAIR_FILE") pairs (all FastANI-aligned)"

# ── Worker function (one pair × one config) ──────────────────────────────────
run_one() {
    local pid="$1" layer="$2" tag="$3" algo="$4" extra="$5" fileA="$6" fileB="$7"
    local out="$TMPDIR_OUT/p${pid}_${tag}.out"
    local t0 t1 dist
    t0=$(date +%s.%N)
    "$BIN" "$algo" -i "$fileA" "$fileB" $extra -d 1.0 -t "$THREADS_PER" \
           -o "$out" >/dev/null 2>&1 || true
    t1=$(date +%s.%N)
    local sec
    sec=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')
    if [[ -s "$out" ]]; then
        dist=$(awk -F'\t' '{print $3; exit}' "$out")
    else
        dist="1.0"
    fi
    [[ -z "$dist" ]] && dist="1.0"
    local ani
    ani=$(awk -v d="$dist" 'BEGIN{printf "%.6f", 1.0 - d}')
    printf "%s,%s,%s,%s,%s,%s,%s\n" "$pid" "$layer" "$tag" \
           "${extra// /_}" "$dist" "$ani" "$sec"
    rm -f "$out"
}
export -f run_one
export BIN TMPDIR_OUT THREADS_PER

# ── Dispatch jobs ────────────────────────────────────────────────────────────
# We stream jobs to stdout via parallel bash &; bounded by PARALLEL jobs.
jobs_file="$TMPDIR_OUT/jobs.txt"
: > "$jobs_file"

# emit one job line per (pair × config) into jobs_file
emit_jobs() {
    local layer="$1"; shift
    local -n configs=$1; shift
    local pair_filter="${1:-}"
    while IFS=$'\t' read -r pid tier accA accB fileA fileB rest; do
        [[ "$pid" == "pair_id" ]] && continue
        [[ -z "$pid" ]] && continue
        if [[ -n "$pair_filter" ]] && ! grep -qx "$pid" "$pair_filter"; then
            continue
        fi
        for cfg in "${configs[@]}"; do
            tag="${cfg%%|*}"; rest1="${cfg#*|}"
            algo="${rest1%%|*}"; extra="${rest1#*|}"
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                   "$pid" "$layer" "$tag" "$algo" "$extra" "$fileA" "$fileB" \
                   >> "$jobs_file"
        done
    done < "$PAIRS_TSV"
}

emit_jobs L1 L1 ""
emit_jobs L2 L2 "$L2_PAIR_FILE"

n_jobs=$(wc -l < "$jobs_file")
echo "Total jobs to run: $n_jobs  (parallelism=$PARALLEL)"
echo "Time started: $(date)"
echo ""

# Bounded-parallel bash loop: each job writes one CSV line to per-pid file,
# we concat at the end (avoids interleaved-line race conditions in append).
results_dir="$TMPDIR_OUT/results"
mkdir -p "$results_dir"
job_idx=0
while IFS=$'\t' read -r pid layer tag algo extra fileA fileB; do
    job_idx=$((job_idx + 1))
    out_line="$results_dir/job_${job_idx}.csv"
    (
        run_one "$pid" "$layer" "$tag" "$algo" "$extra" "$fileA" "$fileB" > "$out_line"
    ) &
    while (( $(jobs -rp | wc -l) >= PARALLEL )); do wait -n; done
    if (( job_idx % 100 == 0 )); then
        echo "  dispatched $job_idx / $n_jobs"
    fi
done < "$jobs_file"
wait

# Concatenate in pair-id, layer order (just numeric pair_id sort)
cat "$results_dir"/job_*.csv 2>/dev/null | sort -t, -k1,1n -k2,2 -k3,3 >> "$OUT_CSV"

echo "Done at $(date)"
echo "Results: $OUT_CSV  ($(($(wc -l < "$OUT_CSV") - 1)) rows)"
rm -rf "$TMPDIR_OUT"
