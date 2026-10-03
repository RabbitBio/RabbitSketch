#!/usr/bin/env bash
# Re-run only SetSketch configurations on the 200 real genome pairs,
# using the updated rabbitsketch binary (jaccard_index = MLE).
# Replaces the old setsketch_* rows in sketch_results.csv.
set -euo pipefail
cd "$(dirname "$0")"

BIN="${BIN:-/home/user_home/tongz/Finish/RabbitSketch/examples/rabbitsketch}"
PAIRS_TSV="pairs.tsv"
OUT_CSV="sketch_results.csv"
TMP_DIR="/tmp/ss_real_rerun_$$"
PARALLEL=64
THREADS_PER=1

[[ -f "$PAIRS_TSV" ]] || { echo "Need $PAIRS_TSV"; exit 1; }
[[ -x "$BIN" ]]      || { echo "Missing binary: $BIN"; exit 1; }

mkdir -p "$TMP_DIR"

# SetSketch configurations (L1 + L2 from run_sketch_real.sh)
declare -a SS_CONFIGS=(
    "setsketch_def|--setsketch|-k 21 -b 13 -a 5 -B 2"
    "setsketch_a0p5|--setsketch|-k 21 -b 13 -a 0.5 -B 2"
    "setsketch_a1p5|--setsketch|-k 21 -b 13 -a 1.5 -B 2"
    "setsketch_k15|--setsketch|-k 15 -b 13 -a 5 -B 2"
    "setsketch_k27|--setsketch|-k 27 -b 13 -a 5 -B 2"
)

# L1 uses all 200 pairs; L2 uses FastANI-aligned pairs
L2_PAIR_FILE="$TMP_DIR/l2_pairs.txt"
if [[ -f fastani_truth.tsv ]]; then
    awk -F'\t' 'NR>1 && $4!="-1.0" {print $1}' fastani_truth.tsv > "$L2_PAIR_FILE"
else
    awk -F'\t' 'NR>1 {print $1}' "$PAIRS_TSV" > "$L2_PAIR_FILE"
fi
echo "L2 pairs: $(wc -l < "$L2_PAIR_FILE")"

run_one() {
    local pid="$1" layer="$2" tag="$3" algo="$4" extra="$5" fileA="$6" fileB="$7"
    local out="$TMP_DIR/p${pid}_${tag}.out"
    local t0 t1
    t0=$(date +%s.%N)
    "$BIN" "$algo" -i "$fileA" "$fileB" $extra -d 1.0 -t "$THREADS_PER" \
           -o "$out" >/dev/null 2>&1 || true
    t1=$(date +%s.%N)
    local sec
    sec=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')
    local dist="1.0"
    [[ -s "$out" ]] && dist=$(awk -F'\t' '{print $3; exit}' "$out")
    [[ -z "$dist" ]] && dist="1.0"
    local ani
    ani=$(awk -v d="$dist" 'BEGIN{printf "%.6f", 1.0 - d}')
    printf "%s,%s,%s,%s,%s,%s,%s\n" "$pid" "$layer" "$tag" \
           "${extra// /_}" "$dist" "$ani" "$sec"
    rm -f "$out"
}
export -f run_one
export BIN TMP_DIR THREADS_PER

# Build job list
JOBS_FILE="$TMP_DIR/jobs.txt"
> "$JOBS_FILE"

while IFS=$'\t' read -r pid tier accA accB fileA fileB rest; do
    [[ "$pid" == "pair_id" || -z "$pid" ]] && continue
    for cfg in "${SS_CONFIGS[@]}"; do
        tag="${cfg%%|*}"; rest1="${cfg#*|}"
        algo="${rest1%%|*}"; extra="${rest1#*|}"

        # Determine layer
        if [[ "$tag" == "setsketch_def" ]]; then
            layer="L1"
        else
            # L2: only if pid is in L2 pair file
            grep -qx "$pid" "$L2_PAIR_FILE" 2>/dev/null || continue
            layer="L2"
        fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
               "$pid" "$layer" "$tag" "$algo" "$extra" "$fileA" "$fileB" >> "$JOBS_FILE"
    done
done < "$PAIRS_TSV"

N_JOBS=$(wc -l < "$JOBS_FILE")
echo "Total jobs: $N_JOBS  (parallelism=$PARALLEL)"
echo "Binary: $BIN (compiled $(date -r "$BIN" '+%Y-%m-%d %H:%M'))"
echo ""

# Run jobs
RESULTS_DIR="$TMP_DIR/results"
mkdir -p "$RESULTS_DIR"
job_idx=0
while IFS=$'\t' read -r pid layer tag algo extra fileA fileB; do
    job_idx=$((job_idx + 1))
    out_line="$RESULTS_DIR/job_${job_idx}.csv"
    ( run_one "$pid" "$layer" "$tag" "$algo" "$extra" "$fileA" "$fileB" > "$out_line" ) &
    while (( $(jobs -rp | wc -l) >= PARALLEL )); do wait -n; done
    if (( job_idx % 200 == 0 )); then
        echo "  dispatched $job_idx / $N_JOBS"
    fi
done < "$JOBS_FILE"
wait

# Collect new setsketch rows
NEW_SS="$TMP_DIR/new_setsketch.csv"
cat "$RESULTS_DIR"/job_*.csv 2>/dev/null | sort -t, -k1,1n -k2,2 -k3,3 > "$NEW_SS"
N_NEW=$(wc -l < "$NEW_SS")
echo ""
echo "Generated $N_NEW new SetSketch rows"

# Replace old setsketch_* rows in sketch_results.csv
python3 - "$OUT_CSV" "$NEW_SS" <<'PYEOF'
import sys, csv

src, new_ss = sys.argv[1:]

with open(src) as f:
    rows = list(csv.reader(f))
header = rows[0]

# Keep non-setsketch rows
kept = [header] + [r for r in rows[1:] if not r[2].startswith("setsketch")]
old_ss_count = len(rows) - 1 - len(kept) + 1

with open(new_ss) as f:
    new_rows = [row.split(',') for row in f.read().splitlines() if row.strip()]

combined = kept + new_rows
combined[1:] = sorted(combined[1:], key=lambda r: (int(r[0]), r[1], r[2]))

with open(src, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(header)
    for r in combined[1:]:
        w.writerow(r)

print(f"  Replaced {old_ss_count} old setsketch rows with {len(new_rows)} new MLE rows")
print(f"  Total rows now: {len(combined)-1}")
PYEOF

echo ""
echo "Updated $OUT_CSV"
rm -rf "$TMP_DIR"
echo "Done at $(date)"
