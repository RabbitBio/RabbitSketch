#!/usr/bin/env bash
# Compute FastANI ground-truth for every pair in pairs.tsv.
# FastANI is invoked once per pair (simple & robust); ~200 pairs * <2 s ≈ 6 min.
#
# Output:  fastani_truth.tsv
#   columns: pair_id, accA, accB, fastani, frac_aligned, total_frags
#   fastani = -1.0 if FastANI returned no result (too divergent).
set -euo pipefail
cd "$(dirname "$0")"

PAIRS_TSV="pairs.tsv"
OUT_TSV="fastani_truth.tsv"
THREADS="${THREADS:-8}"        # how many pair-runs to run in parallel
FASTANI="${FASTANI:-fastANI}"

[[ -f "$PAIRS_TSV" ]] || { echo "Need $PAIRS_TSV first (run select_real_pairs.py)"; exit 1; }
command -v "$FASTANI" >/dev/null || { echo "fastANI not found in PATH"; exit 1; }

n_total=$(($(wc -l < "$PAIRS_TSV") - 1))
echo "Running FastANI on $n_total pairs, parallelism=$THREADS"
echo "FastANI: $($FASTANI --version 2>&1 | head -1)"
echo ""

work_dir="$(mktemp -d -t fastani_work_XXXXXX)"
echo "scratch: $work_dir"

# header
printf "pair_id\taccA\taccB\tfastani\tfrac_aligned\ttotal_frags\n" > "$OUT_TSV"

run_one() {
    local pid="$1" accA="$2" accB="$3" fileA="$4" fileB="$5"
    local out="$work_dir/p${pid}.out"
    # FastANI with single thread per pair; we parallelise across pairs.
    "$FASTANI" -q "$fileA" -r "$fileB" -o "$out" --threads 1 \
               --minFraction 0.0 >/dev/null 2>"$work_dir/p${pid}.err" || true
    if [[ -s "$out" ]]; then
        # FastANI output: query ref ANI frac_aligned total_frags
        awk -v pid="$pid" -v a="$accA" -v b="$accB" \
            '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", pid, a, b, $3, $4, $5}' "$out"
    else
        printf "%s\t%s\t%s\t-1.0\t0\t0\n" "$pid" "$accA" "$accB"
    fi
}
export -f run_one
export work_dir FASTANI

# Read pairs, run in parallel (bounded job count), collect results
tmp_out="$work_dir/all_results.tsv"
: > "$tmp_out"
while IFS=$'\t' read -r pid tier accA accB fileA fileB rest; do
    [[ -z "$pid" ]] && continue
    {
        run_one "$pid" "$accA" "$accB" "$fileA" "$fileB" >> "$tmp_out"
    } &
    while (( $(jobs -rp | wc -l) >= THREADS )); do wait -n; done
done < <(tail -n +2 "$PAIRS_TSV")
wait

sort -k1,1 -n "$tmp_out" >> "$OUT_TSV"

# Summary
n_ok=$(awk -F'\t' 'NR>1 && $4!="-1.0" {n++} END{print n+0}' "$OUT_TSV")
n_no=$(awk -F'\t' 'NR>1 && $4=="-1.0" {n++} END{print n+0}' "$OUT_TSV")
echo ""
echo "Done."
echo "  FastANI succeeded:  $n_ok pairs"
echo "  FastANI no result:  $n_no pairs (too divergent)"
echo "  Output:  $OUT_TSV"

# Cleanup scratch
rm -rf "$work_dir"
