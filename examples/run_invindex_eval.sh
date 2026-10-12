#!/usr/bin/env bash
# run_invindex_eval.sh
# Evaluate inverted-index-based sketch tools (FastKMV, BinDash, SetSketch,
# ProbMinHash) on disk-based benchmark sequences.
#
# Approach:
#   1. Combine A_list.txt + B_list.txt into one file list
#   2. Run each tool in all-to-all inverted-index mode
#   3. Extract matching pair results (pair{id}_A ↔ pair{id}_B)
#   4. Convert Mash distance → Jaccard, compare with ground truth
#
# Usage:
#   bash run_invindex_eval.sh <bench_data_dir> [dist_threshold] [threads]
#
# Example:
#   bash run_invindex_eval.sh /home/bigssd/zt/bench_data/L1M 0.5 64

set -euo pipefail

BENCH_DIR="${1:?usage: $0 <bench_data_dir> [dist_threshold] [threads]}"
DIST="${2:-0.5}"
THREADS="${3:-64}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Verify bench data ────────────────────────────────────────────────────────
[[ -f "$BENCH_DIR/pairs_meta.tsv" ]] || { echo "ERROR: $BENCH_DIR/pairs_meta.tsv not found"; exit 1; }

# ── Prepare combined file list (from actual .fa files, not stale A_list.txt) ─
ALL_LIST="$BENCH_DIR/all_list.txt"
find "$BENCH_DIR" -maxdepth 1 -name 'pair*_A.fa' -o -name 'pair*_B.fa' \
    | sort > "$ALL_LIST"
N_FILES=$(wc -l < "$ALL_LIST")
if [[ "$N_FILES" -eq 0 ]]; then
    echo "ERROR: no pair*_{A,B}.fa files found in $BENCH_DIR"
    exit 1
fi
echo "bench_dir: $BENCH_DIR"
echo "files: $N_FILES  dist_threshold: $DIST  threads: $THREADS"

# ── Tool table: name, executable, kmer_size, ground_truth_column ─────────────
#   pairs_meta.tsv columns (0-indexed):
#     0=pair_id  1=rate  2=len_bp  3=theo_j_k21  4=theo_j_k32  5=theo_j_k20
TOOLS=(
    "FastKMV:exe_test_FastKMV:21:3"
    "BinDash:exe_test_BinDash:21:3"
    "SetSketch:exe_test_SetSketch:32:4"
    "ProbMinHash:exe_test_ProbMinHash:21:3"
)

# Check for pairs_exact_j.tsv (more accurate ground truth for k=21)
EXACT_FILE="$BENCH_DIR/pairs_exact_j.tsv"
HAS_EXACT=0
if [[ -f "$EXACT_FILE" ]]; then
    HAS_EXACT=1
    echo "exact ground truth: $EXACT_FILE"
else
    echo "exact ground truth: not found (using theoretical Jaccard)"
fi

hms() { printf '%dh%02dm%02ds' $(($1/3600)) $(($1%3600/60)) $(($1%60)); }

RESULT_DIR="$BENCH_DIR/invindex_results"
mkdir -p "$RESULT_DIR"

# Collect results: "tool_name  k  time  global_mae"
declare -a RESULTS=()

# ── Run each tool and evaluate ───────────────────────────────────────────────
for entry in "${TOOLS[@]}"; do
    IFS=':' read -r NAME EXE KMER GT_COL <<< "$entry"
    EXE_PATH="$SCRIPT_DIR/$EXE"

    if [[ ! -x "$EXE_PATH" ]]; then
        echo "  [SKIP] $NAME: $EXE_PATH not found" >&2
        continue
    fi

    OUT_FILE="$RESULT_DIR/res.$NAME"

    T0=$SECONDS
    "$EXE_PATH" "$ALL_LIST" "$DIST" "$THREADS" "$OUT_FILE" \
        > "$RESULT_DIR/${NAME}_log.txt" 2>&1
    T1=$(( SECONDS - T0 ))

    # ── Pre-filter: grep matching A↔B pairs ──────────────────────────────
    FILTERED="$RESULT_DIR/${NAME}_pairs.tmp"
    grep -P 'pair(\d+)_A\.fa\t.*pair\1_B\.fa|pair(\d+)_B\.fa\t.*pair\2_A\.fa' \
         "$OUT_FILE" > "$FILTERED" 2>/dev/null || true

    # ── Compute Global MAE ───────────────────────────────────────────────
    MAE=$(python3 - "$FILTERED" "$BENCH_DIR/pairs_meta.tsv" "$KMER" "$GT_COL" \
                    "$HAS_EXACT" "$EXACT_FILE" <<'PYEOF'
import sys, re, math

out_file   = sys.argv[1]
meta_file  = sys.argv[2]
kmer       = int(sys.argv[3])
gt_col_idx = int(sys.argv[4])
has_exact  = sys.argv[5] == '1'
exact_file = sys.argv[6]

gt = {}
with open(meta_file) as f:
    f.readline()
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < gt_col_idx + 1: continue
        pid = int(cols[0]); rate = float(cols[1])
        gt[pid] = (rate, float(cols[gt_col_idx]))

if has_exact and kmer == 21:
    try:
        with open(exact_file) as f:
            f.readline()
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 4: continue
                pid = int(cols[0]); ej = float(cols[3])
                if pid in gt: gt[pid] = (gt[pid][0], ej)
    except: pass

pair_re = re.compile(r'pair(\d+)_([AB])\.fa')
found = {}
with open(out_file) as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 3: continue
        m1 = pair_re.search(cols[0]); m2 = pair_re.search(cols[1])
        if not (m1 and m2): continue
        if m1.group(1) == m2.group(1) and {m1.group(2), m2.group(2)} == {'A','B'}:
            d = float(cols[2]); p = math.exp(-kmer*d)
            found[int(m1.group(1))] = p/(2.0-p) if 2.0-p > 0 else 0.0

errs = []
for pid, (_, j_gt) in gt.items():
    errs.append(abs(found.get(pid, 0.0) - j_gt))
print(f"{sum(errs)/len(errs):.6f}" if errs else "N/A")
PYEOF
    )
    rm -f "$FILTERED"

    RESULTS+=("$(printf '%-14s k=%-4s %6ss  MAE=%s' "$NAME" "$KMER" "$T1" "$MAE")")
done

# ── Final summary ────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════"
echo "  Inverted-Index Accuracy Summary"
echo "  bench: $BENCH_DIR"
echo "  files: $N_FILES  threshold: $DIST"
echo "═══════════════════════════════════════"
for line in "${RESULTS[@]}"; do
    echo "  $line"
done
echo "═══════════════════════════════════════"
