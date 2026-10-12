#!/usr/bin/env bash
# run_pmh_bench.sh
# Run exe_test_ProbMinHash on each bench_data length,
# compute per-rate MAE and output PARSE_RATE lines to stdout.
#
# Usage:
#   bash run_pmh_bench.sh <bench_data_dir> [threads] [dist_threshold]
# Example:
#   bash run_pmh_bench.sh /home/bigssd/zt/bench_data 64 0.99

set -euo pipefail
BENCH="${1:?usage: $0 <bench_data_dir> [threads] [threshold]}"
THREADS="${2:-64}"
DIST="${3:-0.99}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/exe_test_ProbMinHash"
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found"; exit 1; }

LENGTHS=(L100k L500k L1M L2M L4M L8M)

for LEN in "${LENGTHS[@]}"; do
    DIR="$BENCH/$LEN"
    [[ -d "$DIR" ]] || { echo "SKIP $LEN (not found)" >&2; continue; }

    echo "=== $LEN ===" >&2

    # Build file list (all pair*_{A,B}.fa)
    LIST="$DIR/all_pmh.list"
    find "$DIR" -maxdepth 1 \( -name 'pair*_A.fa' -o -name 'pair*_B.fa' \) \
        | sort > "$LIST"
    [[ -s "$LIST" ]] || { echo "  no files in $LEN" >&2; continue; }

    OUT="$DIR/pmh_result.txt"
    "$EXE" "$LIST" "$DIST" "$THREADS" "$OUT" >/dev/null 2>&1 || true

    # Parse per-rate MAE
    python3 - "$DIR" "$OUT" <<'PYEOF'
import sys, re, math, os

bench_dir = sys.argv[1]
out_file  = sys.argv[2]
kmer      = 21

# Load ground truth (exact j or theoretical)
gt = {}  # pair_id -> (rate, len_bp, gt_j)
meta_file = os.path.join(bench_dir, "pairs_meta.tsv")
exact_file = os.path.join(bench_dir, "pairs_exact_j.tsv")
with open(meta_file) as f:
    f.readline()
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 4: continue
        pid = int(cols[0]); rate = float(cols[1]); lbp = int(cols[2]); j = float(cols[3])
        gt[pid] = [rate, lbp, j]
if os.path.exists(exact_file):
    with open(exact_file) as f:
        f.readline()
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 4: continue
            pid = int(cols[0])
            if pid in gt: gt[pid][2] = float(cols[3])

# Parse result file: dist format → jaccard
pair_re = re.compile(r'pair(\d+)_([AB])\.fa')
found = {}
try:
    with open(out_file) as f:
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3: continue
            m1 = pair_re.search(cols[0]); m2 = pair_re.search(cols[1])
            if not (m1 and m2): continue
            if m1.group(1) == m2.group(1) and {m1.group(2),m2.group(2)}=={'A','B'}:
                d = float(cols[2])
                # distance may be Mash distance or Jaccard depending on format
                # test_ProbMinHash outputs Mash distance
                try:
                    j_est = math.exp(-kmer * d) / (2.0 - math.exp(-kmer * d))
                except:
                    j_est = 0.0
                found[int(m1.group(1))] = j_est
except FileNotFoundError:
    pass

# Per-rate accumulation
from collections import defaultdict
rate_errs = defaultdict(list)
rep_len = 0
for pid, (rate, lbp, j_gt) in gt.items():
    if rep_len == 0: rep_len = lbp
    j_est = found.get(pid, 0.0)
    rate_errs[rate].append(abs(j_est - j_gt))

for rate in sorted(rate_errs):
    errs = rate_errs[rate]
    mae = sum(errs)/len(errs) if errs else 0
    print(f"PARSE_RATE len={rep_len} rate={rate:.4f} FreqPMH={mae:.6f}")
PYEOF
done
