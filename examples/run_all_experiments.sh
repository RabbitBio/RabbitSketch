#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_all_experiments.sh  –  master experiment runner
#
# Runs the complete evaluation suite in order:
#
#   Phase 1 (in-memory, no disk space needed):
#     A. Parameter sensitivity sweep (all algos, dense grid, L=1M)
#     B. Parameter × sequence-length cross sweep
#     C. Equal-memory fair comparison
#     D. Low-complexity contamination sweep
#     E. Containment evaluation (length-ratio sweep)
#     F. Existing length sensitivity sweep (exe_eval_sketch_acc)
#
#   Phase 2 (requires iqtree2 + disk space in /home/bigssd/zt):
#     G. AliSim model evaluation (11 models × 5 lengths × 8 divergences)
#        Space-managed: generates & evaluates & deletes FA per combo.
#
# All results stored under  examples/results/
#
# Usage:
#   bash run_all_experiments.sh [THREADS] [PAIRS]
#   THREADS  (default: 64)
#   PAIRS    pairs per rate for param sweeps (default: 300)
#             alisim pairs per combo         (default: 1000)
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"

THREADS="${1:-64}"
PAIRS="${2:-300}"
ALISIM_PAIRS=1000
RESDIR="results"
LOG="$RESDIR/run_all.log"

mkdir -p "$RESDIR"
exec > >(tee -a "$LOG") 2>&1   # mirror all output to log file

stamp() { date '+[%Y-%m-%d %H:%M:%S]'; }

banner() {
    echo ""
    echo "╔══════════════════════════════════════════════════════════════╗"
    printf  "║  %-60s ║\n" "$1"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo "$(stamp) $1"
}

echo "════════════════════════════════════════════════════════════════"
echo " RabbitSketch Comprehensive Evaluation"
echo " THREADS=$THREADS  PAIRS=$PAIRS  ALISIM_PAIRS=$ALISIM_PAIRS"
echo " Results: $RESDIR/"
echo " Log:     $LOG"
echo " Start:   $(date -Iseconds)"
echo "════════════════════════════════════════════════════════════════"

# ── prerequisite checks ───────────────────────────────────────────────────────
for exe in ./exe_eval_param_sweep ./exe_eval_sketch_acc; do
    if [[ ! -x "$exe" ]]; then
        echo "[ERROR] $exe not found. Run 'make' first."
        exit 1
    fi
done

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 1-A: Parameter sensitivity (Sections A–D in run_param_sweep_all.sh)
# ═══════════════════════════════════════════════════════════════════════════════
banner "Phase 1-A: Parameter sensitivity sweep (A/B/C/D)"
bash run_param_sweep_all.sh "$THREADS" "$PAIRS" 1000000
echo "$(stamp) Phase 1-A done"
echo "  → $RESDIR/param_sweep.csv"
echo "  → $RESDIR/crosslen_sweep.csv"
echo "  → $RESDIR/equal_memory.csv"
echo "  → $RESDIR/lowcomp_sweep.csv"

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 1-E: Containment evaluation
# ═══════════════════════════════════════════════════════════════════════════════
banner "Phase 1-E: Containment sweep"
bash run_containment_sweep.sh "$THREADS" "$PAIRS"
echo "$(stamp) Phase 1-E done"
echo "  → $RESDIR/containment_ratio.csv"
echo "  → $RESDIR/containment_size.csv"
echo "  → $RESDIR/containment_rate.csv"

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 1-F: Length sensitivity (existing exe_eval_sketch_acc sweep)
# ═══════════════════════════════════════════════════════════════════════════════
banner "Phase 1-F: Sequence length sensitivity (all algorithms)"
bash run_eval_accuracy_sweep.sh "$THREADS" "$PAIRS" "$RESDIR/eval_accuracy_sweep.log"
echo "$(stamp) Phase 1-F done"
echo "  → $RESDIR/eval_accuracy_sweep.log  (PARSE_FIXED/PARSE_VARLEN lines)"
# Parse sweep log → CSV
python3 - <<'PYEOF' "$RESDIR/eval_accuracy_sweep.log" "$RESDIR/length_sweep.csv"
import sys, re
log, out = sys.argv[1], sys.argv[2]
rows = []
hdr = "len_bp,rate,HLL,SetSketch,KSSD,MinHash,ProbMH,OnePerm,KMV,BinDash"
pat = re.compile(
    r'PARSE_(?:FIXED|RATE)\s+len=(\d+)\s+rate=([0-9.]+)\s+'
    r'HLL=([0-9.]+)\s+SetSketch=([0-9.]+)\s+KSSD=([0-9.]+)\s+'
    r'MinHash=([0-9.]+)\s+ProbMH=([0-9.]+)\s+OnePerm=([0-9.]+)\s+'
    r'KMV=([0-9.]+)\s+BinDash=([0-9.]+)')
with open(log) as f:
    for line in f:
        m = pat.search(line)
        if m: rows.append(','.join(m.groups()))
with open(out, 'w') as f:
    f.write(hdr + '\n')
    f.write('\n'.join(rows) + '\n')
print(f"Parsed {len(rows)} rows → {out}")
PYEOF
echo "$(stamp) Phase 1-F CSV parsed → $RESDIR/length_sweep.csv"

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 2-G: AliSim model evaluation
# ═══════════════════════════════════════════════════════════════════════════════
banner "Phase 2-G: AliSim evolutionary model evaluation"

# Check iqtree2
IQTREE_OK=0
for cand in iqtree2 iqtree; do
    command -v "$cand" &>/dev/null && IQTREE_OK=1 && break
done
if (( !IQTREE_OK )); then
    _it=$(conda run -n cpp20 which iqtree2 2>/dev/null || true)
    [[ -n "$_it" ]] && IQTREE_OK=1
fi

if (( IQTREE_OK )); then
    bash run_alisim_eval.sh "$THREADS" "$ALISIM_PAIRS"
    echo "$(stamp) Phase 2-G done"
    echo "  → $RESDIR/alisim_eval.csv"
else
    echo "$(stamp) [SKIP] iqtree2 not found – Phase 2-G skipped"
    echo "         Install: conda install -c bioconda -c conda-forge iqtree"
    echo "         Then re-run: bash run_alisim_eval.sh $THREADS $ALISIM_PAIRS"
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════════"
echo " ALL EXPERIMENTS COMPLETE"
echo " $(stamp)"
echo ""
echo " Results directory: $RESDIR/"
ls -lh "$RESDIR"/*.csv "$RESDIR"/*.log 2>/dev/null | awk '{print "   "$NF"  "$5}'
echo "════════════════════════════════════════════════════════════════"
