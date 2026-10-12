#!/bin/bash
# Full SetSketch index consistency check (post-fix verification)
# Tests: --index == tile (direct MLE) across
#   (1) 9 alisim evolutionary models × 3 lengths × representative branch lengths
#   (2) 200 real bacterial genome pairs
#
# Each alisim condition: generate N_PAIRS pairs, run PathB (tile) and
# PathC (--index), compare pair-for-pair.

set -euo pipefail
cd "$(dirname "$0")"

IQTREE=/home/user_home/tongz/anaconda3/envs/cpp20/bin/iqtree3
RS=./rabbitsketch
THREADS=64
N_PAIRS=100    # pairs per alisim condition
K=21; BITS=13; A=5.0; BASE=2.0; MAXD=0.05

OUTDIR=results/verify_ss_full
mkdir -p "$OUTDIR"
SUMMARY="$OUTDIR/summary.tsv"
LOG="$OUTDIR/verify.log"
> "$LOG"
echo -e "condition\tn_B\tn_C\tmissed\trecall_pct\tmax_delta\tall_ok" > "$SUMMARY"

# ── 9 evolutionary models (same strings as run_alisim_eval.sh) ─────────────
declare -A MODEL_STR MODEL_INDEL
MODEL_STR[JC69]="JC"
MODEL_STR[HKY_ts2]="HKY{2.0}"
MODEL_STR[HKY_ts5]="HKY{5.0}"
MODEL_STR[GTR_ATrich]="GTR{1.5,4.0,1.2,0.8,3.5,1.0}+F{0.3,0.2,0.2,0.3}"
MODEL_STR[GTR_G02]="GTR+G{0.2}"
MODEL_STR[GTR_G05]="GTR+G{0.5}"
MODEL_STR[GTR_G20]="GTR+G{2.0}"
MODEL_STR[GTR_highGC]="GTR+G{1.0}+F{0.15,0.35,0.35,0.15}"
MODEL_STR[GTR_G_indel1]="GTR+G{1.0}"
MODEL_INDEL[JC69]=""
MODEL_INDEL[HKY_ts2]=""
MODEL_INDEL[HKY_ts5]=""
MODEL_INDEL[GTR_ATrich]=""
MODEL_INDEL[GTR_G02]=""
MODEL_INDEL[GTR_G05]=""
MODEL_INDEL[GTR_G20]=""
MODEL_INDEL[GTR_highGC]=""
MODEL_INDEL[GTR_G_indel1]="0.005,0.005"
MODELS=(JC69 HKY_ts2 HKY_ts5 GTR_ATrich GTR_G02 GTR_G05 GTR_G20 GTR_highGC GTR_G_indel1)

# Lengths: cover short/medium/long
LENGTHS=(10000 100000 1000000)
# Branch lengths that produce Mash distances mostly < MAXD
BRANCH_LENS=(0.005 0.01 0.025)

# ── comparison helper ─────────────────────────────────────────────────────────
compare_and_log() {
    local label="$1" pathB="$2" pathC="$3"
    python3 - "$label" "$pathB" "$pathC" "$SUMMARY" <<'PYEOF'
import sys, os
label, pathB, pathC, summary = sys.argv[1:]

def base(p): return os.path.basename(p.strip())
def load(f):
    d = {}
    try:
        with open(f) as fp:
            for l in fp:
                p = l.split()
                if len(p) < 3: continue
                d[tuple(sorted([base(p[0]), base(p[1])]))] = float(p[2])
    except FileNotFoundError:
        pass
    return d

dB = load(pathB); dC = load(pathC)
missed  = [k for k in dB if k not in dC]
extra   = [k for k in dC if k not in dB]
shared  = [abs(dB[k] - dC[k]) for k in dB if k in dC]
max_d   = max(shared) if shared else 0.0
recall  = 100.0 * len(shared) / max(1, len(dB))
all_ok  = (len(missed) == 0 and len(extra) == 0
           and (not shared or max_d < 1e-5))

flag = "✓" if all_ok else "✗"
print(f"  {flag} {label:50s}  B={len(dB):5,}  C={len(dC):5,}  "
      f"missed={len(missed):3d}  recall={recall:6.1f}%  maxΔ={max_d:.1e}")

with open(summary, 'a') as f:
    f.write(f"{label}\t{len(dB)}\t{len(dC)}\t{len(missed)}\t"
            f"{recall:.2f}\t{max_d:.2e}\t{'1' if all_ok else '0'}\n")
PYEOF
}

run_condition() {
    local label="$1" list_file="$2"
    local outB="$OUTDIR/${label}_B.dist"
    local outC="$OUTDIR/${label}_C.dist"
    $RS --setsketch         -i "$list_file" -l \
        -k $K -b $BITS -a $A -B $BASE -d $MAXD \
        -t $THREADS -o "$outB" 2>>"$LOG"
    $RS --setsketch --index -i "$list_file" -l \
        -k $K -b $BITS -a $A -B $BASE -d $MAXD \
        -t $THREADS -o "$outC" 2>>"$LOG"
    compare_and_log "$label" "$outB" "$outC"
    rm -f "$outB" "$outC"
}

# ─────────────────────────────────────────────────────────────────────────────
echo "==================================================================="
echo " SetSketch Consistency: --index vs tile  (all conditions)"
echo "==================================================================="
echo ""
echo "── Part 1: Alisim simulated (9 models × 3 L × 3 BL) ──────────────"
echo ""

for MODEL in "${MODELS[@]}"; do
    IQMOD="${MODEL_STR[$MODEL]}"
    INDEL="${MODEL_INDEL[$MODEL]}"

    for L in "${LENGTHS[@]}"; do
        for BL in "${BRANCH_LENS[@]}"; do
            LABEL="${MODEL}_L${L}_d${BL}"
            TMPDIR=$(mktemp -d)
            LIST="$TMPDIR/seqs.list"
            > "$LIST"

            printf "(A:%.6f,B:%.6f);\n" "$BL" "$BL" > "$TMPDIR/tree.nwk"

            # spawn N_PAIRS alisim jobs in parallel
            PIDS=()
            for ((i=0; i<N_PAIRS; i++)); do
                CMD=("$IQTREE" --alisim "$TMPDIR/p${i}_sim"
                     -m "$IQMOD" -t "$TMPDIR/tree.nwk"
                     --length "$L" --out-format fasta --seed $((42+i))
                     -nt 1 -quiet)
                [[ -n "$INDEL" ]] && CMD+=(--indel "$INDEL")
                "${CMD[@]}" >/dev/null 2>&1 &
                PIDS+=($!)
                if (( ${#PIDS[@]} >= THREADS )); then
                    wait "${PIDS[0]}" || true; PIDS=("${PIDS[@]:1}")
                fi
            done
            for p in "${PIDS[@]}"; do wait "$p" || true; done

            # split each multi-FASTA into A / B
            for ((i=0; i<N_PAIRS; i++)); do
                FA_RAW="$TMPDIR/p${i}_sim.fa"
                [[ -f "$FA_RAW" ]] || continue
                awk -v oa="$TMPDIR/pair_${i}_A.fa" \
                    -v ob="$TMPDIR/pair_${i}_B.fa" '
                    /^>/ { rec++; next }
                    rec==1 { seqA=seqA $0 }
                    rec==2 { seqB=seqB $0 }
                    END { printf ">A\n%s\n", seqA > oa
                          printf ">B\n%s\n", seqB > ob }
                ' "$FA_RAW"
                echo "$TMPDIR/pair_${i}_A.fa" >> "$LIST"
                echo "$TMPDIR/pair_${i}_B.fa" >> "$LIST"
                rm -f "$FA_RAW"
            done

            N_SEQ=$(wc -l < "$LIST")
            if (( N_SEQ < 4 )); then
                echo "  SKIP $LABEL (only $N_SEQ seq files generated)"
                rm -rf "$TMPDIR"
                continue
            fi
            run_condition "$LABEL" "$LIST"
            rm -rf "$TMPDIR"
        done
    done
done

echo ""
echo "── Part 2: 200 real bacterial genome pairs ─────────────────────────"
echo ""

PAIRS_TSV=real_genome_eval/pairs.tsv
REAL_LIST="$OUTDIR/real_files.list"
awk -F'\t' 'NR>1{print $5; print $6}' "$PAIRS_TSV" | sort -u > "$REAL_LIST"
echo "  Unique genomes: $(wc -l < "$REAL_LIST")"
run_condition "real_200pairs" "$REAL_LIST"

echo ""
echo "── Summary ─────────────────────────────────────────────────────────"
python3 << 'PYEOF2'
import csv

with open("results/verify_ss_full/summary.tsv") as f:
    rows = list(csv.DictReader(f, delimiter='\t'))

total  = len(rows)
passed = sum(1 for r in rows if r['all_ok'] == '1')
failed = [r for r in rows if r['all_ok'] != '1']
# conditions with actual pairs
have_pairs = [r for r in rows if int(r['n_B']) > 0]

print(f"  Total conditions        : {total}")
print(f"  Conditions with pairs   : {len(have_pairs)}")
print(f"  All-consistent ✓        : {passed}")
print(f"  Discrepant ✗            : {len(failed)}")

if failed:
    print(f"\n  FAILURES:")
    for r in failed:
        print(f"    {r['condition']:50s}  "
              f"missed={r['missed']:4s}  recall={r['recall_pct']}%  maxΔ={r['max_delta']}")
else:
    print(f"\n  ✓ SetSketch --index is fully consistent with tile (direct MLE)")
    print(f"    across ALL {len(have_pairs)} conditions that produced pairs.")
PYEOF2
