#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_mash_dist_supplement.sh
#
# Supplement AliSim evaluation for extended Mash-distance coverage (k=21).
# Adds 9 new branch lengths without re-running the original 8 BLs:
#
#   BL      tag     ~target Mash
#   0.008   d0008   0.02
#   0.045   d0045   0.08
#   0.053   d0053   0.10
#   0.11    d011    0.15
#   0.13    d013    0.20
#   0.20    d040    0.30
#   0.28    d056    0.40
#   0.38    d076    0.50
#   0.48    d096    0.60
#
# Fixed: 11 models × 9 lengths × 9 BLs = 891 combos
# Varlen: 3 models × 9 lengths × 3 BLs = 81 combos
#
# Appends rows to results/alisim_sweep.csv (same format as run_alisim_eval.sh).
#
# Usage:
#   bash run_mash_dist_supplement.sh [THREADS] [PAIRS_MAX]
#   nohup bash run_mash_dist_supplement.sh 64 1000 > results/mash_dist_supplement.log 2>&1 &
#
# After completion:
#   python3 build_master_csv.py
#   bash run_post_alisim.sh   # optional full figure rebuild
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"

# IQ-TREE (AliSim) — ensure cpp20 bin is on PATH (nohup shells often lack conda init)
for _iqbin in "${HOME}/anaconda3/envs/cpp20/bin" \
               "/home/user_home/tongz/anaconda3/envs/cpp20/bin"; do
    if [[ -d "$_iqbin" ]]; then
        export PATH="$_iqbin:$PATH"
        break
    fi
done

THREADS="${1:-64}"
PAIRS_MAX="${2:-1000}"
LOG="${LOG:-results/mash_dist_supplement.log}"
COMBO_LIST="results/.mash_dist_supplement_combos.tsv"

# shellcheck source=alisim_bl_config.sh
source "$(dirname "$0")/alisim_bl_config.sh"

echo "══════════════════════════════════════════════════════════════════"
echo " Mash-distance supplement"
echo " New branch lengths: ${SUPPLEMENT_BRANCH_LENS[*]}"
echo " THREADS=$THREADS  PAIRS_MAX=$PAIRS_MAX"
echo " Log: $LOG"
echo "══════════════════════════════════════════════════════════════════"

# Write combo checklist (for progress / resume debugging)
{
    echo -e "dataset\tbranch_len\tbl_tag\ttarget_mash"
    for bl in "${SUPPLEMENT_BRANCH_LENS[@]}"; do
        tag="${BL_TAG[$bl]}"
        tm="${BL_TARGET_MASH[$bl]:-}"
        for model in JC69 HKY_ts2 HKY_ts5 GTR_ATrich GTR_G02 GTR_G05 GTR_G20 GTR_highGC \
                     GTR_G_indel1 GTR_G_indel2 GTR_G_indel3; do
            for lt in L1k L10k L50k L100k L500k L1M L2M L4M L8M; do
                echo -e "${model}/fixed_${lt}/${tag}\t${bl}\t${tag}\t${tm}"
            done
        done
    done
} > "$COMBO_LIST"
echo "Wrote combo list: $COMBO_LIST ($(wc -l < "$COMBO_LIST") lines)"

export MASH_SUPPLEMENT_ONLY=1
mkdir -p results
bash run_alisim_eval.sh "$THREADS" "$PAIRS_MAX"
