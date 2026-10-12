# Shared AliSim branch-length ↔ directory tag mapping.
# Sourced by gen_alisim_data.sh, run_alisim_eval.sh, run_setsketch_mle_eval.sh.
#
# Tree: (A:d, B:d)  →  total divergence ≈ 2d per site.
# Mash distance (k=21):  M = -ln(2J/(1+J)) / 21

# Original AliSim grid (8 points) — MLE backfill target before Mash supplement.
ORIGINAL_BRANCH_LENS=(
    0.0005
    0.0025
    0.005
    0.01
    0.025
    0.05
    0.10
    0.15
)

# Full grid: original 8 + 9 Mash-extension points (17 total, sorted).
BRANCH_LENS_ALL=(
    0.0005
    0.0025
    0.005
    0.008
    0.01
    0.025
    0.045
    0.05
    0.053
    0.10
    0.11
    0.13
    0.15
    0.20
    0.28
    0.38
    0.48
)

# Supplement-only: new points (target Mash ≈ 0.02, 0.08, 0.10†, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60)
SUPPLEMENT_BRANCH_LENS=(
    0.008
    0.045
    0.053
    0.11
    0.13
    0.20
    0.28
    0.38
    0.48
)

# Mash calibration add-on (all 7 algos via run_mash_calib_supplement.sh)
# Chosen from J→M calibration on L10k; verify empirical M after run.
MASH_CALIB_BRANCH_LENS=(
    0.038
    0.080
    0.116
    0.31
)
MASH_CALIB_TARGET_M=(0.075 0.150 0.200 0.400)

declare -A BL_TAG
BL_TAG[0.0005]="d0001"
BL_TAG[0.0025]="d0005"
BL_TAG[0.005]="d001"
BL_TAG[0.008]="d0008"
BL_TAG[0.01]="d002"
BL_TAG[0.025]="d005"
BL_TAG[0.045]="d0045"
BL_TAG[0.05]="d010"
BL_TAG[0.053]="d0053"
BL_TAG[0.10]="d020"
BL_TAG[0.11]="d011"
BL_TAG[0.13]="d013"
BL_TAG[0.15]="d030"
BL_TAG[0.20]="d040"
BL_TAG[0.28]="d056"
BL_TAG[0.38]="d076"
BL_TAG[0.48]="d096"
BL_TAG[0.038]="d0038"
BL_TAG[0.080]="d0080"
BL_TAG[0.116]="d0116"
BL_TAG[0.31]="d0310"

# Default alias used by main eval scripts
BRANCH_LENS=("${BRANCH_LENS_ALL[@]}")

# Varlen experiment: original + extension representatives
VARLEN_BLS=(0.005 0.05 0.15)
VARLEN_BLS_SUPPLEMENT=(0.008 0.13 0.48)

# Target Mash (k=21) for documentation / plotting (supplement branch lengths)
declare -A BL_TARGET_MASH
BL_TARGET_MASH[0.008]="0.02"
BL_TARGET_MASH[0.045]="0.08"
BL_TARGET_MASH[0.053]="0.10"
BL_TARGET_MASH[0.11]="0.15"
BL_TARGET_MASH[0.13]="0.20"
BL_TARGET_MASH[0.20]="0.30"
BL_TARGET_MASH[0.28]="0.40"
BL_TARGET_MASH[0.38]="0.50"
BL_TARGET_MASH[0.48]="0.60"
BL_TARGET_MASH[0.038]="0.075"
BL_TARGET_MASH[0.080]="0.150"
BL_TARGET_MASH[0.116]="0.200"
BL_TARGET_MASH[0.31]="0.400"

# Landscape column axis: design Mash → AliSim BL (x tick = Mash distance only)
LANDSCAPE_MASH_BL=(
    "0.001:0.0005"
    "0.005:0.0025"
    "0.01:0.005"
    "0.05:0.025"
    "0.075:0.038"
    "0.10:0.053"
    "0.15:0.080"
    "0.20:0.116"
    "0.25:0.150"
    "0.30:0.200"
    "0.35:0.280"
    "0.40:0.31"
)

bl_tag_for() {
    local bl="$1"
    echo "${BL_TAG[$bl]:-}"
}
