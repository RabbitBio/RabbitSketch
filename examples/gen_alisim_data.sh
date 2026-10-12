#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# gen_alisim_data.sh
#
# Generate pairwise sequence datasets with IQ-TREE 2 AliSim under multiple
# evolutionary models, sequence lengths, and divergence levels.
#
# Output layout:
#   /home/bigssd/zt/ref/
#     {MODEL}/{LEN_TAG}/{RATE_TAG}/
#       pair_{i}_A.fa   (sequence A)
#       pair_{i}_B.fa   (sequence B)
#       pairs_meta.tsv  (pair_id, rate, len_A, len_B)
#
# Usage:
#   bash gen_alisim_data.sh [THREADS] [PAIRS_PER_COMBO]
#
# Prerequisites:
#   conda install -c bioconda -c conda-forge iqtree   (provides iqtree2)
#   or: mamba install -c bioconda iqtree
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

THREADS="${1:-8}"
PAIRS="${2:-100}"
OUTBASE="/home/bigssd/zt/ref"

# ── locate iqtree2 ─────────────────────────────────────────────────────────
IQTREE=""
for candidate in iqtree2 iqtree "$(conda run -n cpp20 which iqtree2 2>/dev/null)" "$(which iqtree2 2>/dev/null)"; do
    if command -v "$candidate" &>/dev/null 2>&1; then
        IQTREE="$candidate"; break
    fi
done
if [[ -z "$IQTREE" ]]; then
    echo "[ERROR] iqtree2 not found."
    echo "Install with: conda install -c bioconda -c conda-forge iqtree"
    echo "Then re-run this script."
    exit 1
fi
echo "Using IQ-TREE: $IQTREE  ($(${IQTREE} --version 2>&1 | head -1))"

# ── parameters ─────────────────────────────────────────────────────────────
# Evolutionary models
# Pure substitution models
#
# Coverage rationale:
#   JC69           : baseline, equal rates & frequencies (control)
#   HKY_ts2/ts5    : transition/transversion bias (κ=2 human-like, κ=5 mitochondria-like)
#   GTR            : fully general substitution, AT-rich composition (~40% GC)
#   GTR_G05        : rate heterogeneity α=0.5 – highly variable (16S rRNA-like)
#   GTR_G20        : rate heterogeneity α=2.0 – mildly variable
#   GTR_G02        : extreme rate heterogeneity α=0.2 – conserved core genomes [NEW]
#   GTR_highGC     : high-GC composition F≈70% GC (Streptomyces / M.tuberculosis) [NEW]
declare -A MODELS
MODELS["JC69"]="JC"
MODELS["HKY_ts2"]="HKY{2.0}"
MODELS["HKY_ts5"]="HKY{5.0}"
MODELS["GTR"]="GTR{1.5,4.0,1.2,0.8,3.5,1.0}+F{0.3,0.2,0.2,0.3}"
MODELS["GTR_G05"]="GTR+G{0.5}"
MODELS["GTR_G20"]="GTR+G{2.0}"
MODELS["GTR_G02"]="GTR+G{0.2}"
MODELS["GTR_highGC"]="GTR+G{1.0}+F{0.15,0.35,0.35,0.15}"

# Indel models: GTR+G with insertion and deletion rates
# These test how k-mer sketches degrade under length-changing mutations.
# Format for alisim: model + --indel ins_rate del_rate
declare -A INDEL_MODELS
declare -A INDEL_RATES
INDEL_MODELS["GTR_G_indel1"]="GTR+G{1.0}"
INDEL_RATES["GTR_G_indel1"]="0.005 0.005"   # 0.5% ins + 0.5% del
INDEL_MODELS["GTR_G_indel2"]="GTR+G{1.0}"
INDEL_RATES["GTR_G_indel2"]="0.01 0.01"     # 1% ins + 1% del
INDEL_MODELS["GTR_G_indel3"]="GTR+G{1.0}"
INDEL_RATES["GTR_G_indel3"]="0.03 0.03"     # 3% ins + 3% del

# Sequence lengths (bp)
declare -a LENGTHS=(5000 50000 500000 1000000 4000000)
declare -A LEN_TAGS
LEN_TAGS[5000]="L5k"
LEN_TAGS[50000]="L50k"
LEN_TAGS[500000]="L500k"
LEN_TAGS[1000000]="L1M"
LEN_TAGS[4000000]="L4M"

# Branch lengths — see alisim_bl_config.sh (17 levels incl. Mash extension)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=alisim_bl_config.sh
source "$SCRIPT_DIR/alisim_bl_config.sh"
declare -A BL_TAGS=()
for bl in "${BRANCH_LENS[@]}"; do
    BL_TAGS[$bl]="${BL_TAG[$bl]}"
done

TMPDIR_BASE="/tmp/alisim_tmp_$$"
mkdir -p "$TMPDIR_BASE"
trap 'rm -rf "$TMPDIR_BASE"' EXIT

# ── helper: generate one batch of PAIRS pairs ──────────────────────────────
# Args: model_name iqtree_model len bl_tag branch_len outdir
gen_batch() {
    local model_name="$1"
    local iqtree_model="$2"
    local len="$3"
    local bl_tag="$4"
    local bl="$5"
    local outdir="$6"
    local indel_rates="${7:-}"   # optional: "ins_rate del_rate" string

    mkdir -p "$outdir"
    local meta_file="$outdir/pairs_meta.tsv"
    echo -e "pair_id\tbranch_len\tlen_A\tlen_B\tmodel" > "$meta_file"

    local tmpd="$TMPDIR_BASE/${model_name}_${len}_${bl_tag}"
    mkdir -p "$tmpd"

    # Write the 2-taxon Newick tree
    local tree_file="$tmpd/tree2.nwk"
    printf "(A:%.6f,B:%.6f);\n" "$bl" "$bl" > "$tree_file"

    # Generate PAIRS simulation runs in parallel
    local pids=()
    local completed=0
    local slot=0

    for ((i=0; i<PAIRS; i++)); do
        local pfx="$tmpd/pair_${i}"

        # Build alisim command (add --indel if rates provided)
        local alisim_cmd=("$IQTREE"
            --alisim "${pfx}_sim"
            -m "$iqtree_model"
            -t "$tree_file"
            --length "$len"
            --out-format fasta
            --seed $((42 + i))
            -nt 1
            -quiet)
        if [[ -n "$indel_rates" ]]; then
            read -r ins_rate del_rate <<< "$indel_rates"
            alisim_cmd+=(--indel "$ins_rate" "$del_rate")
        fi

        # Run alisim in background
        "${alisim_cmd[@]}" >/dev/null 2>&1 &

        pids+=($!)
        ((slot++))

        # Throttle: wait when THREADS jobs are running
        if (( slot >= THREADS )); then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
            ((slot--))
        fi
    done
    # Wait remaining
    for pid in "${pids[@]}"; do wait "$pid"; done

    # Split the output FASTA (alisim writes all taxa into one .fa file)
    # and copy to final output directory
    for ((i=0; i<PAIRS; i++)); do
        local sim_fa="$tmpd/pair_${i}_sim.fa"
        if [[ ! -f "$sim_fa" ]]; then
            echo "[WARN] Missing $sim_fa" >&2
            continue
        fi

        # Extract sequence A (first record) and B (second record)
        awk '
            /^>/ { rec++; next }
            rec==1 { seqA = seqA $0 }
            rec==2 { seqB = seqB $0 }
            END {
                printf ">A\n%s\n", seqA > outA
                printf ">B\n%s\n", seqB > outB
            }
        ' outA="$outdir/pair_${i}_A.fa" outB="$outdir/pair_${i}_B.fa" "$sim_fa"

        local lenA lenB
        lenA=$(grep -v "^>" "$outdir/pair_${i}_A.fa" | tr -d '\n' | wc -c)
        lenB=$(grep -v "^>" "$outdir/pair_${i}_B.fa" | tr -d '\n' | wc -c)
        printf "%d\t%.6f\t%d\t%d\t%s\n" "$i" "$bl" "$lenA" "$lenB" "$model_name" >> "$meta_file"
    done

    rm -rf "$tmpd"
    echo "  [OK] $outdir  ($PAIRS pairs)"
}

# ── main loop ─────────────────────────────────────────────────────────────
echo "=== AliSim data generation ==="
echo "Output base: $OUTBASE"
echo "Models: ${!MODELS[*]}"
echo "Lengths: ${LENGTHS[*]}"
echo "Branch lengths: ${BRANCH_LENS[*]}"
echo "Pairs per combo: $PAIRS"
echo "Threads: $THREADS"
echo "Start: $(date -Iseconds)"
echo ""

total=0

# ── pure substitution models ─────────────────────────────────────────────────
for model_name in "${!MODELS[@]}"; do
    iqtree_model="${MODELS[$model_name]}"
    for len in "${LENGTHS[@]}"; do
        len_tag="${LEN_TAGS[$len]}"
        for bl in "${BRANCH_LENS[@]}"; do
            bl_tag="${BL_TAGS[$bl]}"
            outdir="$OUTBASE/${model_name}/${len_tag}/${bl_tag}"
            if [[ -f "$outdir/pairs_meta.tsv" ]]; then
                echo "  [SKIP] $model_name/$len_tag/$bl_tag"
                continue
            fi
            printf "Generating: %-22s %-8s %-8s ... " "$model_name" "$len_tag" "$bl_tag"
            gen_batch "$model_name" "$iqtree_model" "$len" "$bl_tag" "$bl" "$outdir"
            ((total++))
        done
    done
done

# ── indel models ─────────────────────────────────────────────────────────────
echo ""
echo "=== Generating indel models ==="
for model_name in "${!INDEL_MODELS[@]}"; do
    iqtree_model="${INDEL_MODELS[$model_name]}"
    indel_rates="${INDEL_RATES[$model_name]}"
    for len in "${LENGTHS[@]}"; do
        len_tag="${LEN_TAGS[$len]}"
        for bl in "${BRANCH_LENS[@]}"; do
            bl_tag="${BL_TAGS[$bl]}"
            outdir="$OUTBASE/${model_name}/${len_tag}/${bl_tag}"
            if [[ -f "$outdir/pairs_meta.tsv" ]]; then
                echo "  [SKIP] $model_name/$len_tag/$bl_tag"
                continue
            fi
            printf "Generating: %-22s %-8s %-8s ... " "$model_name" "$len_tag" "$bl_tag"
            gen_batch "$model_name" "$iqtree_model" "$len" "$bl_tag" "$bl" "$outdir" "$indel_rates"
            ((total++))
        done
    done
done

echo ""
echo "Done: $(date -Iseconds)"
echo "Total combos generated: $total"
echo "Data in: $OUTBASE"
