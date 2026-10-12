#!/bin/bash
# Verify that the three SetSketch computation paths give consistent results:
#   PATH A: single-pair mode (-i fileA fileB)  → direct jaccard_index_mle()
#   PATH B: all-pairs tile mode (-i list -l)   → two-stage incl-excl + MLE
#   PATH C: all-pairs index mode (--index)     → two-stage CSR index + MLE
#
# Uses -d 0.10 so we get a reasonably sized result set covering pairs
# near the 95%+ ANI range (Mash d ≤ 0.10).

set -euo pipefail
cd "$(dirname "$0")"

RABBITSKETCH=./rabbitsketch
N=20          # number of genomes to test (more overlap chances)
LIST=233.list
OUTDIR=results/verify_ss
K=21; BITS=13; A=5.0; BASE=2.0
THREADS=4
MAXD=0.10     # Mash distance cutoff

mkdir -p "$OUTDIR"

# ── Pick first N genomes ──────────────────────────────────────────────────
head -${N} "$LIST" > "$OUTDIR/test.list"
echo "Using $N genomes, max-dist=$MAXD"

# ── Path B: all-pairs tile (no index) ─────────────────────────────────────
echo ""
echo "=== PATH B: all-pairs tile (no --index) ==="
$RABBITSKETCH --setsketch -i "$OUTDIR/test.list" -l \
    -k $K -b $BITS -a $A -B $BASE -d $MAXD \
    -t $THREADS -o "$OUTDIR/pathB.dist" 2>"$OUTDIR/pathB.log"
echo "rows: $(wc -l < "$OUTDIR/pathB.dist")"

# ── Path C: all-pairs inverted index ──────────────────────────────────────
echo ""
echo "=== PATH C: all-pairs --index ==="
$RABBITSKETCH --setsketch --index -i "$OUTDIR/test.list" -l \
    -k $K -b $BITS -a $A -B $BASE -d $MAXD \
    -t $THREADS -o "$OUTDIR/pathC.dist" 2>"$OUTDIR/pathC.log"
echo "rows: $(wc -l < "$OUTDIR/pathC.dist")"

# ── Path A: single-pair loop (pairwise, write to a temp file per pair) ────
echo ""
echo "=== PATH A: single-pair mode (loop over all C(N,2) pairs) ==="
FILES=($(cat "$OUTDIR/test.list"))
TMPF=$(mktemp)
> "$OUTDIR/pathA.dist"
for ((i=0; i<${#FILES[@]}; i++)); do
    for ((j=i+1; j<${#FILES[@]}; j++)); do
        # write each pair to its own temp file to avoid /dev/stdout buffering
        $RABBITSKETCH --setsketch \
            -i "${FILES[$i]}" "${FILES[$j]}" \
            -k $K -b $BITS -a $A -B $BASE -d $MAXD \
            -t 1 -o "$TMPF" 2>/dev/null
        cat "$TMPF" >> "$OUTDIR/pathA.dist"
    done
done
rm -f "$TMPF"
echo "rows: $(wc -l < "$OUTDIR/pathA.dist")"

# ── Comparison via Python ──────────────────────────────────────────────────
echo ""
echo "=== Comparing paths ==="
python3 - "$OUTDIR" <<'PYEOF'
import sys, os

def basename(p):
    return os.path.basename(p.strip())

def load_dist(path):
    d = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if len(parts) < 3: continue
            a, b, dist = basename(parts[0]), basename(parts[1]), float(parts[2])
            key = tuple(sorted([a, b]))
            d[key] = dist
    return d

BASE = sys.argv[1]
dA = load_dist(f"{BASE}/pathA.dist")
dB = load_dist(f"{BASE}/pathB.dist")
dC = load_dist(f"{BASE}/pathC.dist")

all_keys = set(dA) | set(dB) | set(dC)
print(f"Pair counts:  A={len(dA)}  B={len(dB)}  C={len(dC)}  union={len(all_keys)}")

def compare(ref, tgt, name_ref, name_tgt):
    shared_diffs = []
    only_ref, only_tgt = [], []
    for k in sorted(all_keys):
        vr = ref.get(k)
        vt = tgt.get(k)
        if vr is None and vt is None: continue
        if vr is None: only_tgt.append((k, vt)); continue
        if vt is None: only_ref.append((k, vr)); continue
        shared_diffs.append((abs(vr - vt), k, vr, vt))

    print(f"\n  {name_ref}  vs  {name_tgt}")
    print(f"  shared pairs: {len(shared_diffs)}")
    if shared_diffs:
        diffs = [x[0] for x in shared_diffs]
        print(f"    max  |Δdist| = {max(diffs):.8f}")
        print(f"    mean |Δdist| = {sum(diffs)/len(diffs):.8f}")
        tols = [1e-6, 1e-5, 1e-4]
        for tol in tols:
            n_ok = sum(1 for d in diffs if d < tol)
            print(f"    match |Δ|<{tol:.0e}: {n_ok}/{len(diffs)}", "✓ ALL" if n_ok==len(diffs) else "")
        if max(diffs) > 1e-4:
            print(f"    top-5 discordant pairs:")
            for diff, (fa, fb), vr, vt in sorted(shared_diffs, reverse=True)[:5]:
                print(f"      {fa} vs {fb}")
                print(f"        {name_ref}={vr:.8f}  {name_tgt}={vt:.8f}  Δ={diff:.8f}")
    if only_ref:
        print(f"  pairs ONLY in {name_ref}: {len(only_ref)}")
        for (fa,fb),v in only_ref[:3]:
            print(f"    {fa} vs {fb}  dist={v:.6f}")
    if only_tgt:
        print(f"  pairs ONLY in {name_tgt}: {len(only_tgt)}")
        for (fa,fb),v in only_tgt[:3]:
            print(f"    {fa} vs {fb}  dist={v:.6f}")

compare(dA, dB, "PathA(single-pair MLE)", "PathB(tile two-stage)")
compare(dA, dC, "PathA(single-pair MLE)", "PathC(index two-stage)")
compare(dB, dC, "PathB(tile two-stage)", "PathC(index two-stage)")
PYEOF
