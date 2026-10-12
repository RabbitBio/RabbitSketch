#!/usr/bin/env python3
"""
Evaluate all 7 rabbitsketch algorithms on bench_data using pairwise mode.
Uses the actual rabbitsketch binary: rabbitsketch --<algo> -d 1.0 -i A B
Outputs PARSE_RATE lines for downstream analysis.
"""
import subprocess, os, re, math, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

RABBITSKETCH = os.path.join(os.path.dirname(__file__), "rabbitsketch")
BENCH_BASE   = "/home/bigssd/zt/bench_data"
OUT_DIR      = "/home/bigssd/zt/bench_data/pairwise_eval_results"
LENGTHS      = ["L100k", "L500k", "L1M", "L2M", "L4M", "L8M"]
ALGOS        = ["minhash", "kssd", "hll", "setsketch", "probminhash", "fastkmv", "bindash"]
DISPLAY      = {
    "minhash":"MinHash","kssd":"KSSD","hll":"HLL","setsketch":"SetSketch",
    "probminhash":"PMH-freq","fastkmv":"FastKMV","bindash":"BinDash"
}
K = 21
WORKERS = 32

os.makedirs(OUT_DIR, exist_ok=True)

def mash_to_jaccard(d):
    """Mash distance -> Jaccard: J = exp(-dk) / (2 - exp(-dk))"""
    r = math.exp(-d * K)
    if r >= 2.0: return 1.0
    return r / (2.0 - r)

def run_pair(rs_bin, algo, fa, fb, tmp_dir):
    """Run rabbitsketch pairwise, return Jaccard estimate or None on error."""
    out = os.path.join(tmp_dir, f"p_{os.getpid()}.dist")
    try:
        cmd = [rs_bin, f"--{algo}", "-d", "1.0", "-i", fa, fb, "-o", out]
        res = subprocess.run(cmd, capture_output=True, timeout=30)
        if res.returncode != 0:
            return None
        if os.path.exists(out):
            with open(out) as f:
                line = f.readline().strip()
            if line:
                parts = line.split()
                if len(parts) >= 3:
                    return mash_to_jaccard(float(parts[2]))
        return 0.0
    except Exception:
        return None
    finally:
        try: os.remove(out)
        except: pass

def eval_length(bench_dir, algo, log_fh):
    """Evaluate one algorithm on one length directory. Returns per-rate MAE."""
    # Load pairs metadata
    meta_file = os.path.join(bench_dir, "pairs_meta.tsv")
    exact_file = os.path.join(bench_dir, "pairs_exact_j.tsv")
    if not os.path.exists(meta_file):
        return {}

    pairs = []  # (pair_id, rate, len_bp, theo_j)
    with open(meta_file) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4: continue
            pairs.append((int(parts[0]), float(parts[1]), int(parts[2]), float(parts[3])))

    # Override with exact J if available
    exact_j = {}
    if os.path.exists(exact_file):
        with open(exact_file) as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    exact_j[int(parts[0])] = float(parts[3])

    gt = {pid: exact_j.get(pid, theo) for pid, rate, lbp, theo in pairs}

    # Build file paths
    def fa_path(pid): return os.path.join(bench_dir, f"pair{pid:06d}_A.fa")
    def fb_path(pid): return os.path.join(bench_dir, f"pair{pid:06d}_B.fa")

    tmp_dir = "/tmp/rs_paireval"
    os.makedirs(tmp_dir, exist_ok=True)

    results = {}  # pair_id -> j_est
    futures = {}
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for pid, rate, lbp, theo in pairs:
            fa, fb = fa_path(pid), fb_path(pid)
            if not os.path.exists(fa) or not os.path.exists(fb):
                continue
            fut = ex.submit(run_pair, RABBITSKETCH, algo, fa, fb, tmp_dir)
            futures[fut] = pid
        for fut in as_completed(futures):
            pid = futures[fut]
            j = fut.result()
            results[pid] = j if j is not None else 0.0

    # Aggregate per rate
    rate_data = defaultdict(list)
    lbp_ref = pairs[0][2] if pairs else 0
    for pid, rate, lbp, theo in pairs:
        if pid not in results: continue
        j_est = results[pid]
        j_gt  = gt[pid]
        rate_data[rate].append(abs(j_est - j_gt))

    return lbp_ref, rate_data

def main():
    import time
    print(f"Evaluating {len(ALGOS)} algorithms × {len(LENGTHS)} lengths...\n", flush=True)

    # Collect all PARSE_RATE lines
    all_lines = []

    for algo in ALGOS:
        disp = DISPLAY[algo]
        for lname in LENGTHS:
            bench_dir = os.path.join(BENCH_BASE, lname)
            if not os.path.isdir(bench_dir):
                print(f"  SKIP {lname} (not found)", flush=True)
                continue
            t0 = time.time()
            print(f"  Running {disp:12s} {lname} ... ", end="", flush=True)
            lbp, rate_data = eval_length(bench_dir, algo, sys.stdout)
            elapsed = time.time() - t0
            if not rate_data:
                print(f"NO DATA", flush=True)
                continue

            for rate in sorted(rate_data):
                errs = rate_data[rate]
                mae  = sum(errs)/len(errs) if errs else float('nan')
                line = f"PARSE_RATE len={lbp} rate={rate:.4f} {disp}={mae:.6f}"
                all_lines.append(line)
                print(f"\n    rate={rate:.4f} n={len(errs)} MAE={mae:.5f}", end="", flush=True)
            print(f"\n  done in {elapsed:.1f}s\n", flush=True)

    # Write combined log
    log_path = os.path.join(OUT_DIR, "pairwise_combined.log")
    with open(log_path, "w") as f:
        for line in all_lines:
            f.write(line + "\n")
    print(f"\nLog written to: {log_path}", flush=True)

if __name__ == "__main__":
    main()
