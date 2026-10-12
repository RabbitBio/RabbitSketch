#!/usr/bin/env python3
"""
Fast parallel FASTA splitter – replaces the _split_and_meta bash loop.

Each IQ-TREE AliSim output file (p{i}_sim.fa) contains two sequences (A & B).
This script reads all of them in parallel, writes pair_{i}_A.fa / pair_{i}_B.fa,
and emits pairs_meta.tsv.

Usage:
    python3 split_pairs_fast.py <tmpd> <outdir> <branch_len> <n_pairs> [n_workers]
"""
import sys
import os
from concurrent.futures import ThreadPoolExecutor, as_completed


def split_one(i: int, fa: str, oa: str, ob: str, bl: float):
    if not os.path.exists(fa):
        return None
    seqA: list[str] = []
    seqB: list[str] = []
    rec = 0
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                rec += 1
            elif rec == 1:
                seqA.append(line.rstrip("\r\n"))
            elif rec == 2:
                seqB.append(line.rstrip("\r\n"))
    sA = "".join(seqA)
    sB = "".join(seqB)
    with open(oa, "w") as fh:
        fh.write(f">A\n{sA}\n")
    with open(ob, "w") as fh:
        fh.write(f">B\n{sB}\n")
    return (i, bl, len(sA), len(sB))


def main():
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} <tmpd> <outdir> <branch_len> <n_pairs> [n_workers]",
              file=sys.stderr)
        sys.exit(1)

    tmpd = sys.argv[1]
    outdir = sys.argv[2]
    bl = float(sys.argv[3])
    n_pairs = int(sys.argv[4])
    n_workers = int(sys.argv[5]) if len(sys.argv) > 5 else 16

    tasks = [
        (i,
         os.path.join(tmpd, f"p{i}_sim.fa"),
         os.path.join(outdir, f"pair_{i}_A.fa"),
         os.path.join(outdir, f"pair_{i}_B.fa"),
         bl)
        for i in range(n_pairs)
    ]

    results = [None] * n_pairs
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        futures = {ex.submit(split_one, *t): t[0] for t in tasks}
        for fut in as_completed(futures):
            idx = futures[fut]
            results[idx] = fut.result()

    meta_path = os.path.join(outdir, "pairs_meta.tsv")
    with open(meta_path, "w") as fh:
        fh.write("pair_id\tbranch_len\tlen_A\tlen_B\n")
        for r in results:
            if r is not None:
                fh.write(f"{r[0]}\t{r[1]:.6f}\t{r[2]}\t{r[3]}\n")


if __name__ == "__main__":
    main()
