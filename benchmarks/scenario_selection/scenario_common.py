"""Paths and data invariants for the independent scenario-selection extension."""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import re

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = Path('/home/bigssd/zt/rabbitsketch2_benchmark')
TAG = 'scenario_selection_20261105'
DATA, RUN, OUT = (ROOT / parent / TAG for parent in ['data', 'runs', 'analysis'])
PUBLISHED = HERE.parent / 'results' / TAG
OLD = ROOT / 'data/robustness_20261001'
BASE = Path('/home/bigssd/zt/ref')
CONFIG = json.loads((HERE / 'config.json').read_text())
BINARY = ROOT / 'build_comprehensive_20261030/cmake/rabbitsketch_accuracy'
ALGORITHMS = ['minhash', 'fracminhash', 'hll', 'setsketch', 'kssd', 'probminhash_unweighted', 'bindash']
NAMES = dict(minhash='FastKMV', fracminhash='FracMinHash', hll='HLL', setsketch='SetSketch-fixed',
             kssd='KSSD', probminhash_unweighted='ProbMinHash', bindash='BinDash', exact='Exact sets')
RC = bytes.maketrans(b'ACGTN', b'TGCAN')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def key(value):
    return hashlib.sha256(('rs2-scenario-20261105:' + str(value)).encode()).hexdigest()


def atomic(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temp.replace(path)


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def tsv(path, rows):
    if not rows: raise ValueError('empty table: ' + str(path))
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def accession_root(text):
    m = re.search(r'GC[AF]_(\d+)', str(text))
    return m[1] if m else None


def fasta_records(path):
    raw = gzip.decompress(Path(path).read_bytes()) if str(path).endswith('.gz') else Path(path).read_bytes()
    result = []
    for block in raw.split(b'>')[1:]:
        header, sep, seq = block.partition(b'\n')
        if not sep: continue
        seq = re.sub(rb'[^ACGT]', b'N', re.sub(rb'\s+', b'', seq).upper())
        if seq: result.append(seq)
    if not result: raise ValueError('empty FASTA: ' + str(path))
    return result


def fingerprint(records):
    parts = sorted(hashlib.sha256(min(s, s.translate(RC)[::-1])).digest() for s in records)
    return hashlib.sha256(b''.join(parts)).hexdigest()


def inspect(path):
    records = fasta_records(path)
    lengths = sorted(map(len, records), reverse=True); total = sum(lengths)
    n50 = lengths[int(np.searchsorted(np.cumsum(lengths), total / 2))]
    return dict(fingerprint=fingerprint(records), source_sha256=sha(path), bases=total, n50=n50,
                acgt_fraction=sum(len(s) - s.count(b'N') for s in records) / total,
                records=len(records), max_record_length=max(lengths))


def bounded_containment(j, ca, cb):
    j, ca, cb = map(lambda x: np.asarray(x, dtype=float), [j, ca, cb])
    valid = np.isfinite(j) & np.isfinite(ca) & np.isfinite(cb) & (ca > 0) & (cb > 0)
    out = np.full(j.shape, np.nan)
    out[valid] = np.clip(j[valid] * (ca[valid] + cb[valid]) / (1 + j[valid]),
                         0, np.minimum(ca[valid], cb[valid])) / ca[valid]
    return out


def native_grid():
    b = CONFIG['pair_budgets']
    out = [(a, x // 16, 6) for a in ['minhash', 'probminhash_unweighted'] for x in b]
    out += [('bindash', x // 4, 6) for x in b]
    out += [(a, (x // 2).bit_length() - 1, 6) for a in ['hll', 'setsketch'] for x in b]
    out += [('fracminhash', p, 6) for p in CONFIG['frac_scaled']]
    out += [('kssd', p, h) for p, h in CONFIG['kssd_resolutions']]
    return out


def configuration(algorithm, parameter, half_subk=6):
    return f'{algorithm}:{int(parameter)}:{int(half_subk)}'


def save_csv(path, frame):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None)
