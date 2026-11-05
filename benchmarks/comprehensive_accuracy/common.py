"""Shared accuracy-only analysis, with explicit missing-output handling."""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = Path('/home/bigssd/zt/rabbitsketch2_benchmark')
RUN = ROOT / 'runs/comprehensive_accuracy_20261030'
OUT = ROOT / 'analysis/comprehensive_accuracy_20261030'
PUBLISHED = HERE.parent / 'results/comprehensive_accuracy_20261030'
CONFIG = json.loads((HERE / 'config.json').read_text())
ALGORITHMS = ['minhash', 'fracminhash', 'hll', 'setsketch', 'kssd', 'probminhash', 'bindash']
NAMES = dict(minhash='FastKMV', fracminhash='FracMinHash', hll='HLL',
             setsketch='SetSketch-fixed', kssd='KSSD', probminhash='ProbMinHash-unweighted',
             bindash='BinDash', exact='Exact sets')
COLORS = dict(minhash='#0072B2', fracminhash='#CC79A7', hll='#56B4E9',
              setsketch='#009E73', kssd='#E69F00', probminhash='#332288', bindash='#D55E00')
METRICS = ['j_error', 'c_error', 'card_relative_error']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def atomic(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    tmp.replace(path)


def save(name, frame):
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / (name + '.csv'), index=False)


def save_gzip(path, frame):
    frame.to_csv(path, index=False, compression={'method': 'gzip', 'mtime': 0})


def zone(value):
    if value == 0: return 'zero'
    if value <= .1: return 'low'
    if value < .9: return 'middle'
    if value < 1: return 'high'
    return 'full'


def bounded_containment(j, ca, cb):
    j, ca, cb = map(np.asarray, [j, ca, cb])
    valid = np.isfinite(j) & np.isfinite(ca) & np.isfinite(cb) & (ca > 0) & (cb > 0)
    result = np.full(j.shape, np.nan)
    result[valid] = np.clip(j[valid] * (ca[valid] + cb[valid]) / (1 + j[valid]),
                            0, np.minimum(ca[valid], cb[valid])) / ca[valid]
    return result


def add_errors(d):
    d = d.copy()
    d['pair_bytes'] = d.payload_bytes_a + d.payload_bytes_b
    d['empty'] = (d.empty_a != 0) | (d.empty_b != 0)
    good = np.isfinite(d.jaccard) & d.jaccard.between(0, 1) & ~d['empty']
    d['j_error'] = np.where(good, abs(d.jaccard - d.truth_j), np.nan)
    d['j_bias'] = np.where(good, d.jaccard - d.truth_j, np.nan)
    cgood = np.isfinite(d.containment_ab) & d.containment_ab.between(0, 1) & ~d['empty']
    d['c_error'] = np.where(cgood & ~d.algorithm.isin(['kssd', 'probminhash']),
                            abs(d.containment_ab - d.truth_c), np.nan)
    ca, cb = d.cardinality_a, d.cardinality_b
    cardgood = np.isfinite(ca) & np.isfinite(cb) & (ca > 0) & (cb > 0) & ~d['empty']
    cardgood &= ~d.algorithm.isin(['kssd', 'probminhash'])
    d['card_relative_error'] = np.where(cardgood, (abs(ca / d.n_a - 1) + abs(cb / d.n_b - 1)) / 2, np.nan)
    d['native_c_error'] = np.where(np.isfinite(d.native_containment_ab) & ~d['empty'] & ~d.algorithm.isin(['kssd', 'probminhash']),
                                   abs(d.native_containment_ab - d.truth_c), np.nan)
    return d


def summaries(frame, keys, metrics=METRICS):
    frame = frame.copy()
    frame['_zero'] = frame.jaccard == 0
    frame['_budget_bad'] = frame.pair_bytes > frame.budget if 'budget' in frame else np.nan
    grouped = frame.groupby(keys, observed=True, dropna=False, sort=True)
    base = grouped.agg(observations=('jaccard', 'size'), mean_pair_bytes=('pair_bytes', 'mean'),
        max_pair_bytes=('pair_bytes', 'max'), empty_fraction=('empty', 'mean'),
        zero_fraction=('_zero', 'mean'), budget_violation_fraction=('_budget_bad', 'mean'))
    truth = grouped.truth_j.agg(['first', 'nunique'])
    rows = []
    for metric in metrics:
        values = grouped[metric]
        result = base.copy()
        result['metric'] = metric
        result['valid'] = values.count()
        result['availability'] = result.valid / result.observations
        result['mae'] = values.mean()
        result['p95'] = values.quantile(.95)
        squared = frame[metric] ** 2
        result['rmse'] = np.sqrt(frame.assign(_squared=squared).groupby(keys, observed=True, dropna=False)._squared.mean())
        result['bias'] = grouped.j_bias.mean() if metric == 'j_error' else np.nan
        result['relative_mae'] = (result.mae / truth['first'].where((truth['nunique'] == 1) & (truth['first'] > 0))) if metric == 'j_error' else np.nan
        rows.append(result.reset_index())
    return pd.concat(rows, ignore_index=True)


def paired_family(values_by_budget, seed, samples=None):
    """Paired input bootstrap, sharing draws across every budget in one regime.

    Input maps budget -> {algorithm: 1D losses on the same independent units}.
    Eligibility must have been checked first. No mean imputation for missing data.
    Output intervals are simultaneous across all contrasts in this family.
    """
    samples = samples or CONFIG['bootstrap_samples']
    lengths = {len(v) for algorithms in values_by_budget.values() for v in algorithms.values()}
    if not lengths: return [], np.nan
    if len(lengths) != 1 or next(iter(lengths)) < 2: raise ValueError('unpaired input cohorts')
    n = next(iter(lengths))
    rng = np.random.default_rng(seed)
    weights = rng.multinomial(n, np.repeat(1 / n, n), size=samples) / n
    records, differences = [], []
    for budget, algorithms in sorted(values_by_budget.items()):
        for left, right in itertools.combinations(sorted(algorithms), 2):
            a, b = np.asarray(algorithms[left]), np.asarray(algorithms[right])
            if not np.isfinite(a).all() or not np.isfinite(b).all():
                raise ValueError('missing output in paired family')
            differences.append(a - b)
            records.append(dict(budget=budget, left=left, right=right, mae_left=a.mean(), mae_right=b.mean()))
    if not differences: return [], np.nan
    differences = np.array(differences)
    means = differences.mean(axis=1)
    ses = differences.std(axis=1, ddof=1) / np.sqrt(n)
    boots = weights @ differences.T
    active = ses > 0
    critical = float(np.quantile(np.max(abs((boots[:, active] - means[active]) / ses[active]), axis=1),
                                 1 - CONFIG['alpha'])) if active.any() else 0.
    low, high = np.quantile(boots, [CONFIG['alpha'] / 2, 1 - CONFIG['alpha'] / 2], axis=0)
    for i, record in enumerate(records):
        record.update(difference=means[i], se=ses[i], pointwise_low=low[i], pointwise_high=high[i],
                      simultaneous_low=means[i] - critical * ses[i],
                      simultaneous_high=means[i] + critical * ses[i], family_critical=critical,
                      independent_units=n, family_contrasts=len(records))
    return records, critical


def preferred_left(row):
    """Never turn a degenerate sample or a tiny significant difference into a win."""
    required = max(CONFIG['minimum_absolute_advantage'],
                   CONFIG['minimum_relative_advantage'] * row['mae_right'])
    return row['se'] > 0 and row['simultaneous_high'] < 0 and row['mae_right'] - row['mae_left'] >= required


def candidates(algorithms, comparisons):
    beaten = set()
    defeats = {a: set() for a in algorithms}
    for row in comparisons:
        if preferred_left(row):
            defeats[row['left']].add(row['right']); beaten.add(row['right'])
        reverse = dict(row, mae_left=row['mae_right'], mae_right=row['mae_left'],
                       simultaneous_high=-row['simultaneous_low'])
        if preferred_left(reverse):
            defeats[row['right']].add(row['left']); beaten.add(row['left'])
    remaining = sorted(set(algorithms) - beaten)
    winners = [a for a in algorithms if len(algorithms) > 1 and defeats[a] == set(algorithms) - {a}]
    return remaining, winners[0] if len(winners) == 1 else ''


def bootstrap_clusters(values, seed=20261030, samples=4000):
    values = np.asarray(values, dtype=float)
    if len(values) < 2 or not np.isfinite(values).all(): return np.nan, np.nan
    rng = np.random.default_rng(seed)
    means = values[rng.integers(0, len(values), (samples, len(values)))].mean(axis=1)
    return tuple(np.quantile(means, [.025, .975]))
