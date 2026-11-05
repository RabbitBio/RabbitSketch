#!/usr/bin/env python3
"""Analyze registered synthetic batches without using validation to choose parameters."""
import argparse
import gzip
import json

from common import *

META = ['regime_id', 'k', 'n_a', 'n_b', 'size_ratio', 'nominal_containment', 'truth_j', 'truth_c']


def load_batch(directory, units):
    complete = json.loads((directory / 'completion.json').read_text())
    path = directory / 'raw.jsonl.gz'
    if sha(path) != complete['raw_sha256']: raise ValueError('raw hash mismatch')
    with gzip.open(path, 'rt') as stream:
        rows = [json.loads(line) for line in stream]
    d = pd.DataFrame(rows)
    if len(d) != len(units) * 29: raise ValueError('wrong registered row count')
    if set(d.unit_id) != set(units.unit_id): raise ValueError('missing or unknown input')
    if d.duplicated(['unit_id', 'algorithm', 'parameter', 'requested_budget', 'half_subk']).any():
        raise ValueError('duplicate native result')
    if not (d.groupby('unit_id').size() == 29).all(): raise ValueError('incomplete input cohort')
    if not d.algorithm.isin(ALGORITHMS).all(): raise ValueError('unknown backend')
    fixed = d[d.algorithm != 'kssd']
    if not fixed.requested_budget.isin(CONFIG['pair_budgets']).all(): raise ValueError('unregistered budget')
    expected = np.where(fixed.algorithm.isin(['minhash', 'probminhash']), fixed.requested_budget // 16,
        np.where(fixed.algorithm == 'bindash', fixed.requested_budget // 4,
        np.where(fixed.algorithm.isin(['hll', 'setsketch']), np.log2(fixed.requested_budget / 2),
        np.ceil(8 * (fixed.n_a + fixed.n_b) / (CONFIG['frac_expected_budget_fraction'] * fixed.requested_budget)))))
    if not np.array_equal(fixed.parameter.to_numpy(), expected): raise ValueError('unregistered native parameter')
    kssd = d[d.algorithm == 'kssd']
    if set(zip(kssd.parameter, kssd.half_subk)) != set(map(tuple, CONFIG['kssd_resolutions'])):
        raise ValueError('unregistered KSSD resolution')
    if not (kssd.requested_budget == 0).all(): raise ValueError('KSSD was given a per-input oracle budget')
    identity = units.set_index('unit_id')
    for name in ['k', 'n_a', 'n_b', 'intersection', 'replicate', 'data_seed']:
        if not (d[name].to_numpy() == d.unit_id.map(identity[name]).to_numpy()).all():
            raise ValueError(f'registered {name} differs')
    expected_seed = d.unit_id.map(identity.hash_seed).where(d.algorithm != 'kssd', 42)
    if not (d.hash_seed == expected_seed).all(): raise ValueError('unregistered hash seed')
    if (d.groupby('unit_id')[['input_digest_a', 'input_digest_b']].nunique() != 1).any().any():
        raise ValueError('algorithms consumed different input')
    if not d.status.isin(['ok', 'invalid_estimate']).all(): raise ValueError('native exception')
    for name in ['payload_bytes_a', 'payload_bytes_b']:
        if (d[name] < 0).any(): raise ValueError('negative payload')
    metadata = units[['unit_id', 'regime_id', 'split', 'size_ratio', 'nominal_containment', 'truth_j', 'truth_c']]
    d = d.merge(metadata, on='unit_id', validate='many_to_one')
    d['zone'] = d.truth_c.map(zone)
    return add_errors(d)


def choose_kssd(d):
    """Resolution fixed by calibration bytes, never by validation error or payload."""
    keys = ['regime_id', 'parameter', 'half_subk']
    kssd = d[d.algorithm == 'kssd']
    calibration = kssd[kssd.split == 'calibration']
    sizes = calibration.groupby(keys).pair_bytes.agg(['max', 'mean']).reset_index()
    chosen, records = [], []
    for budget in CONFIG['pair_budgets']:
        fits = sizes[sizes['max'] <= budget].sort_values(['mean', 'parameter'], ascending=[False, True])
        best = fits.drop_duplicates('regime_id')
        chosen.append(kssd.merge(best[keys], on=keys, validate='many_to_one').assign(budget=budget))
        for row in best.to_dict('records'):
            records.append(dict(row, budget=budget))
    fixed = d[d.algorithm != 'kssd'].copy(); fixed['budget'] = fixed.requested_budget
    return pd.concat([fixed] + chosen, ignore_index=True), pd.DataFrame(records)


def decisions(selected, summary):
    rows, intervals = [], []
    for regime, d in selected.groupby('regime_id', sort=True):
        metadata = {k: d.iloc[0][k] for k in META}
        s = summary[summary.regime_id == regime]
        for metric in METRICS:
            base = s[s.metric == metric]
            cal = base[base.split == 'calibration'].set_index(['budget', 'algorithm'])
            val = base[base.split == 'validation'].set_index(['budget', 'algorithm'])
            allowed = {}; nominations = {}
            for budget in CONFIG['pair_budgets']:
                cs = cal.loc[budget]
                cs = cs[(cs.availability >= 1 - CONFIG['maximum_unavailable_fraction']) &
                        (cs.budget_violation_fraction <= CONFIG['maximum_budget_violation_fraction'])]
                nominations[budget] = cs.mae.idxmin() if len(cs) else ''
                values = {}
                for algorithm in cs.index:
                    if (budget, algorithm) not in val.index: continue
                    r = val.loc[budget, algorithm]
                    if r.availability < 1 - CONFIG['maximum_unavailable_fraction'] or r.budget_violation_fraction > CONFIG['maximum_budget_violation_fraction']: continue
                    observations = d[(d.split == 'validation') & (d.budget == budget) & (d.algorithm == algorithm)].sort_values('replicate')
                    if observations.replicate.tolist() != list(range(CONFIG['replicates'] // 2, CONFIG['replicates'])):
                        raise ValueError('unaligned validation inputs')
                    values[algorithm] = observations[metric].to_numpy()
                allowed[budget] = values
            seed = (int(regime[:12], 16) + CONFIG['bootstrap_seed'] + METRICS.index(metric)) % 2**32
            contrasts, critical = paired_family(allowed, seed)
            intervals.extend(dict(metadata, metric=metric, **row) for row in contrasts)
            for budget, values in allowed.items():
                relevant = [row for row in contrasts if row['budget'] == budget]
                surviving, winner = candidates(values, relevant)
                nominee = nominations[budget]
                point_best = min(values, key=lambda a: np.mean(values[a])) if values else ''
                best_error = np.mean(values[point_best]) if point_best else np.nan
                nominee_error = np.mean(values[nominee]) if nominee in values else np.nan
                rows.append(dict(metadata, metric=metric, budget=budget,
                    calibration_nominee=nominee, nominee_validates=int(nominee in values),
                    validation_point_best=point_best, validation_best_mae=best_error,
                    nominee_validation_mae=nominee_error,
                    nominee_regret=nominee_error - best_error,
                    eligible_count=len(values), eligible_candidates=';'.join(sorted(values)),
                    candidate_set=';'.join(surviving), supported_winner=winner,
                    status='supported_winner' if winner else ('no_candidate' if not values else
                        'single_eligible' if len(values) == 1 else 'unresolved'),
                    family_critical=critical))
    return pd.DataFrame(rows), pd.DataFrame(intervals)


def batch_analysis(directory, units):
    target = OUT / 'synthetic_batches' / directory.name
    target.mkdir(parents=True, exist_ok=True)
    identity = dict(raw_sha256=sha(directory / 'raw.jsonl.gz'),
        analyzer_sha256=sha(__file__), common_sha256=sha(HERE / 'common.py'),
        config_sha256=sha(HERE / 'config.json'))
    done = target / 'analysis_complete.json'
    if done.exists() and json.loads(done.read_text()) == identity:
        return
    d = load_batch(directory, units)
    selected, kssd = choose_kssd(d)
    s = summaries(selected, META + ['split', 'algorithm', 'parameter', 'half_subk', 'budget'])
    ds, comparisons = decisions(selected, s)
    s.to_csv(target / 'accuracy.csv', index=False)
    kssd.to_csv(target / 'kssd_calibration_choice.csv', index=False)
    ds.to_csv(target / 'decisions.csv', index=False)
    comparisons.to_csv(target / 'paired_intervals.csv', index=False)
    # Compact per-input validation evidence supports independent reanalysis.
    keep = ['unit_id', 'regime_id', 'replicate', 'hash_seed', 'split', 'algorithm', 'parameter', 'half_subk',
            'budget', 'pair_bytes', 'empty', 'jaccard', 'truth_j', 'truth_c', 'j_error', 'j_bias',
            'c_error', 'native_c_error', 'card_relative_error']
    save_gzip(target / 'selected_rows.csv.gz', selected[keep])
    # Marginal sensitivity tables describe the balanced grid, not new independent data.
    sensitivity = summaries(selected, ['k', 'n_b', 'budget', 'algorithm', 'split', 'hash_seed'])
    sensitivity.to_csv(target / 'seed_sensitivity.csv', index=False)
    raw_configuration = summaries(d.assign(budget=d.requested_budget.replace(0, np.nan)),
        META + ['split', 'algorithm', 'parameter', 'half_subk', 'requested_budget'])
    raw_configuration.to_csv(target / 'all_configurations.csv', index=False)
    # Support rows expose why small query / large reference conditions can fail.
    f = selected[(selected.algorithm == 'fracminhash') & (selected.split == 'validation')].copy()
    f['retained_query'] = f.payload_bytes_a / 8
    support = f.groupby(META + ['budget'], observed=True).agg(mean_query_hashes=('retained_query', 'mean'),
        min_query_hashes=('retained_query', 'min'), empty_fraction=('empty', 'mean'),
        mae_c=('c_error', 'mean'), mae_j=('j_error', 'mean')).reset_index()
    support.to_csv(target / 'frac_support.csv', index=False)
    # Per-object cardinality evidence keeps actual object size visible.
    card = selected[~selected.algorithm.isin(['kssd', 'probminhash'])].copy()
    object_frames = []
    for side in ['a', 'b']:
        estimate = card['cardinality_' + side]
        valid = np.isfinite(estimate) & (estimate > 0) & (card['empty_' + side] == 0)
        objects = card[['regime_id', 'unit_id', 'k', 'split', 'algorithm', 'parameter', 'budget']].copy()
        objects['side'] = side; objects['n'] = card['n_' + side]
        objects['object_payload'] = card['payload_bytes_' + side]
        objects['relative_error'] = np.where(valid, abs(estimate / objects.n - 1), np.nan)
        object_frames.append(objects)
    objects = pd.concat(object_frames, ignore_index=True)
    object_summary = objects.groupby(['k', 'split', 'algorithm', 'budget', 'n']).agg(
        observations=('relative_error', 'size'), valid=('relative_error', 'count'),
        mean_object_bytes=('object_payload', 'mean'), max_object_bytes=('object_payload', 'max'),
        relative_mae=('relative_error', 'mean'), relative_p95=('relative_error', lambda x: x.quantile(.95))).reset_index()
    object_summary.to_csv(target / 'object_cardinality.csv', index=False)
    atomic(done, identity)
    print(json.dumps(dict(batch=directory.name, rows=len(d), decisions=len(ds),
                          states=ds.status.value_counts().to_dict())), flush=True)


def consolidate():
    directories = sorted((RUN / 'synthetic').glob('*'))
    expected = len(CONFIG['k']) * len(CONFIG['larger_cardinalities'])
    if len(directories) != expected or any(not (p / 'completion.json').exists() for p in directories):
        raise ValueError('synthetic native grid incomplete')
    tables = ['accuracy', 'kssd_calibration_choice', 'decisions', 'paired_intervals',
              'seed_sensitivity', 'all_configurations', 'frac_support', 'object_cardinality']
    for name in tables:
        frames = [pd.read_csv(OUT / 'synthetic_batches' / p.name / (name + '.csv')) for p in directories]
        save('synthetic_' + name, pd.concat(frames, ignore_index=True))
    s = pd.read_csv(OUT / 'synthetic_accuracy.csv')
    s['zone'] = s.truth_c.map(zone)
    # A nominal-grid descriptive average; record truth duplicates explicitly.
    good = s[(s.split == 'validation') & (s.availability >= .99) & (s.budget_violation_fraction <= .01)]
    keys = ['k', 'budget', 'size_ratio', 'zone', 'algorithm', 'metric']
    pooled = good.groupby(keys).agg(cells=('regime_id', 'nunique'), macro_mae=('mae', 'mean'),
        median_cell_mae=('mae', 'median'), worst_cell_p95=('p95', 'max'),
        minimum_availability=('availability', 'min')).reset_index()
    save('synthetic_pooled_accuracy', pooled)
    repeats = []
    keys = ['metric', 'budget', 'size_ratio', 'nominal_containment']
    decision_table = pd.read_csv(OUT / 'synthetic_decisions.csv')
    for group, cells in decision_table.groupby(keys):
        for algorithm in ALGORITHMS:
            wins = cells[cells.supported_winner == algorithm]
            repeats.append(dict(zip(keys, group), algorithm=algorithm, tested_cells=len(cells),
                supported_cells=len(wins), supported_k=','.join(map(str, sorted(wins.k.unique()))),
                supported_cardinalities=','.join(map(str, sorted(wins.n_b.unique()))),
                supported_everywhere=bool(len(wins) == len(cells))))
    save('synthetic_repeatability', pd.DataFrame(repeats))
    decisions = pd.read_csv(OUT / 'synthetic_decisions.csv')
    save('synthetic_decision_counts', decisions.groupby(['metric', 'budget', 'status', 'supported_winner'], dropna=False).size().rename('cells').reset_index())
    atomic(OUT / 'synthetic_complete.json', dict(status='complete', native_batches=expected,
        native_measurements=107520 * 29, nominal_regimes=840,
        realized_conditions=len(s.drop_duplicates(['k', 'n_a', 'n_b', 'truth_j'])),
        intervals_scope='simultaneous candidate contrasts across budgets within each regime and metric; not across regimes',
        analysis_sha256=sha(__file__), common_sha256=sha(HERE / 'common.py')))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--completed-only', action='store_true'); args = ap.parse_args()
    units = pd.read_csv(RUN / 'synthetic_units.tsv', sep='\t', dtype={'data_seed': 'uint64'})
    for (k, nb), group in units.groupby(['k', 'n_b'], sort=True):
        directory = RUN / 'synthetic' / f'k{k}_n{nb}'
        if not (directory / 'completion.json').exists():
            if args.completed_only: continue
            raise ValueError(f'pending native batch: {directory}')
        batch_analysis(directory, group)
    if not args.completed_only: consolidate()


if __name__ == '__main__': main()
