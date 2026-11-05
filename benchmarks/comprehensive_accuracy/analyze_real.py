#!/usr/bin/env python3
"""Complete real/task panels, fixed calibration parameters, genus-aware evidence."""
import argparse
import gzip
import json

from common import *


def raw_job(job):
    directory = RUN / 'real' / job['job_id']
    completion = json.loads((directory / 'completion.json').read_text())
    path = directory / 'raw.jsonl.gz'
    if sha(path) != completion['raw_sha256']: raise ValueError('real raw digest changed')
    rows = {'pair': [], 'object': []}
    counts = {}
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            row = json.loads(line); typ = row['type']
            counts[typ] = counts.get(typ, 0) + 1
            if typ in rows:
                row = {k: v for k, v in row.items() if 'seconds' not in k and 'rss' not in k and k != 'input_bytes'}
                row.update(phase=job['phase'], half_subk=job['half_subk'], job_id=job['job_id'])
                if job['phase'] != 'weighted' and row['algorithm'] == 'probminhash_unweighted':
                    row['algorithm'] = 'probminhash'
                rows[typ].append(row)
    if counts != completion['rows']: raise ValueError('real archive row counts differ')
    return {typ: pd.DataFrame(values) for typ, values in rows.items()}


def metadata(phase):
    d = pd.read_csv(RUN / 'data' / (phase + '_pairs.tsv'), sep='\t')
    if phase == 'mixed':
        refs = sorted(d.reference.unique(), key=lambda x: sha_text('mixed-calibration:' + x))
        split = {r: 'calibration' if i < len(refs) // 2 else 'validation' for i, r in enumerate(refs)}
        d['split'] = d.reference.map(split); d['cluster'] = d.reference
    if phase == 'containment': d['cluster'] = d.reference
    return d


def sha_text(text):
    import hashlib
    return hashlib.sha256(text.encode()).hexdigest()


def compare_frozen_truth(truth, phase, k):
    if phase in ['natural', 'mixed']:
        path = ROOT / 'analysis/robustness_20261001/exact_truth_and_models.csv'
        prior = pd.read_csv(path)
        prior = prior[(prior.panel == phase) & (prior.k == k)]
    else:
        path = ROOT / 'analysis/regimes_20261001/containment_exact_truth.csv'
        prior = pd.read_csv(path); prior = prior[prior.k == k]
    if prior.empty: return dict(previous_pairs=0, matched_pairs=0, note='new k, no earlier exact panel')
    columns = ['pair_id', 'n_a', 'n_b', 'intersection']
    joined = prior[columns].merge(truth[columns], on='pair_id', suffixes=('_old', '_new'), validate='one_to_one')
    if len(joined) != len(prior): raise ValueError('earlier pairs missing from expanded panel')
    for reducer in [np.minimum, np.maximum]:
        old = reducer(joined.n_a_old, joined.n_b_old)
        new = reducer(joined.n_a_new, joined.n_b_new)
        if not np.array_equal(old, new): raise ValueError('exact cardinality changed across frozen runs')
    if not np.array_equal(joined.intersection_old, joined.intersection_new):
        raise ValueError('exact intersection changed across frozen runs')
    return dict(previous_pairs=len(prior), matched_pairs=len(joined), exact_counts_match=True, prior_sha256=sha(path))


def load_panel(jobs, phase, k):
    chosen = [j for j in jobs if j['phase'] == phase and j['k'] == k]
    frames = [raw_job(j) for j in chosen]
    pairs = pd.concat([f['pair'] for f in frames], ignore_index=True)
    objects = pd.concat([f['object'] for f in frames], ignore_index=True)
    keys = ['algorithm', 'parameter', 'half_subk', 'seed']
    if pairs.duplicated(keys + ['pair_id']).any() or objects.duplicated(keys + ['object_id']).any():
        raise ValueError('duplicate real observation')
    exact_objects = objects[objects.algorithm == 'exact'].set_index('object_id')
    nonexact = objects[objects.algorithm != 'exact']
    # int conversion occurs only on integer-typed fingerprints, never iterrows float coercion.
    observed = nonexact.input_digest.to_numpy(dtype=np.uint64)
    expected = nonexact.object_id.map(exact_objects.input_digest).to_numpy(dtype=np.uint64)
    if not np.array_equal(observed, expected): raise ValueError('normalized real inputs differ')
    truth = pairs[pairs.algorithm == 'exact'][['pair_id', 'n_a', 'n_b', 'intersection']].copy()
    truth['truth_j'] = truth.intersection / (truth.n_a + truth.n_b - truth.intersection)
    truth['truth_c'] = truth.intersection / truth.n_a
    meta = metadata(phase)
    truth = truth.merge(meta, on='pair_id', validate='one_to_one')
    native = pairs[pairs.algorithm != 'exact'].drop(columns=['n_a', 'n_b', 'intersection'], errors='ignore')
    native = native.merge(truth.drop(columns=['a', 'b']), on='pair_id', validate='many_to_one')
    ca = nonexact[keys + ['object_id', 'cardinality']]
    native = native.merge(ca.rename(columns={'object_id': 'a', 'cardinality': 'cardinality_a'}),
                          on=keys + ['a'], validate='many_to_one')
    native = native.merge(ca.rename(columns={'object_id': 'b', 'cardinality': 'cardinality_b'}),
                          on=keys + ['b'], validate='many_to_one')
    native['containment_ab'] = bounded_containment(native.jaccard, native.cardinality_a, native.cardinality_b)
    native['native_containment_ab'] = native.containment_native_ab
    d = add_errors(native)
    # Preserve alignment-model residuals separately; they are not Jaccard truth.
    d['ani_from_j'] = np.power(2 * d.jaccard / (1 + d.jaccard), 1 / k)
    d['exact_ani_from_j'] = np.power(2 * d.truth_j / (1 + d.truth_j), 1 / k)
    d['ani_sampling_error'] = np.where(d.j_error.notna(), abs(d.ani_from_j - d.exact_ani_from_j), np.nan)
    if phase == 'mixed':
        d['column_ani_error'] = np.where(d.j_error.notna(), abs(d.ani_from_j - d.column_identity), np.nan)
        d['exact_model_error'] = abs(d.exact_ani_from_j - d.column_identity)
    return d, truth, objects


def fixed_parameters(d):
    """One native parameter per panel/k/algorithm/budget, frozen on calibration."""
    keys = ['algorithm', 'parameter', 'half_subk']
    cal = d[d.split == 'calibration']
    sizes = cal.groupby(keys).pair_bytes.agg(['max', 'mean']).reset_index()
    selected, choices = [], []
    for budget in CONFIG['pair_budgets']:
        fits = sizes[sizes['max'] <= budget].sort_values(['mean', 'parameter', 'half_subk'], ascending=[False, True, True])
        fixed = fits.drop_duplicates('algorithm')
        selected.append(d.merge(fixed[keys], on=keys, validate='many_to_one').assign(budget=budget))
        choices.append(fixed.assign(budget=budget))
    return pd.concat(selected, ignore_index=True), pd.concat(choices, ignore_index=True)


def stratify(d, phase):
    rows = [d.assign(stratum='all')]
    if phase == 'natural':
        within = d[d.genus_a == d.genus_b]
        rows.append(within.assign(stratum='within_genus'))
        rows.append(within[within.species_a == within.species_b].assign(stratum='within_species'))
        rows.append(within[within.species_a != within.species_b].assign(stratum='different_species_same_genus'))
        rows.append(d[d.genus_a != d.genus_b].assign(stratum='cross_genus'))
        for label, lo, hi in [('J_0_001', 0., .01), ('J_001_01', .01, .1), ('J_01_05', .1, .5), ('J_05_1', .5, 1.000001)]:
            rows.append(d[(d.truth_j >= lo) & (d.truth_j < hi)].assign(stratum=label))
    elif phase == 'mixed':
        for (target, indel), group in d.groupby(['target_identity', 'indel_share']):
            rows.append(group.assign(stratum=f'identity{target:g}_indel{indel:g}'))
    elif phase == 'containment':
        for (ratio, c), group in d.groupby(['nominal_ratio', 'nominal_containment']):
            rows.append(group.assign(stratum=f'ratio{ratio:g}_containment{c:g}'))
    return pd.concat(rows, ignore_index=True)


def cluster_intervals(d, phase):
    """Conditional on tested seeds; contrasts use common biological pairs."""
    import itertools
    rows = []
    metrics = METRICS + ['ani_sampling_error'] + (['column_ani_error'] if phase == 'mixed' else [])
    validation = stratify(d[d.split == 'validation'], phase)
    for (budget, stratum), group in validation.groupby(['budget', 'stratum']):
        eligible = []
        for alg, g in group.groupby('algorithm'):
            if (g.pair_bytes > g.budget).mean() <= .01 and g.j_error.notna().mean() >= .99:
                eligible.append(alg)
        for metric in metrics:
            # Conditional MAE alone must not hide missing-output problems.
            usable = [a for a in eligible if group[group.algorithm == a][metric].notna().mean() >= .99]
            for left, right in itertools.combinations(sorted(usable), 2):
                a = group[group.algorithm == left].groupby(['pair_id', 'cluster'])[metric].mean().rename('a')
                b = group[group.algorithm == right].groupby(['pair_id', 'cluster'])[metric].mean().rename('b')
                paired = pd.concat([a, b], axis=1).dropna().reset_index()
                if paired.empty: continue
                # Cross-genus pairs are not an independent cluster called "cross_genus".
                # Keep them in descriptive tables; formal natural intervals use within-genus pairs.
                if phase == 'natural': paired = paired[paired.cluster != 'cross_genus']
                if paired.empty: continue
                paired['difference'] = paired.a - paired.b
                cluster = paired.groupby('cluster').difference.mean()
                lo, hi = bootstrap_clusters(cluster.to_numpy())
                rows.append(dict(phase=phase, budget=budget, stratum=stratum, metric=metric,
                    left=left, right=right, common_pairs=len(paired), clusters=len(cluster),
                    cluster_macro_difference=cluster.mean(), ci_low=lo, ci_high=hi,
                    scope='pointwise cluster macro; within-genus only for natural; no multiplicity adjustment'))
    return pd.DataFrame(rows)


def analyze_panel(jobs, phase, k):
    target = OUT / 'real_panels' / f'{phase}_k{k}'; target.mkdir(parents=True, exist_ok=True)
    identity = dict(analyzer_sha256=sha(__file__), common_sha256=sha(HERE / 'common.py'),
        completions={j['job_id']: sha(RUN / 'real' / j['job_id'] / 'completion.json')
                     for j in jobs if j['phase'] == phase and j['k'] == k})
    if (target / 'analysis_complete.json').exists() and json.loads((target / 'analysis_complete.json').read_text()) == identity:
        return
    d, truth, objects = load_panel(jobs, phase, k)
    atomic(target / 'frozen_truth_comparison.json', compare_frozen_truth(truth, phase, k))
    selected, choices = fixed_parameters(d)
    metrics = METRICS + ['native_c_error', 'ani_sampling_error'] + (['column_ani_error'] if phase == 'mixed' else [])
    summaries(stratify(selected, phase), ['k', 'budget', 'algorithm', 'parameter', 'half_subk', 'split', 'stratum'], metrics).to_csv(target / 'budget_accuracy.csv', index=False)
    summaries(d, ['k', 'algorithm', 'parameter', 'half_subk', 'split'], metrics).to_csv(target / 'configuration_accuracy.csv', index=False)
    choices.to_csv(target / 'calibration_parameters.csv', index=False)
    truth.assign(k=k).to_csv(target / 'exact_truth.csv', index=False)
    save_gzip(target / 'selected_rows.csv.gz', selected)
    save_gzip(target / 'objects.csv.gz', objects)
    save_gzip(target / 'all_pair_rows.csv.gz', d)
    cluster_intervals(selected, phase).to_csv(target / 'cluster_intervals.csv', index=False)
    coverage = []
    for split, base in truth.groupby('split'):
        for budget in CONFIG['pair_budgets']:
            for alg in ALGORITHMS:
                g = selected[(selected.split == split) & (selected.budget == budget) & (selected.algorithm == alg)]
                coverage.append(dict(phase=phase, k=k, split=split, budget=budget, algorithm=alg,
                    total_pairs=len(base), measured_pairs=g.pair_id.nunique(),
                    parameter_chosen=int(len(g) > 0), valid_observations=g.j_error.notna().sum(),
                    observations=len(g), budget_violations=(g.pair_bytes > budget).sum()))
    pd.DataFrame(coverage).to_csv(target / 'coverage.csv', index=False)
    atomic(target / 'analysis_complete.json', identity)
    print(json.dumps(dict(phase=phase, k=k, measurements=len(d), biological_pairs=len(truth))), flush=True)


def analyze_weighted(jobs):
    frames = [raw_job(j)['pair'] for j in jobs if j['phase'] == 'weighted']
    d = pd.concat(frames, ignore_index=True)
    meta = metadata('weighted')
    d = d.merge(meta.drop(columns=['a', 'b']), on='pair_id', validate='many_to_one')
    d['error'] = abs(d.jaccard - d.probability_jaccard)
    d['bias'] = d.jaccard - d.probability_jaccard
    d['squared'] = d.error ** 2
    keys = ['family', 'model', 'perturbation', 'split', 'parameter']
    summary = d.groupby(keys, dropna=False).agg(pairs=('pair_id', 'nunique'), observations=('pair_id', 'size'),
        seeds=('seed', 'nunique'), mae=('error', 'mean'), p95=('error', lambda x: x.quantile(.95)),
        mse=('squared', 'mean'), bias=('bias', 'mean')).reset_index()
    summary['rmse'] = np.sqrt(summary.pop('mse')); summary['pair_payload'] = summary.parameter * 16
    save('weighted_accuracy', summary)
    save_gzip(OUT / 'weighted_rows.csv.gz', d)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--completed-only', action='store_true'); args = ap.parse_args()
    jobs = json.loads((RUN / 'plan.json').read_text())['real_jobs']
    for phase in ['natural', 'mixed', 'containment']:
        for k in (CONFIG['real_k'] if phase == 'natural' else [24]):
            panel = [j for j in jobs if j['phase'] == phase and j['k'] == k]
            if any(not (RUN / 'real' / j['job_id'] / 'completion.json').exists() for j in panel):
                if args.completed_only: continue
                raise ValueError(f'incomplete real panel {phase} k={k}')
            analyze_panel(jobs, phase, k)
    weighted = [j for j in jobs if j['phase'] == 'weighted']
    if all((RUN / 'real' / j['job_id'] / 'completion.json').exists() for j in weighted): analyze_weighted(jobs)
    elif not args.completed_only: raise ValueError('incomplete weighted panel')
    if not args.completed_only:
        for table in ['budget_accuracy', 'configuration_accuracy', 'calibration_parameters', 'cluster_intervals', 'coverage']:
            frames = [pd.read_csv(p).assign(phase=p.parent.name.split('_k')[0], k=int(p.parent.name.split('_k')[1]))
                      for p in sorted((OUT / 'real_panels').glob('*/' + table + '.csv'))]
            save('real_' + table, pd.concat(frames, ignore_index=True))
        atomic(OUT / 'real_complete.json', dict(status='complete', native_jobs=len(jobs),
            analyzer_sha256=sha(__file__), common_sha256=sha(HERE / 'common.py')))


if __name__ == '__main__': main()
