#!/usr/bin/env python3
"""A frozen, observable, precision-only sketch advisor; no performance claim."""
import argparse
import datetime

from scenario_common import *
from measure import META_FEATURES, PILOT_FEATURES, METRICS, STATISTICS

MODEL = OUT / 'model'


def nanmean(values, axis=0):
    valid = np.isfinite(values)
    count = valid.sum(axis=axis)
    return np.divide(np.where(valid, values, 0).sum(axis=axis), count,
                     out=np.full(count.shape, np.nan, dtype=float), where=count > 0)


def read_scores(split):
    frames = [pd.read_csv(p) for p in sorted((OUT / split / 'scores').glob('*.csv.gz'))]
    if not frames: raise ValueError('no scored configurations')
    return pd.concat(frames, ignore_index=True)


def payload_guard(features, k, algorithm, parameter, metric, policy):
    if algorithm in ['minhash', 'probminhash_unweighted']: return float(16 * parameter)
    if algorithm == 'bindash': return float(4 * parameter)
    if algorithm in ['hll', 'setsketch']: return float(2 * (2 ** parameter))
    if policy == 'pilot':
        nq = 10 ** float(features['log_pilot_query_n_max'])
        nr = 10 ** float(features['log_pilot_reference_n_max'])
        safety = CONFIG['pilot_cardinality_safety_factor']
    else:
        nq = 10 ** float(features['log_query_bases_max'])
        nr = 10 ** float(features['log_reference_bases_max'])
        safety = 1.  # Input bases upper-bound distinct k-mers.
    total = 2 * nq if metric in ['cardinality', 'clustering'] else nq + nr
    if algorithm == 'fracminhash': return float(8 * safety * total / parameter)
    if algorithm == 'kssd':
        item_bytes = 8 if k // 2 - parameter > 8 else 4
        return float(item_bytes * safety * total / 16 ** parameter)
    raise ValueError('unknown backend')


def capability(algorithm, metric):
    return not (metric in ['cardinality', 'containment'] and algorithm in ['kssd', 'probminhash_unweighted'])


def nearest_distinct(features, training, fields, scale, k, kind, groups, train_k, train_kind, count=8):
    x = np.array([float(features[f]) for f in fields])
    if not np.isfinite(x).all(): raise ValueError('non-finite observable feature')
    candidate = np.flatnonzero((train_k == k) & (train_kind == kind))
    if not len(candidate): raise ValueError('unsupported k/input kind')
    distance = np.sqrt(np.mean(((training[candidate] - x) / scale) ** 2, axis=1))
    order = np.argsort(distance, kind='stable')
    selected, seen, distances = [], set(), []
    for position in order:
        i = candidate[position]
        if groups[i] in seen: continue
        seen.add(groups[i]); selected.append(i); distances.append(float(distance[position]))
        if len(selected) == count: break
    minimum, maximum = training[candidate].min(axis=0), training[candidate].max(axis=0)
    outside = (x < minimum - .5 * scale) | (x > maximum + .5 * scale)
    supported = len(selected) == count and not outside.any() and float(np.mean(distances)) <= 3
    return np.array(selected), supported, float(np.mean(distances)), [f for f, b in zip(fields, outside) if b]


def fit():
    if not (OUT / 'development_scoring_complete.json').exists(): raise ValueError('finish development scoring first')
    if (RUN / 'validation-outcomes_plan.json').exists(): raise ValueError('cannot refit after validation outcomes start')
    MODEL.mkdir(parents=True, exist_ok=True)
    features = pd.read_csv(OUT / 'development_features.csv').sort_values(['k', 'workload_id']).reset_index(drop=True)
    d = read_scores('development'); d = d[d.algorithm != 'exact']
    cfgs = [configuration(*v) for v in native_grid()]
    cfg_index = {c: i for i, c in enumerate(cfgs)}
    work_index = {(r.workload_id, int(r.k)): i for i, r in enumerate(features.itertuples())}
    budget_index = {b: i for i, b in enumerate(CONFIG['pair_budgets'])}
    metric_index = {m: i for i, m in enumerate(METRICS)}
    scores = np.full((len(features), len(budget_index), len(cfgs), len(METRICS), len(STATISTICS)), np.nan)
    wi = np.array([work_index[v] for v in zip(d.workload_id, d.k)])
    bi = d.budget.map(budget_index).to_numpy(); ci = d.config.map(cfg_index).to_numpy(); mi = d.metric.map(metric_index).to_numpy()
    if d.duplicated(['workload_id', 'k', 'budget', 'config', 'metric']).any(): raise ValueError('duplicate calibration score')
    scores[wi, bi, ci, mi] = d[STATISTICS].to_numpy()
    if not np.isfinite(scores[..., 0]).all(): raise ValueError('missing calibration loss; failures must be explicit')
    scalings = {}
    arrays = dict(scores=scores.astype(np.float32), train_k=features.k.to_numpy(),
                  train_kind=features.kind.to_numpy(dtype=str), groups=features.genus.to_numpy(dtype=str),
                  workload_ids=features.workload_id.to_numpy(dtype=str))
    for policy, fields in [('metadata', META_FEATURES), ('pilot', PILOT_FEATURES)]:
        x = features[fields].to_numpy()
        floor = np.array([.1 if f.startswith('log_') else .01 for f in fields])
        scale = np.maximum(np.quantile(x, .75, axis=0) - np.quantile(x, .25, axis=0), floor)
        scalings[policy] = dict(fields=fields, scale=scale.tolist(), median=np.median(x, axis=0).tolist())
        arrays[policy + '_features'] = x
    baselines, baseline_rows = {}, []
    for k in CONFIG['k']:
        for kind in sorted(features.kind.unique()):
            use = (features.k == k) & (features.kind == kind)
            average = nanmean(scores[use], axis=0)
            max_bytes = np.max(scores[use, ..., STATISTICS.index('max_payload')], axis=0)
            for b, budget in enumerate(CONFIG['pair_budgets']):
                for m, metric in enumerate(METRICS):
                    for family in ['fixed'] + ALGORITHMS:
                        options = []
                        fallback = []
                        for i, (alg, p, h) in enumerate(native_grid()):
                            if family != 'fixed' and alg != family: continue
                            if not capability(alg, metric): continue
                            fallback.append(i)
                            if max_bytes[b, i, m] <= budget: options.append(i)
                        # An unavailable family remains explicitly unsupported.
                        candidates = options or fallback
                        selected = min(candidates, key=lambda i: (average[b, i, m, 0], cfgs[i])) if candidates else None
                        record = dict(k=k, kind=kind, budget=budget, metric=metric, policy=family,
                            config_index=selected, config=cfgs[selected] if selected is not None else '',
                            calibration_budget_fits=bool(options),
                            statistics=average[b, selected, m].tolist() if selected is not None else [1., None, None, 0., 0., 0., None])
                        record['statistics'] = [float(v) if v is not None and np.isfinite(v) else None for v in record['statistics']]
                        baselines[f'{k}/{kind}/{budget}/{metric}/{family}'] = record
                        baseline_rows.append({k: v for k, v in record.items() if k != 'statistics'} | {
                            'predicted_' + s: v for s, v in zip(STATISTICS, record['statistics'])})
    np.savez_compressed(MODEL / 'calibration.npz', **arrays)
    metadata = dict(schema_version=1, configurations=[dict(algorithm=a, parameter=p, half_subk=h, config=c)
        for (a, p, h), c in zip(native_grid(), cfgs)], metrics=METRICS, statistics=STATISTICS,
        budgets=CONFIG['pair_budgets'], k=CONFIG['k'], scalings=scalings, baselines=baselines,
        neighbor_genera=CONFIG['neighbor_genera'], data_sha256=sha(MODEL / 'calibration.npz'),
        pilot=dict(algorithm='FastKMV', K=32, seed=919, max_bytes_per_object=256,
                   lifecycle='discard before final sketch construction; time cost not evaluated'),
        notes=['Precision-only experimental advisor', 'Known-genus candidate banks',
               'SetSketch requires isolated corrected implementation', 'Probability Jaccard and order similarity are separate objectives'])
    atomic(MODEL / 'model.json', metadata)
    save_csv(OUT / 'development_fixed_choices.csv', pd.DataFrame(baseline_rows))
    files = [MODEL / 'model.json', MODEL / 'calibration.npz', OUT / 'task_thresholds.json',
             OUT / 'development_features.csv', OUT / 'development_fixed_choices.csv',
             HERE / 'advisor.py', HERE / 'measure.py', HERE / 'scenario_common.py', HERE / 'config.json',
             HERE / 'PROTOCOL.md', HERE / 'ANALYSIS_ADDENDUM.md', HERE / 'validation_extra_conditions.json',
             HERE / 'extend_validation.py', DATA / 'validation_extension_complete.json',
             DATA / 'validation_preparation.json', DATA / 'validation_objects.tsv',
             DATA / 'validation_pairs.tsv', DATA / 'validation_workloads.json', RUN / 'input_sha256.json', BINARY]
    atomic(RUN / 'model_frozen.json', dict(frozen_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        files={str(p): sha(p) for p in files}, development_genera=features.genus.nunique(),
        development_workloads=features.workload_id.nunique(), validation_outcomes_started=False))
    print('Frozen model', MODEL, flush=True)


class Advisor:
    def __init__(self, directory=MODEL):
        directory = Path(directory); self.metadata = json.loads((directory / 'model.json').read_text())
        if sha(directory / 'calibration.npz') != self.metadata['data_sha256']: raise ValueError('model data changed')
        with np.load(directory / 'calibration.npz', allow_pickle=False) as values:
            self.data = {k: values[k] for k in values.files}

    def neighbors(self, features, k, kind, policy):
        scaling = self.metadata['scalings'][policy]
        return nearest_distinct(features, self.data[policy + '_features'], scaling['fields'],
            np.array(scaling['scale']), k, kind, self.data['groups'], self.data['train_k'], self.data['train_kind'],
            count=self.metadata['neighbor_genera'])

    def choose(self, features, k, kind, metric, budget, policy='pilot', neighbors=None):
        if metric not in self.metadata['metrics']:
            raise ValueError('unsupported estimand; weighted probability Jaccard and order similarity require separate evidence')
        if k not in self.metadata['k'] or budget not in self.metadata['budgets']:
            raise ValueError('k or budget outside the validated grid')
        b = self.metadata['budgets'].index(budget); m = self.metadata['metrics'].index(metric)
        cfgs = self.metadata['configurations']
        if policy not in ['pilot', 'metadata']:
            baseline = self.metadata['baselines'].get(f'{k}/{kind}/{budget}/{metric}/{policy}')
            if baseline is None: raise ValueError('unsupported input kind/policy')
            i = baseline['config_index']
            result = dict(config=cfgs[i]['config'] if i is not None else '', supported=bool(i is not None),
                mean_neighbor_distance=None, outside_features='', neighbors=0, candidates=cfgs[i]['config'] if i is not None else '',
                payload_guard_bytes=None, status='calibration_fixed' if i is not None else 'unsupported')
            result.update({'predicted_' + s: v for s, v in zip(STATISTICS, baseline['statistics'])})
            return result
        nearest, supported, distance, outside = neighbors if neighbors is not None else self.neighbors(features, k, kind, policy)
        prediction = nanmean(self.data['scores'][nearest, b, :, m, :], axis=0)
        options, guards = [], {}
        for i, c in enumerate(cfgs):
            if not capability(c['algorithm'], metric): continue
            guard = payload_guard(features, k, c['algorithm'], c['parameter'], metric, policy)
            if guard <= budget:
                options.append(i); guards[i] = guard
        if not options: raise ValueError('no configuration within the observable payload guard')
        chosen = min(options, key=lambda i: (prediction[i, 0], cfgs[i]['config']))
        loss = float(prediction[chosen, 0]); tolerance = max(.0001, .1 * loss)
        candidates = [cfgs[i]['config'] for i in options if prediction[i, 0] <= loss + tolerance]
        result = dict(config=cfgs[chosen]['config'], supported=bool(supported),
            mean_neighbor_distance=distance, outside_features=';'.join(outside), neighbors=len(nearest),
            candidates=';'.join(candidates), payload_guard_bytes=guards[chosen],
            status='candidate' if supported else 'out_of_support_candidate')
        result.update({'predicted_' + s: float(v) if np.isfinite(v) else None
                       for s, v in zip(STATISTICS, prediction[chosen])})
        return result


def freeze_decisions():
    if (RUN / 'validation-outcomes_plan.json').exists(): raise ValueError('cannot change choices after validation outcomes start')
    frozen = json.loads((RUN / 'model_frozen.json').read_text())
    for p, value in frozen['files'].items():
        if sha(p) != value: raise ValueError('frozen development artifact changed')
    features = pd.read_csv(OUT / 'validation_features.csv')
    model = Advisor(); rows = []
    for f in features.to_dict('records'):
        neighbors = {policy: model.neighbors(f, int(f['k']), f['kind'], policy) for policy in ['pilot', 'metadata']}
        for budget in CONFIG['pair_budgets']:
            for metric in METRICS:
                for policy in ['pilot', 'metadata', 'fixed'] + ALGORITHMS:
                    choice = model.choose(f, int(f['k']), f['kind'], metric, budget, policy, neighbors.get(policy))
                    rows.append(dict(workload_id=f['workload_id'], genus=f['genus'], condition=f['condition'],
                        kind=f['kind'], k=f['k'], metric=metric, budget=budget, policy=policy, **choice))
    target = OUT / 'validation_decisions.csv.gz'; save_csv(target, pd.DataFrame(rows))
    atomic(RUN / 'decisions_frozen.json', dict(frozen_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        rows=len(rows), files={str(p): sha(p) for p in [target, OUT / 'validation_features.csv']},
        model_freeze_sha256=sha(RUN / 'model_frozen.json'), validation_outcomes_started=False))
    print('Frozen validation decisions:', len(rows), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage', choices=['fit', 'freeze-decisions', 'recommend'], required=True)
    ap.add_argument('--model', type=Path, default=MODEL)
    ap.add_argument('--features-json', type=Path)
    ap.add_argument('--metric', default='jaccard')
    ap.add_argument('--budget', type=int, default=8192)
    ap.add_argument('--k', type=int, default=24)
    ap.add_argument('--kind', choices=['assembly', 'contig', 'reads'], default='assembly')
    args = ap.parse_args()
    if args.stage == 'fit': fit()
    elif args.stage == 'freeze-decisions': freeze_decisions()
    else:
        if args.features_json is None: ap.error('--features-json is required for recommend')
        result = Advisor(args.model).choose(json.loads(args.features_json.read_text()), args.k, args.kind,
                                            args.metric, args.budget)
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__': main()
