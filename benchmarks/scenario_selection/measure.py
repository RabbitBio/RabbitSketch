#!/usr/bin/env python3
"""Extract observable pilot features and score native outputs without hidden failures."""
import argparse
import importlib.util
import sys
import time

from scenario_common import *

# Reuse the previously tested ARI and linkage implementations without invoking
# their old pipeline or changing the frozen Oct30 files.
sys.path.insert(0, str(HERE.parent / 'comprehensive_accuracy'))
spec = importlib.util.spec_from_file_location('oct30_downstream_math', HERE.parent / 'comprehensive_accuracy/downstream.py')
math_module = importlib.util.module_from_spec(spec); spec.loader.exec_module(math_module)
sys.path.pop(0)
linkage_curve = math_module.linkage_curve
adjusted_rand = math_module.adjusted_rand
single_linkage = math_module.single_linkage
METRICS = ['jaccard', 'containment', 'cardinality', 'neighbor', 'species', 'clustering']
STATISTICS = ['loss', 'mae', 'p95', 'availability', 'budget_violation', 'max_payload', 'zero_fraction']
META_FEATURES = ['log_query_bases_median', 'log_query_bases_max', 'log_reference_bases_max',
                 'log_query_records_median', 'log_query_record_length_median']
PILOT_FEATURES = META_FEATURES + ['log_pilot_query_n_median', 'log_pilot_query_n_max',
    'log_pilot_reference_n_max', 'log_pilot_ratio_median', 'log_pilot_ratio_p90',
    'pilot_j_mean', 'pilot_j_p90', 'pilot_c_mean', 'pilot_c_p90', 'pilot_zero_fraction',
    'pilot_top_margin_mean']
THRESHOLDS = np.arange(1, 100) / 100
WAIT_FOR_JOBS = False


def plan_jobs(split, purpose):
    if split == 'development': name = 'development'
    else: name = 'validation-pilot' if purpose == 'pilot' else 'validation-outcomes'
    plan = json.loads((RUN / (name + '_plan.json')).read_text())
    return [j for j in plan['jobs'] if j['purpose'] == purpose]


def load_job(job):
    directory = RUN / 'native' / job['job_id']
    started = time.monotonic()
    while WAIT_FOR_JOBS and not (directory / 'completion.json').exists():
        if time.monotonic() - started > 12 * 3600: raise TimeoutError('native job did not finish')
        time.sleep(2)
    completion = json.loads((directory / 'completion.json').read_text())
    raw = directory / 'raw.jsonl.gz'
    if sha(raw) != completion['raw_sha256']: raise ValueError('native archive hash mismatch')
    pairs, objects, counts = [], [], {}
    with gzip.open(raw, 'rt') as f:
        for line in f:
            row = json.loads(line); typ = row['type']; counts[typ] = counts.get(typ, 0) + 1
            if typ in ['pair', 'object']:
                row = {k: v for k, v in row.items() if 'seconds' not in k and 'rss' not in k and k != 'input_bytes'}
                (pairs if typ == 'pair' else objects).append(row)
    if counts != completion['rows']: raise ValueError('native archive counts changed')
    p, o = pd.DataFrame(pairs), pd.DataFrame(objects)
    if p.pair_id.duplicated().any() or o.object_id.duplicated().any(): raise ValueError('duplicate native observation')
    return p.set_index('pair_id'), o.set_index('object_id')


def normalized_digest_check(objects, expected):
    if set(objects.index) != set(expected.index): raise ValueError('native object inventory changed')
    a = objects.input_digest.to_numpy(dtype=np.uint64)
    b = expected.loc[objects.index, 'input_digest'].to_numpy(dtype=np.uint64)
    if not np.array_equal(a, b): raise ValueError('normalized inputs differ between native jobs')


def observable_features(workload, metadata, pairs, objects):
    """Only the explicit metadata and pilot columns below are ever accessed."""
    queries, refs = workload['queries'], workload['references']
    q, r = metadata.loc[queries], metadata.loc[refs]
    p = pairs.loc[workload['retrieval_pairs']]
    qa = objects.loc[p.a, 'cardinality'].to_numpy(dtype=float)
    rb = objects.loc[p.b, 'cardinality'].to_numpy(dtype=float)
    j = p.jaccard.to_numpy(dtype=float)
    c = bounded_containment(j, qa, rb)
    ratio = np.maximum(qa, rb) / np.minimum(qa, rb)
    matrix = j.reshape(len(queries), len(refs))
    ordered = np.sort(matrix, axis=1)
    qn = objects.loc[queries, 'cardinality'].to_numpy(dtype=float)
    rn = objects.loc[refs, 'cardinality'].to_numpy(dtype=float)
    log = lambda x: float(np.log10(max(float(x), 1)))
    values = dict(log_query_bases_median=log(q.bases.median()), log_query_bases_max=log(q.bases.max()),
        log_reference_bases_max=log(r.bases.max()), log_query_records_median=log(q.records.median()),
        log_query_record_length_median=log(q.median_record_length.median()),
        log_pilot_query_n_median=log(np.median(qn)), log_pilot_query_n_max=log(qn.max()),
        log_pilot_reference_n_max=log(rn.max()), log_pilot_ratio_median=log(np.median(ratio)),
        log_pilot_ratio_p90=log(np.quantile(ratio, .9)), pilot_j_mean=float(j.mean()),
        pilot_j_p90=float(np.quantile(j, .9)), pilot_c_mean=float(np.mean(c)),
        pilot_c_p90=float(np.quantile(c, .9)), pilot_zero_fraction=float(np.mean(j == 0)),
        pilot_top_margin_mean=float(np.mean(ordered[:, -1] - ordered[:, -2])),
        pilot_payload_bytes=int(objects.loc[queries + refs, 'payload_bytes'].sum()))
    if not all(np.isfinite(v) for v in values.values()): raise ValueError('pilot features unavailable')
    return values


def extract_features(split):
    workloads = json.loads((DATA / (split + '_workloads.json')).read_text())
    metadata = pd.read_csv(DATA / (split + '_objects.tsv'), sep='\t').set_index('object_id')
    rows = []
    for job in plan_jobs(split, 'pilot'):
        p, o = load_job(job)
        save_csv(OUT / split / f'pilot_k{job["k"]}_pairs.csv.gz', p.reset_index())
        save_csv(OUT / split / f'pilot_k{job["k"]}_objects.csv.gz', o.reset_index())
        for w in workloads:
            rows.append(dict(workload_id=w['workload_id'], split=split, genus=w['genus'],
                condition=w['condition'], kind=w['kind'], k=job['k'], **observable_features(w, metadata, p, o)))
    frame = pd.DataFrame(rows)
    if len(frame) != len(workloads) * len(CONFIG['k']): raise ValueError('missing pilot workload')
    save_csv(OUT / (split + '_features.csv'), frame)
    return frame


def square(values, n):
    if len(values) != n * (n - 1) // 2: raise ValueError('missing query-query edge')
    out = np.eye(n); i, j = np.triu_indices(n, 1); out[i, j] = values; out[j, i] = values
    return out


def truth_workloads(split, k, jobs, workloads, metadata):
    job = next(j for j in jobs if j['k'] == k and j['algorithm'] == 'exact')
    p, o = load_job(job)
    pilot = pd.read_csv(OUT / split / f'pilot_k{k}_objects.csv.gz', dtype={'input_digest': 'uint64'}).set_index('object_id')
    normalized_digest_check(o, pilot)
    p['truth_j'] = p.intersection / (p.n_a + p.n_b - p.intersection)
    p['truth_c'] = p.intersection / p.n_a
    save_csv(OUT / split / f'exact_k{k}_pairs.csv.gz', p.reset_index())
    save_csv(OUT / split / f'exact_k{k}_objects.csv.gz', o.reset_index())
    result = {}
    for w in workloads:
        qr = p.loc[w['retrieval_pairs']]; qq = p.loc[w['cluster_pairs']]
        result[w['workload_id']] = dict(j=qr.truth_j.to_numpy().reshape(9, 6),
            c=qr.truth_c.to_numpy().reshape(9, 6), n=o.loc[w['queries'], 'n'].to_numpy(),
            qq=square(qq.truth_j.to_numpy(), 9),
            query_species=metadata.loc[w['queries'], 'species'].to_numpy(),
            reference_species=metadata.loc[w['references'], 'species'].to_numpy())
    return result, o


def native_workloads(pairs, objects, workloads, algorithm, seed):
    records = []
    for w in workloads:
        qr = pairs.loc[w['retrieval_pairs']]; qq = pairs.loc[w['cluster_pairs']]
        j = qr.jaccard.to_numpy(dtype=float).reshape(9, 6)
        ca = objects.loc[qr.a, 'cardinality'].to_numpy(dtype=float)
        cb = objects.loc[qr.b, 'cardinality'].to_numpy(dtype=float)
        c = bounded_containment(qr.jaccard, ca, cb).reshape(9, 6)
        available = ((qr.empty_a.to_numpy() == 0) & (qr.empty_b.to_numpy() == 0)).reshape(9, 6)
        available &= np.isfinite(j) & (j >= 0) & (j <= 1)
        qqvalues = qq.jaccard.to_numpy(dtype=float)
        qqvalid = (np.isfinite(qqvalues) & (qqvalues >= 0) & (qqvalues <= 1) &
                   (qq.empty_a.to_numpy() == 0) & (qq.empty_b.to_numpy() == 0))
        cards = objects.loc[w['queries'], 'cardinality'].to_numpy(dtype=float)
        cardgood = np.isfinite(cards) & (cards > 0) & (objects.loc[w['queries'], 'empty'].to_numpy() == 0)
        if algorithm in ['kssd', 'probminhash_unweighted']:
            c[:] = np.nan; cardgood[:] = False
        records.append(dict(workload=w, seed=seed, j=j, c=c, available=available,
            pair_bytes=(qr.payload_bytes_a + qr.payload_bytes_b).to_numpy().reshape(9, 6),
            qq=square(qqvalues, 9), qqvalid=bool(qqvalid.all()),
            qqbytes=(qq.payload_bytes_a + qq.payload_bytes_b).to_numpy(),
            cards=cards, cardgood=cardgood,
            cardbytes=2 * objects.loc[w['queries'], 'payload_bytes'].to_numpy()))
    return records


def fit_task_thresholds(records, truths):
    result = {}
    for kind in sorted({r['workload']['kind'] for r in records}):
        rs = [r for r in records if r['workload']['kind'] == kind]
        # Every workload has the same nine queries; average seeds within a
        # workload before averaging genera/conditions.
        curves, tops = {}, {}
        for r in rs:
            wid = r['workload']['workload_id']; t = truths[wid]
            if r['qqvalid']:
                curves.setdefault(wid, []).append(linkage_curve(r['qq'], t['query_species'], THRESHOLDS))
            valid_queries = r['available'].all(axis=1)
            tops.setdefault(wid, []).extend(r['j'].max(axis=1)[valid_queries].tolist())
        valid_curves = [np.mean(v, axis=0) for v in curves.values()]
        threshold = float(THRESHOLDS[int(np.argmax(np.mean(valid_curves, axis=0)))]) if valid_curves else None
        scores = np.array([v for values in tops.values() for v in values])
        accept = float(np.quantile(scores, .1, method='lower')) if len(scores) else None
        result[kind] = dict(clustering_threshold=threshold, acceptance_threshold=accept,
            clustering_workloads=len(curves), acceptance_queries=len(scores),
            target_development_known_recall=.9)
    return result


def statistics(errors, available, budget_bad, payload, zeros=None, clip=False):
    errors = np.asarray(errors, dtype=float).ravel()
    available = np.asarray(available, dtype=bool).ravel() & np.isfinite(errors)
    budget_bad = np.broadcast_to(budget_bad, errors.shape).ravel()
    good = available & ~budget_bad
    loss = np.where(good, np.minimum(errors, 1) if clip else errors, 1.)
    values = errors[good]
    return dict(loss=float(np.mean(loss)), mae=float(np.mean(values)) if len(values) else np.nan,
        p95=float(np.quantile(values, .95)) if len(values) else np.nan,
        availability=float(np.mean(good)), budget_violation=float(np.mean(budget_bad)),
        max_payload=float(np.max(payload)), zero_fraction=float(np.mean(zeros)) if zeros is not None else np.nan)


def score_group(records, truths, thresholds, budget):
    errors = {m: [] for m in METRICS}; avail = {m: [] for m in METRICS}
    violations = {m: [] for m in METRICS}; payload = {m: [] for m in METRICS}; zeros = []
    known, unknown, retrieval_valid, unknown_valid, aris, oracle_species = [], [], [], [], [], []
    truths_for_diagnostics = {'jaccard': [], 'containment': []}
    for r in records:
        w = r['workload']; t = truths[w['workload_id']]; fit = thresholds[w['kind']]
        bad = r['pair_bytes'] > budget
        for metric, error, valid in [('jaccard', abs(r['j'] - t['j']), r['available']),
                ('containment', abs(r['c'] - t['c']), r['available'] & np.isfinite(r['c']))]:
            errors[metric].extend(error.ravel()); avail[metric].extend(valid.ravel())
            violations[metric].extend(bad.ravel()); payload[metric].extend(r['pair_bytes'].ravel())
            truths_for_diagnostics[metric].extend(t['j' if metric == 'jaccard' else 'c'].ravel())
        errors['cardinality'].extend(abs(r['cards'] / t['n'] - 1)); avail['cardinality'].extend(r['cardgood'])
        violations['cardinality'].extend(r['cardbytes'] > budget); payload['cardinality'].extend(r['cardbytes'])
        zeros.extend((r['j'] == 0).ravel())
        valid_query = r['available'].all(axis=1)
        query_budget_bad = bad.any(axis=1)
        for i in range(9):
            top = r['j'][i] == np.nanmax(r['j'][i])
            exact_top = t['j'][i] == t['j'][i].max()
            species_matches = t['reference_species'] == t['query_species'][i]
            nn_credit = float(np.mean(exact_top[top])) if top.any() else 0.
            species_credit = float(np.mean(species_matches[top])) if top.any() else 0.
            for metric, error in [('neighbor', 1 - nn_credit), ('species', 1 - species_credit)]:
                errors[metric].append(error); avail[metric].append(valid_query[i]);
                violations[metric].append(query_budget_bad[i]); payload[metric].append(r['pair_bytes'][i].max())
            threshold = fit['acceptance_threshold']
            good = bool(valid_query[i] and not query_budget_bad[i] and threshold is not None)
            unknown_good = bool(r['available'][i][~species_matches].all() and
                not bad[i][~species_matches].any() and threshold is not None)
            known.append(float(good and r['j'][i].max() >= threshold))
            unknown.append(float(unknown_good and r['j'][i][~species_matches].max() >= threshold))
            retrieval_valid.append(float(good))
            unknown_valid.append(float(unknown_good))
            oracle_species.append(float(np.mean(species_matches[exact_top])))
        ct = fit['clustering_threshold']
        good = r['qqvalid'] and ct is not None
        ari = adjusted_rand(t['query_species'], single_linkage(r['qq'], ct)) if good else np.nan
        errors['clustering'].append((1 - ari) / 2); avail['clustering'].append(good)
        violations['clustering'].append(bool((r['qqbytes'] > budget).any())); payload['clustering'].append(r['qqbytes'].max())
        aris.append(ari if good and not (r['qqbytes'] > budget).any() else np.nan)
    result = {}
    for metric in METRICS:
        result[metric] = statistics(errors[metric], avail[metric], violations[metric], payload[metric],
                                    zeros if metric in ['jaccard', 'containment'] else None, clip=metric == 'cardinality')
        if metric in truths_for_diagnostics:
            truth = np.asarray(truths_for_diagnostics[metric]); error = np.asarray(errors[metric])
            good = np.asarray(avail[metric]) & ~np.asarray(violations[metric]) & np.isfinite(error) & (truth > 0)
            relative = error[good] / truth[good]
            result[metric].update(constant0_mae=float(truth.mean()), constant1_mae=float((1-truth).mean()),
                truth_mean=float(truth.mean()), positive_truth_fraction=float(np.mean(truth > 0)),
                positive_pair_relative_mae=float(relative.mean()) if len(relative) else np.nan)
    extra = dict(known_acceptance=float(np.mean(known)), unknown_false_acceptance=float(np.mean(unknown)),
        open_availability=float(np.mean(retrieval_valid)), ari=float(np.nanmean(aris)) if np.isfinite(aris).any() else np.nan,
        unknown_availability=float(np.mean(unknown_valid)),
        exact_species_accuracy=float(np.mean(oracle_species)))
    return result, extra


def exact_records(workloads, truths):
    records = []
    for w in workloads:
        t = truths[w['workload_id']]
        records.append(dict(workload=w, seed=0, j=t['j'], c=t['c'], available=np.ones((9, 6), bool),
            pair_bytes=np.zeros((9, 6)), qq=t['qq'], qqvalid=True, qqbytes=np.zeros(36),
            cards=t['n'], cardgood=np.ones(9, bool), cardbytes=np.zeros(9)))
    return records


def score_panel(split):
    jobs = plan_jobs(split, 'outcomes')
    workloads = json.loads((DATA / (split + '_workloads.json')).read_text())
    metadata = pd.read_csv(DATA / (split + '_objects.tsv'), sep='\t').set_index('object_id')
    threshold_file = OUT / 'task_thresholds.json'
    all_thresholds = json.loads(threshold_file.read_text()) if threshold_file.exists() else {}
    for k in CONFIG['k']:
        truths, exact_objects = truth_workloads(split, k, jobs, workloads, metadata)
        configurations = [('exact', 1, 6)] + native_grid()
        for number, (algorithm, parameter, half) in enumerate(configurations, 1):
            cfg = configuration(algorithm, parameter, half)
            target = OUT / split / 'scores' / f'k{k}_{algorithm}_{parameter}_{half}.csv.gz'
            detail_target = OUT / split / 'tasks' / target.name
            done = target.with_suffix('.complete.json')
            if done.exists():
                identity = json.loads(done.read_text())
                if sha(__file__) != identity['measure_sha256']: raise ValueError('scoring code changed after completed panel')
                if sha(target) != identity['scores_sha256'] or sha(detail_target) != identity['tasks_sha256']:
                    raise ValueError('scored panel changed')
                continue
            if algorithm == 'exact': records = exact_records(workloads, truths)
            else:
                chosen = [j for j in jobs if j['k'] == k and j['algorithm'] == algorithm and
                          j['parameter'] == parameter and j['half_subk'] == half]
                records = []
                for job in chosen:
                    p, o = load_job(job); normalized_digest_check(o, exact_objects)
                    records.extend(native_workloads(p, o, workloads, algorithm, job['seed']))
                if not chosen: raise ValueError('missing configuration')
            threshold_key = f'{k}/{cfg}'
            if split == 'development':
                all_thresholds[threshold_key] = fit_task_thresholds(records, truths)
                atomic(threshold_file, all_thresholds)
            elif threshold_key not in all_thresholds: raise ValueError('no frozen task thresholds')
            by_workload = {}
            for record in records: by_workload.setdefault(record['workload']['workload_id'], []).append(record)
            rows, details = [], []
            for w in workloads:
                rs = by_workload[w['workload_id']]
                for budget in CONFIG['pair_budgets']:
                    scores, extra = score_group(rs, truths, all_thresholds[threshold_key], budget)
                    base = dict(workload_id=w['workload_id'], split=split, genus=w['genus'],
                        condition=w['condition'], kind=w['kind'], k=k, budget=budget,
                        algorithm=algorithm, parameter=parameter, half_subk=half, config=cfg, seeds=len(rs))
                    rows.extend(dict(base, metric=m, **v) for m, v in scores.items())
                    details.append(dict(base, **extra))
            save_csv(target, pd.DataFrame(rows)); save_csv(detail_target, pd.DataFrame(details))
            atomic(done, dict(measure_sha256=sha(__file__), scores_sha256=sha(target), tasks_sha256=sha(detail_target)))
            print(split, 'scored k', k, number, '/', len(configurations), cfg, flush=True)
    atomic(OUT / (split + '_scoring_complete.json'), dict(measure_sha256=sha(__file__),
        task_math_sha256=sha(HERE.parent / 'comprehensive_accuracy/downstream.py'),
        configurations=(len(native_grid()) + 1) * len(CONFIG['k'])))


def main():
    global WAIT_FOR_JOBS
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--split', choices=['development', 'validation'], required=True)
    ap.add_argument('--stage', choices=['features', 'scores'], required=True)
    ap.add_argument('--wait', action='store_true', help='wait for registered native jobs to finish')
    args = ap.parse_args()
    WAIT_FOR_JOBS = args.wait
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == 'features': extract_features(args.split)
    else: score_panel(args.split)


if __name__ == '__main__': main()
