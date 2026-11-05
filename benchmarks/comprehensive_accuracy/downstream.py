#!/usr/bin/env python3
"""Species retrieval and calibrated single-linkage clustering on held-out genera."""
import itertools

from common import *


def adjusted_rand(labels, clusters):
    labels, clusters = np.asarray(labels), np.asarray(clusters)
    n = len(labels)
    if n != len(clusters) or n < 2: raise ValueError('ARI needs matching nontrivial partitions')
    _, a = np.unique(labels, return_inverse=True)
    _, b = np.unique(clusters, return_inverse=True)
    table = np.zeros((a.max() + 1, b.max() + 1), dtype=int)
    np.add.at(table, (a, b), 1)
    choose2 = lambda x: np.asarray(x) * (np.asarray(x) - 1) / 2
    observed = choose2(table).sum()
    rows, columns = choose2(table.sum(axis=1)).sum(), choose2(table.sum(axis=0)).sum()
    expected = rows * columns / choose2(n)
    maximum = (rows + columns) / 2
    return 1. if maximum == expected else (observed - expected) / (maximum - expected)


def single_linkage(scores, threshold):
    scores = np.asarray(scores)
    n = len(scores); parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a in range(n):
        for b in range(a):
            if scores[a, b] >= threshold:
                ra, rb = find(a), find(b); parent[ra] = rb
    return np.array([find(i) for i in range(n)])


def linkage_curve(scores, labels, thresholds):
    """Reuse nested single-linkage partitions across the registered thresholds."""
    n = len(scores); parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    edges = sorted(((scores[i, j], i, j) for i in range(n) for j in range(i)), reverse=True)
    edge = 0; previous = None; score = np.nan; result = np.empty(len(thresholds))
    for t in np.argsort(thresholds)[::-1]:
        threshold = thresholds[t]
        while edge < len(edges) and edges[edge][0] >= threshold:
            _, i, j = edges[edge]; parent[find(i)] = find(j); edge += 1
        partition = tuple(find(i) for i in range(n))
        if partition != previous:
            score = adjusted_rand(labels, partition); previous = partition
        result[t] = score
    return result


def nearest_scores(scores, labels, exact_scores):
    """Fractional expected credit under uniform tie-breaking; self is excluded."""
    scores, exact_scores = np.array(scores, copy=True), np.array(exact_scores, copy=True)
    np.fill_diagonal(scores, -np.inf); np.fill_diagonal(exact_scores, -np.inf)
    species, recovery, ties = [], [], []
    labels = np.asarray(labels)
    for i in range(len(scores)):
        maximum = scores[i].max(); truth_maximum = exact_scores[i].max()
        # Exact equality keeps the same rule for discrete and continuous backends.
        top = scores[i] == maximum; exact_top = exact_scores[i] == truth_maximum
        species.append(np.mean(labels[top] == labels[i]))
        recovery.append(np.mean(exact_top[top])); ties.append(top.sum())
    return float(np.mean(species)), float(np.mean(recovery)), float(np.mean(ties))


def matrix(pairs, ids, column):
    n = len(ids); lookup = {name: i for i, name in enumerate(ids)}
    if len(pairs) != n * (n - 1) // 2: raise ValueError('missing within-genus edges')
    out = np.full((n, n), np.nan); np.fill_diagonal(out, 1.)
    seen = set()
    for row in pairs[['a', 'b', column]].itertuples(index=False, name=None):
        a, b, value = row; i, j = lookup[a], lookup[b]
        edge = tuple(sorted([i, j]))
        if i == j or edge in seen: raise ValueError('duplicate or diagonal edge')
        seen.add(edge); out[i, j] = out[j, i] = value
    if not np.isfinite(out).all(): raise ValueError('unavailable edge in distance matrix')
    return out


def calibrate_threshold(curves, splits, thresholds):
    calibration = np.asarray(curves)[np.asarray(splits) == 'calibration']
    if len(calibration) == 0 or not np.isfinite(calibration).all():
        raise ValueError('threshold calibration has missing genera')
    mean = calibration.mean(axis=0)
    # Registered threshold order resolves exact ties deterministically.
    return int(np.argmax(mean))


def analyze(k):
    directory = OUT / 'real_panels' / f'natural_k{k}'
    d = pd.read_csv(directory / 'selected_rows.csv.gz', low_memory=False)
    truth = pd.read_csv(directory / 'exact_truth.csv')
    taxonomy = pd.read_csv(RUN / 'data/natural_taxonomy.csv')
    truth = truth[truth.genus_a == truth.genus_b]
    d = d[d.genus_a == d.genus_b]
    genus_info, exact_matrices = {}, {}
    thresholds = np.arange(1, 100) / 100
    for genus, taxa in taxonomy.groupby('genus'):
        taxa = taxa.sort_values('object_id')
        ids, labels = taxa.object_id.tolist(), taxa.species.to_numpy()
        if len(ids) != 15 or len(set(labels)) != 3: raise ValueError('unexpected genus composition')
        e = matrix(truth[truth.genus_a == genus], ids, 'truth_j')
        genus_info[genus] = (ids, labels, taxa.split.iloc[0])
        exact_matrices[genus] = e
    rows, calibration_records = [], []
    experiments = []
    for budget in CONFIG['pair_budgets']:
        experiments.append(('exact', budget, 0, None))
        sub = d[d.budget == budget]
        for algorithm, frame in sub.groupby('algorithm'):
            for seed, group in frame.groupby('seed'):
                experiments.append((algorithm, budget, int(seed), group))
    # A threshold is calibrated jointly over the available seeds, not separately per seed.
    records_by_config = {}
    for algorithm, budget, seed, group in experiments:
        for genus, (ids, labels, split) in genus_info.items():
            base = dict(k=k, algorithm=algorithm, budget=budget, seed=seed, genus=int(genus), split=split)
            if algorithm == 'exact':
                scores = exact_matrices[genus]; valid = True
            else:
                g = group[group.genus_a == genus]
                valid = (len(g) == 105 and g.j_error.notna().all() and (g.pair_bytes <= budget).all())
                scores = matrix(g, ids, 'jaccard') if valid else None
            if valid:
                species, recovery, ties = nearest_scores(scores, labels, exact_matrices[genus])
                curve = linkage_curve(scores, labels, thresholds)
            else:
                species = recovery = ties = np.nan; curve = np.full(len(thresholds), np.nan)
            record = dict(base, valid=int(valid), species_accuracy=species,
                          exact_neighbor_recovery=recovery, mean_tie_size=ties, curve=curve)
            records_by_config.setdefault((algorithm, budget), []).append(record)
    for (algorithm, budget), records in records_by_config.items():
        cal = [r for r in records if r['split'] == 'calibration']
        calibration_valid = all(r['valid'] for r in cal)
        if calibration_valid:
            chosen = calibrate_threshold(np.array([r['curve'] for r in records]),
                                         [r['split'] for r in records], thresholds)
            threshold = float(thresholds[chosen])
            calibration_ari = np.mean([r['curve'][chosen] for r in cal])
        else:
            chosen = None; threshold = calibration_ari = np.nan
        calibration_records.append(dict(k=k, algorithm=algorithm, budget=budget,
            threshold=threshold, calibration_ari=calibration_ari,
            calibration_complete=int(calibration_valid), calibration_genus_seed_rows=len(cal),
            calibration_valid_rows=sum(r['valid'] for r in cal)))
        for r in records:
            curve = r.pop('curve')
            r.update(ari=curve[chosen] if chosen is not None else np.nan, threshold=threshold,
                     calibration_complete=int(calibration_valid))
            rows.append(r)
    return pd.DataFrame(rows), pd.DataFrame(calibration_records)


def summarize(rows, calibrations):
    # Repeated seeds are nested in a genus, not independent biological units.
    keys = ['k', 'budget', 'algorithm', 'genus', 'split']
    per_genus = rows.groupby(keys).agg(species_accuracy=('species_accuracy', 'mean'),
        exact_neighbor_recovery=('exact_neighbor_recovery', 'mean'), ari=('ari', 'mean'),
        all_seeds_valid=('valid', 'min'), calibration_complete=('calibration_complete', 'min'),
        mean_tie_size=('mean_tie_size', 'mean')).reset_index()
    # An incomplete seed panel remains unavailable instead of averaging only successes.
    for metric in ['species_accuracy', 'exact_neighbor_recovery', 'ari']:
        per_genus.loc[per_genus.all_seeds_valid == 0, metric] = np.nan
    save('downstream_genus_scores', per_genus)
    records = []
    for group, g in per_genus.groupby(['k', 'budget', 'algorithm', 'split']):
        for metric in ['species_accuracy', 'exact_neighbor_recovery', 'ari']:
            valid = g[metric].dropna()
            lo, hi = bootstrap_clusters(valid.to_numpy())
            records.append(dict(zip(['k', 'budget', 'algorithm', 'split'], group), metric=metric,
                total_genera=len(g), valid_genera=len(valid), mean=valid.mean(), ci_low=lo, ci_high=hi,
                calibration_complete=int(g.calibration_complete.min())))
    summary = pd.DataFrame(records); save('downstream_summary', summary)
    differences = []
    for (k, budget), group in per_genus[per_genus.split == 'validation'].groupby(['k', 'budget']):
        for left, right in itertools.combinations(sorted(group.algorithm.unique()), 2):
            a = group[group.algorithm == left].set_index('genus')
            b = group[group.algorithm == right].set_index('genus')
            for metric in ['species_accuracy', 'exact_neighbor_recovery', 'ari']:
                pair = pd.concat([a[metric].rename('a'), b[metric].rename('b')], axis=1).dropna()
                delta = pair.a - pair.b; lo, hi = bootstrap_clusters(delta.to_numpy())
                differences.append(dict(k=k, budget=budget, metric=metric, left=left, right=right,
                    common_genera=len(pair), difference=delta.mean(), ci_low=lo, ci_high=hi,
                    scope='pointwise paired genus bootstrap; no multiplicity correction'))
    save('downstream_paired_intervals', pd.DataFrame(differences))
    # Choose a fixed algorithm for each task on calibration, then report untouched validation.
    choices = []
    for (k, budget, metric), group in summary.groupby(['k', 'budget', 'metric']):
        calibration = group[(group.split == 'calibration') & (group.algorithm != 'exact') &
                            (group.valid_genera == 30) & (group.calibration_complete == 1)]
        if calibration.empty:
            choices.append(dict(k=k, budget=budget, metric=metric, algorithm='no_candidate')); continue
        best = calibration.sort_values(['mean', 'algorithm'], ascending=[False, True]).iloc[0]
        validation = group[(group.split == 'validation') & (group.algorithm == best.algorithm)].iloc[0]
        exact = group[(group.split == 'validation') & (group.algorithm == 'exact')].iloc[0]
        choices.append(dict(k=k, budget=budget, metric=metric, algorithm=best.algorithm,
            calibration_mean=best['mean'], validation_mean=validation['mean'],
            validation_valid_genera=int(validation.valid_genera), exact_validation_mean=exact['mean'],
            scope='one fixed task-level choice, no per-genus oracle'))
    save('downstream_calibrated_choices', pd.DataFrame(choices))


def main():
    all_rows, all_calibration = [], []
    for k in CONFIG['real_k']:
        rows, calibration = analyze(k)
        all_rows.append(rows); all_calibration.append(calibration)
        print(f'Downstream k={k}: {len(rows)} genus/seed observations', flush=True)
    rows, calibration = pd.concat(all_rows, ignore_index=True), pd.concat(all_calibration, ignore_index=True)
    save('downstream_seed_scores', rows); save('downstream_thresholds', calibration)
    summarize(rows, calibration)
    atomic(OUT / 'downstream_complete.json', dict(status='complete', calibration_genera=30,
        validation_genera=30, genomes=900, species=180,
        source_sha256=sha(__file__), common_sha256=sha(HERE / 'common.py')))


if __name__ == '__main__': main()
