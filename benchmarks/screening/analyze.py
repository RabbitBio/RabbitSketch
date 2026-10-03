#!/usr/bin/env python3
"""Aggregate paired replicates, compute Pareto fronts, and flag review candidates."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from run import atomic_json, digest, file_sha


OBJECTIVES = ('mae', 'payload_bytes', 'build_seconds', 'query_seconds')


def dominates(a, b, keys=OBJECTIVES):
    return all(a[k] <= b[k] for k in keys) and any(a[k] < b[k] for k in keys)


def crossover(a, b):
    denominator = a['query_seconds'] - b['query_seconds']
    if denominator == 0:
        return None
    q = (b['build_seconds'] - a['build_seconds']) / denominator
    return float(q) if np.isfinite(q) and q > 0 else None


def paired_error_upper(a, b, indices):
    """Upper two-sided 95% bootstrap bound for MAE(a) - MAE(b)."""
    differences = np.asarray(a) - np.asarray(b)
    return float(np.quantile(differences[indices].mean(axis=1), 0.975))


def aggregate(plan, rows):
    profiles = {(r.get('measurement_profile', 'serial_pair'), r.get('batch_threads', 1))
                for r in rows if r['status'] == 'ok'}
    if len(profiles) > 1:
        raise ValueError('cannot mix measurement profiles or thread counts in one frontier')
    settings = plan['config']['analysis']
    expected = {j['job_id']: j for j in plan['jobs']}
    by_job = {}
    input_checksums = {}
    for row in rows:
        job_id = row['job_id']
        if job_id not in expected or job_id in by_job:
            raise ValueError('unexpected or duplicate job ID')
        job = expected[job_id]
        if any(row.get(key) != value for key, value in job.items()):
            raise ValueError('row no longer matches its registered plan')
        if row['status'] == 'ok':
            for field in ('estimate', 'sketch_payload_bytes', 'build_seconds', 'query_seconds'):
                value = row.get(field)
                if value is None or not np.isfinite(value) or value < 0:
                    raise ValueError(f'invalid successful measurement: {field}')
            if row['estimate'] > 1:
                raise ValueError('Jaccard outside [0,1]')
            pair_key = (job['regime_id'], job['replicate'])
            checksum = (row['input_digest_a'], row['input_digest_b'])
            if pair_key in input_checksums and input_checksums[pair_key] != checksum:
                raise ValueError('algorithms did not receive the same input sets')
            input_checksums[pair_key] = checksum
        by_job[job_id] = row
    groups = {}
    for job in plan['jobs']:
        groups.setdefault((job['regime_id'], job['config_id']), []).append(job)
    summaries, samples, bootstrap = [], {}, {}
    nrep = plan['config']['replicates']
    for (regime, config_id), jobs in sorted(groups.items()):
        jobs.sort(key=lambda j: j['replicate'])
        present = [by_job[j['job_id']] for j in jobs if j['job_id'] in by_job]
        good = [r for r in present if r['status'] == 'ok']
        first = jobs[0]
        summary = {k: first[k] for k in ('regime_id', 'config_id', 'algorithm', 'parameter', 'k',
                                         'n_a', 'n_b', 'size_ratio', 'target_jaccard', 'ground_truth',
                                         'intersection', 'bbits', 'half_subk')}
        empty = sum(bool(r.get('empty_sample_a') or r.get('empty_sample_b')) for r in good)
        summary.update(expected=nrep, observed=len(present), valid=len(good), missing=nrep-len(present),
                       failed=len(present)-len(good), empty_samples=empty,
                       eligible=len(good) == nrep and empty == 0,
                       error_space_frontier=False, pareto_frontier=False, robustly_dominated=False)
        if good:
            summary.update(measurement_profile=good[0].get('measurement_profile', 'serial_pair'),
                           batch_threads=good[0].get('batch_threads', 1),
                           rss_scope=good[0].get('rss_scope', 'isolated_pair_process'))
            errors = np.array([r['estimate'] - r['ground_truth'] for r in good])
            absolute = np.abs(errors)
            if regime not in bootstrap:
                seed = int(digest([settings['bootstrap_seed'], regime])[:16], 16)
                bootstrap[regime] = np.random.default_rng(seed).integers(
                    0, nrep, size=(settings['bootstrap_samples'], nrep))
            if len(good) == nrep:
                boot = absolute[bootstrap[regime]].mean(axis=1)
                lo, hi = np.quantile(boot, [0.025, 0.975])
                samples[(regime, config_id)] = absolute
            else:
                lo = hi = np.nan
            summary.update(mae=float(absolute.mean()), rmse=float(np.sqrt(np.mean(errors**2))),
                           bias=float(errors.mean()), mae_ci_lower=float(lo), mae_ci_upper=float(hi),
                           error_p95=float(np.quantile(absolute, 0.95)),
                           zero_estimate_fraction=float(np.mean([r['estimate'] == 0 for r in good])),
                           payload_bytes=float(np.mean([r['sketch_payload_bytes'] for r in good])),
                           build_seconds=float(np.median([r['build_seconds'] for r in good])),
                           setup_seconds=float(np.median([r['setup_seconds'] for r in good])),
                           query_seconds=float(np.median([r['query_seconds'] for r in good])),
                           query_cold_seconds=float(np.median([r['query_cold_seconds'] for r in good])),
                           process_peak_rss_bytes=max(r['process_peak_rss_bytes'] for r in good))
            for side in ('a', 'b'):
                cards = [abs(r[f'cardinality_estimate_{side}'] - r[f'n_{side}']) / r[f'n_{side}']
                         for r in good if r.get(f'cardinality_estimate_{side}') is not None]
                summary[f'cardinality_relative_mae_{side}'] = float(np.mean(cards)) if cards else np.nan
        summaries.append(summary)
    for regime in (r['regime_id'] for r in plan['regimes']):
        local = [s for s in summaries if s['regime_id'] == regime and s['eligible']]
        for target in local:
            dominators = [s for s in local if dominates(s, target)]
            target['pareto_frontier'] = not dominators
            target['error_space_frontier'] = not any(dominates(s, target, ('mae', 'payload_bytes')) for s in local)
            robust = []
            for other in dominators:
                if other['algorithm'] == target['algorithm']:
                    continue
                upper = paired_error_upper(samples[(regime, other['config_id'])],
                                           samples[(regime, target['config_id'])], bootstrap[regime])
                useful_gain = any(other[k] <= target[k] * (1-settings['resource_advantage_fraction'])
                                  for k in OBJECTIVES[1:])
                if upper <= settings['error_equivalence_margin'] and useful_gain:
                    robust.append(other['config_id'])
            target['robustly_dominated'] = bool(robust)
            target['dominator_config_ids'] = ';'.join(robust)
    return summaries


def decisions(plan, summaries):
    settings = plan['config']['analysis']
    result = []
    for algorithm in plan['config']['algorithms']:
        local = [s for s in summaries if s['algorithm'] == algorithm]
        frontier_count = robust_count = 0
        for regime in plan['regimes']:
            points = [s for s in local if s['regime_id'] == regime['regime_id'] and s['eligible']]
            frontier_count += any(s['pareto_frontier'] for s in points)
            robust_count += bool(points) and all(s['robustly_dominated'] for s in points)
        gains = []
        for point in (s for s in local if s['eligible']):
            others = [s for s in summaries if s['regime_id'] == point['regime_id'] and s['eligible']
                      and s['algorithm'] != algorithm]
            for epsilon in settings['error_constraints']:
                if point['mae_ci_upper'] > epsilon:
                    continue
                feasible = [s for s in others if s['mae_ci_upper'] <= epsilon]
                for resource in OBJECTIVES[1:]:
                    if not feasible or point[resource] <= min(s[resource] for s in feasible) * (1-settings['resource_advantage_fraction']):
                        gains.append((point['regime_id'], epsilon, resource))
        regimes = len(plan['regimes'])
        complete = all(s['missing'] == 0 and s['failed'] == 0 for s in summaries)
        tested_every_regime = all(any(s['regime_id'] == r['regime_id'] and s['eligible'] for s in local)
                                  for r in plan['regimes'])
        mature = (plan['config']['purpose'] != 'smoke_only' and
                  plan['config']['replicates'] >= settings['minimum_replicates'])
        candidate = (complete and tested_every_regime and mature and not gains and
                     frontier_count / regimes <= 1-settings['dominated_regime_fraction'] + 1e-12 and
                     robust_count / regimes >= settings['dominated_regime_fraction'])
        status = 'review_pruning_candidate' if candidate else 'retain_pending_evidence'
        if not mature:
            status = 'insufficient_replicates_or_smoke'
        elif not complete or not tested_every_regime:
            status = 'incomplete_or_unavailable'
        result.append(dict(algorithm=algorithm, status=status, regimes=regimes,
                           frontier_regimes=frontier_count, frontier_fraction=frontier_count/regimes,
                           robustly_dominated_regimes=robust_count,
                           meaningful_resource_advantage=bool(gains),
                           scope='general_set_jaccard_only',
                           specialized_metric_protection=algorithm in ('probminhash', 'hll', 'setsketch')))
    return result


def reuse_tables(plan, summaries):
    winners, crossovers = [], []
    for regime in plan['regimes']:
        points = [s for s in summaries if s['regime_id'] == regime['regime_id'] and s['eligible']]
        for epsilon, q in itertools.product(plan['config']['analysis']['error_constraints'],
                                            plan['config']['analysis']['query_reuse']):
            possible = [s for s in points if s['mae_ci_upper'] <= epsilon]
            if not possible:
                winners.append(dict(regime_id=regime['regime_id'], epsilon=epsilon, query_reuse=q,
                                    algorithm='no_feasible_configuration'))
                continue
            winner = min(possible, key=lambda s: s['build_seconds'] + q*s['query_seconds'])
            winners.append(dict(regime_id=regime['regime_id'], epsilon=epsilon, query_reuse=q,
                                algorithm=winner['algorithm'], parameter=winner['parameter'],
                                config_id=winner['config_id'], payload_bytes=winner['payload_bytes'],
                                predicted_seconds=winner['build_seconds'] + q*winner['query_seconds']))
        for a, b in itertools.combinations(points, 2):
            if a['algorithm'] == b['algorithm']:
                continue
            q = crossover(a, b)
            if q is not None:
                crossovers.append(dict(regime_id=regime['regime_id'], config_a=a['config_id'],
                                       config_b=b['config_id'], algorithm_a=a['algorithm'],
                                       algorithm_b=b['algorithm'], query_reuse_crossover=q,
                                       large_q_winner=a['algorithm'] if a['query_seconds'] < b['query_seconds'] else b['algorithm']))
    return winners, crossovers


def plot(summaries, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    df = pd.DataFrame(summaries)
    valid = df[df['valid'] > 0]
    if valid.empty:
        return
    colors = dict(zip(sorted(df.algorithm.unique()), plt.get_cmap('tab10').colors))
    with PdfPages(output / 'pareto_by_regime.pdf') as pdf:
        for page, (_, points) in enumerate(valid.groupby('regime_id', sort=True)):
            first = points.iloc[0]
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
            for name, group in points.groupby('algorithm', sort=True):
                for ax, column in zip(axes, ['mae', 'query_seconds']):
                    ax.scatter(group.payload_bytes, group[column], label=name, color=colors[name], s=28)
                    frontier = group[group.pareto_frontier]
                    ax.scatter(frontier.payload_bytes, frontier[column], s=65, facecolors='none',
                               edgecolors=colors[name])
            for ax in axes:
                ax.set_xscale('log'); ax.grid(alpha=.2)
                ax.set_xlabel('Pair payload (bytes; A + B)')
            axes[0].set_yscale('symlog', linthresh=1e-5)
            axes[0].set_ylim(bottom=0)
            axes[0].set_ylabel('Jaccard MAE')
            axes[1].set_yscale('log')
            axes[1].set_ylabel('Warm comparison time (s/pair, batch load)' if
                              first.get('measurement_profile') == 'main_openmp_batch' else
                              'Warm comparison time (s/pair)')
            axes[1].legend(fontsize=8)
            fig.suptitle(f'k={first.k}; |A|={first.n_a:,}; |B|={first.n_b:,}; '
                         f'J={first.ground_truth:.5g}; n={first.expected}; exploratory frontier')
            fig.tight_layout(); pdf.savefig(fig)
            if page == 0:
                fig.savefig(output / 'example_regime.png', dpi=170)
            plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('run', type=Path)
    p.add_argument('--no-plots', action='store_true')
    args = p.parse_args()
    plan = json.loads((args.run / 'plan.json').read_text())
    rows = [json.loads(path.read_text()) for path in sorted((args.run / 'rows').glob('*.json'))]
    summary = aggregate(plan, rows)
    decision = decisions(plan, summary)
    winners, crossovers = reuse_tables(plan, summary)
    output = args.run / 'analysis'; output.mkdir(exist_ok=True)
    for name, data in [('configurations', summary), ('algorithm_decisions', decision),
                       ('query_reuse_winners', winners), ('crossovers', crossovers)]:
        pd.DataFrame(data).to_csv(output / (name + '.csv'), index=False)
    batches = [json.loads(path.read_text()) for path in sorted((args.run / 'batches').glob('*.json'))]
    if batches:
        pd.DataFrame(batches).to_csv(output / 'batches.csv', index=False)
    profiles = sorted({(r.get('measurement_profile', 'serial_pair'), r.get('batch_threads', 1))
                       for r in rows if r['status'] == 'ok'})
    atomic_json(output / 'analysis_manifest.json', dict(plan_sha256=digest(plan),
                raw_rows_sha256=digest(rows), analysis_settings=plan['config']['analysis'],
                analyzer_sha256=file_sha(__file__), numpy_version=np.__version__, pandas_version=pd.__version__,
                execution_profiles=profiles, batch_metadata_sha256=digest(batches)))
    missing = len(plan['jobs']) - len(rows)
    failed = sum(r['status'] != 'ok' for r in rows)
    empty = sum(bool(r.get('empty_sample_a') or r.get('empty_sample_b')) for r in rows)
    lines = ['# Synthetic-set screening report', '',
             f"Purpose: `{plan['config']['purpose']}`. {len(rows)}/{len(plan['jobs'])} jobs present; "
             f'{missing} missing; {failed} failed; {empty} with empty sampled sets.', '',
             'MinHash is the current FastKMV backend. No separate Legacy MinHash candidate was tested.', '',
             'MAE is against exact constructed canonical-k-mer sets. It is not ANI error. '
             'Payload is A+B query-array bytes; RSS includes inputs and generator allocations. '
             'Native hashing and N-separated carrier parsing are included in build time.', '',
             f'Execution profile and threads: {profiles}. '
             'For main_openmp_batch, pair times are elapsed times under concurrent sample-list load; '
             'RSS is the shared process high-water mark, not per-pair memory. '
             'KSSD dictionary setup is paid once per batch and recorded in batches.csv. '
             'Batch wall time includes generation, construction and repeated queries. '
             'Timing frontiers and reuse predictions are specific to this execution profile.', '',
             'Errors for failed/missing configurations are conditional on valid observations and '
             'cannot enter the frontier. Empty-sample configurations also cannot enter the frontier.', '',
             '| Algorithm | Status | Frontier regimes |', '|---|---|---|']
    lines += [f"| {d['algorithm']} | {d['status']} | {d['frontier_regimes']}/{d['regimes']} |" for d in decision]
    lines += ['', 'Frontiers use MAE, mean pair payload, median build time and median warm query time. '
              'Pruning candidates additionally require paired error bootstrap evidence, a practical '
              'resource gain by a competitor, complete coverage, and no budget-constrained advantage. '
              'They require independent genome validation before a selection guide changes.', '',
              'Query reuse predictions describe repeated comparisons of this pair with warm caches. '
              'They do not measure all-vs-all throughput or validate an index.', '']
    (output / 'report.md').write_text('\n'.join(lines))
    if not args.no_plots:
        plot(summary, output)
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
