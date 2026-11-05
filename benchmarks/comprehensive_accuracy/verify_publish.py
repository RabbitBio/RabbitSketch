#!/usr/bin/env python3
"""Audit completed experiments and publish compact, fingerprinted evidence."""
import ast
import gzip
import re
import shutil
import subprocess
import sys

from common import *


def copy_compact(source, target):
    if source.suffix == '.csv' and source.stat().st_size > 2 * 1024 * 1024:
        target = target.with_suffix('.csv.gz')
        with source.open('rb') as src, target.open('wb') as dest:
            with gzip.GzipFile(fileobj=dest, mode='wb', mtime=0, compresslevel=6) as encoded:
                shutil.copyfileobj(src, encoded)
    else: shutil.copyfile(source, target)
    return target


def verify():
    plan = json.loads((RUN / 'plan.json').read_text())
    identity = json.loads((RUN / 'identity.json').read_text())
    if plan['config'] != CONFIG: raise ValueError('registered config differs')
    units = pd.read_csv(RUN / 'synthetic_units.tsv', sep='\t', dtype={'data_seed': 'uint64'})
    if len(units) != 107520 or units.unit_id.nunique() != 107520: raise ValueError('synthetic unit inventory differs')
    if units.groupby('regime_id').replicate.nunique().ne(128).any(): raise ValueError('independent replication incomplete')
    if units.groupby(['regime_id', 'split', 'hash_seed']).size().ne(8).any(): raise ValueError('seed/split balance differs')
    native_sources = ROOT / 'source/RabbitSketch-validated'
    for name, expected in identity['native_source_hashes'].items():
        if sha(native_sources / name) != expected: raise ValueError('frozen native source changed')
    measurements, invalid, archive_hashes = 0, 0, {}
    for k in CONFIG['k']:
        for nb in CONFIG['larger_cardinalities']:
            directory = RUN / 'synthetic' / f'k{k}_n{nb}'
            completion = json.loads((directory / 'completion.json').read_text())
            if completion['synthetic_rows'] != 259840 or completion['synthetic_units'] != 8960:
                raise ValueError('synthetic batch incomplete')
            if completion['binary_sha256'] != identity['binary_sha256']: raise ValueError('native binary identity differs')
            path = directory / 'raw.jsonl.gz'
            if sha(path) != completion['raw_sha256']: raise ValueError('native archive hash differs')
            archive_hashes[str(path)] = completion['raw_sha256']
            measurements += completion['synthetic_rows']; invalid += completion['invalid_estimates']
            analysis = json.loads((OUT / 'synthetic_batches' / directory.name / 'analysis_complete.json').read_text())
            if analysis['analyzer_sha256'] != sha(HERE / 'analyze_synthetic.py') or analysis['common_sha256'] != sha(HERE / 'common.py'):
                raise ValueError('synthetic summary is stale')
    real_counts = {}
    for job in plan['real_jobs']:
        directory = RUN / 'real' / job['job_id']
        completion = json.loads((directory / 'completion.json').read_text())
        if completion['binary_sha256'] != identity['binary_sha256']: raise ValueError('native binary identity differs')
        path = directory / 'raw.jsonl.gz'
        if sha(path) != completion['raw_sha256']: raise ValueError('real archive hash differs')
        archive_hashes[str(path)] = completion['raw_sha256']
        kind = job['phase'] + ('_exact' if job['algorithm'] == 'exact' else '_native')
        real_counts[kind] = real_counts.get(kind, 0) + completion['rows'].get('pair', 0)
    if measurements != 3118080 or len(plan['real_jobs']) != 505: raise ValueError('registered experiment not completed')
    expected_counts = dict(natural_exact=6426 * 3, natural_native=6426 * 92 * 3,
        mixed_exact=320, mixed_native=320 * 92, containment_exact=700,
        containment_native=700 * 92, weighted_native=765 * 40)
    if real_counts != expected_counts: raise ValueError('real measurement counts differ from registered panels')
    for p in (OUT / 'real_panels').glob('*/analysis_complete.json'):
        x = json.loads(p.read_text())
        if x['analyzer_sha256'] != sha(HERE / 'analyze_real.py') or x['common_sha256'] != sha(HERE / 'common.py'):
            raise ValueError('real analysis is stale')
    summary = pd.read_csv(OUT / 'synthetic_accuracy.csv')
    if summary.regime_id.nunique() != 840: raise ValueError('summary omitted regimes')
    if summary[['regime_id', 'split', 'algorithm', 'budget', 'metric']].duplicated().any(): raise ValueError('duplicate summaries')
    if summary.observations.ne(64).any(): raise ValueError('wrong independent validation count')
    decisions = pd.read_csv(OUT / 'synthetic_decisions.csv')
    if len(decisions) != 840 * 4 * 3: raise ValueError('decision grid incomplete')
    for row in decisions[decisions.supported_winner.notna()].itertuples():
        if row.candidate_set != row.supported_winner or row.eligible_count < 2:
            raise ValueError('invalid unique-winner declaration')
    profiles = pd.read_csv(OUT / 'synthetic_constraint_profiles.csv')
    if len(profiles) != 50400: raise ValueError('constraint profile grid incomplete')
    relative_zero = (profiles.scale == 'relative') & np.where(profiles.metric == 'j_error', profiles.truth_j == 0, profiles.truth_c == 0)
    if not profiles.loc[relative_zero, 'status'].eq('inapplicable_zero_truth').all(): raise ValueError('undefined relative error used')
    # Recompute reported MAE/P95 from compact independent-input rows, including NaN handling.
    for directory in sorted((OUT / 'synthetic_batches').glob('*')):
        d = pd.read_csv(directory / 'selected_rows.csv.gz', low_memory=False)
        s = pd.read_csv(directory / 'accuracy.csv')
        for metric in METRICS:
            actual = d.groupby(['regime_id', 'split', 'algorithm', 'budget'])[metric].agg(['mean', 'count'])
            recorded = s[s.metric == metric].set_index(['regime_id', 'split', 'algorithm', 'budget'])
            actual = actual.reindex(recorded.index)
            if not np.allclose(actual['mean'], recorded.mae, equal_nan=True, atol=1e-14, rtol=1e-12):
                raise ValueError('published MAE does not match independent-input rows')
            if not np.array_equal(actual['count'], recorded.valid): raise ValueError('unavailable outputs were dropped incorrectly')
            p95 = d.groupby(['regime_id', 'split', 'algorithm', 'budget'])[metric].quantile(.95).reindex(recorded.index)
            if not np.allclose(p95, recorded.p95, equal_nan=True, atol=1e-14, rtol=1e-12): raise ValueError('P95 differs')
    # No timing columns may enter the accuracy or downstream publication tables.
    for path in OUT.glob('*.csv'):
        columns = pd.read_csv(path, nrows=0).columns
        if any('seconds' in name or 'rss' in name for name in columns): raise ValueError(f'timing column in {path}')
    for path in HERE.glob('*.py'): ast.parse(path.read_text())
    tests = subprocess.run([sys.executable, '-m', 'pytest', '-q', str(HERE / 'test_analysis.py')],
                           capture_output=True, text=True, check=True)
    (OUT / 'analysis_tests.log').write_text(tests.stdout + tests.stderr)
    verified = dict(status='passed', synthetic_native_measurements=measurements,
        synthetic_independent_inputs=len(units), synthetic_invalid_jaccard_outputs=invalid,
        native_real_jobs=len(plan['real_jobs']), real_pair_counts=real_counts,
        total_native_pair_measurements=measurements + sum(v for k, v in real_counts.items() if k.endswith('_native')),
        exact_pair_measurements=sum(v for k, v in real_counts.items() if k.endswith('_exact')),
        all_raw_archive_hashes_verified=True, all_registered_jobs_complete=True,
        synthetic_statistics_recomputed_from_input_rows=True,
        synthetic_regimes=840, realized_count_conditions=832,
        calibration_inputs_per_regime=64, validation_inputs_per_regime=64,
        no_timing_columns_in_accuracy_tables=True, test_output=tests.stdout.strip(),
        archive_hashes=archive_hashes, sources_sha256={p.name: sha(p) for p in HERE.iterdir() if p.is_file()})
    atomic(OUT / 'verification.json', verified)
    return verified


def publish(verified):
    PUBLISHED.mkdir(parents=True, exist_ok=True)
    names = ['report.md', 'verification.json', 'analysis_tests.log', 'synthetic_complete.json',
             'real_complete.json', 'downstream_complete.json', 'constraint_profiles_complete.json']
    names += [p.name for p in OUT.glob('*.csv')]
    names += [p.name for p in OUT.glob('*.png')] + [p.name for p in OUT.glob('*.pdf')]
    mapping = {}
    for name in sorted(set(names)):
        dest = copy_compact(OUT / name, PUBLISHED / name); mapping[name] = dest.name
    report_path = PUBLISHED / 'report.md'
    report = report_path.read_text()
    for old, new in mapping.items(): report = report.replace('](' + old + ')', '](' + new + ')')
    report_path.write_text(report)
    for source in ['config.json', 'PROTOCOL.md', 'ANALYSIS_ADDENDUM.md']:
        shutil.copyfile(HERE / source, PUBLISHED / source)
    for source in ['identity.json', 'real_input_sha256.json']:
        shutil.copyfile(RUN / source, PUBLISHED / source)
    shutil.copyfile(RUN / 'preflight/verification.json', PUBLISHED / 'native_preflight.json')
    content = '''# Comprehensive accuracy evidence

Read the [Chinese report](report.md) and [workflow](../../comprehensive_accuracy/README.md).
The original registered protocol, separate post-registration diagnostic addendum,
configuration and source/input identities are included here. Native raw JSONL,
complete per-input tables and the frozen genome inventory remain in the external
workspace documented by the workflow. Compressed CSV files are ordinary CSV
encoded with gzip; compression does not change the tabular values.

The report distinguishes registered inferential winners, empirical error
requirements, and downstream outcomes. Synthetic intervals are simultaneous
within each regime/metric; real and downstream cluster intervals are pointwise.
This is accuracy and query-payload evidence, not a performance comparison.

Published artifacts and omissions are fingerprinted in
[publication_manifest.json](publication_manifest.json); completion and independent
recalculation checks are in [verification.json](verification.json).
'''
    (PUBLISHED / 'README.md').write_text(content)
    artifacts = {}
    for path in sorted(PUBLISHED.iterdir()):
        if path.is_file() and path.name != 'publication_manifest.json':
            item = dict(bytes=path.stat().st_size, sha256=sha(path))
            if '.csv' in path.suffixes: item['rows'] = len(pd.read_csv(path, low_memory=False))
            artifacts[path.name] = item
    atomic(PUBLISHED / 'publication_manifest.json', dict(schema_version=1,
        registered_date='2026-10-30', resumed_date='2026-10-31', complete=True,
        native_pair_measurements=verified['total_native_pair_measurements'],
        artifacts=artifacts, source_sha256=verified['sources_sha256'],
        omissions=['Raw native JSONL archives and full per-input analysis (external workspace)',
                   'Frozen genome/weighted input inventory and native executable (external workspace)'],
        timings_used_for_selection=False))
    for document in [PUBLISHED / 'report.md', PUBLISHED / 'README.md']:
        for link in re.findall(r'\]\(([^)]+)\)', document.read_text()):
            if not link.startswith(('http:', 'https:', '#')):
                if not (document.parent / link.split('#')[0]).resolve().exists():
                    raise ValueError(f'broken evidence link: {document}: {link}')
    atomic(RUN / 'complete.json', dict(status='complete', verified_native_pair_measurements=verified['total_native_pair_measurements'],
        evidence=str(PUBLISHED), verification_sha256=sha(PUBLISHED / 'verification.json')))
    print(json.dumps(dict(published=str(PUBLISHED), files=len(artifacts), bytes=sum(x['bytes'] for x in artifacts.values()),
                          verified_native_pair_measurements=verified['total_native_pair_measurements'])), flush=True)


if __name__ == '__main__': publish(verify())
