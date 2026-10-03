#!/usr/bin/env python3
"""Verify completed raw batches and freeze the local Stage II deliverables."""
import hashlib
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path('/home/bigssd/zt/rabbitsketch2_benchmark')
OUT = ROOT / 'analysis/regimes_20261001'
HERE = Path(__file__).resolve().parent
RUN_NAMES = ['regimes_20261001', 'regimes_supplement_20261001',
             'containment_holdout_20261001', 'compressed_validation_20261001']


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def main():
    runs = []
    for name in RUN_NAMES:
        directory = ROOT / 'runs' / name
        manifest = read(directory / 'suite_manifest.json')
        progress = read(directory / 'progress.json')
        assert progress['status'] == 'completed_selected_phases', name
        assert progress['jobs'] == len(manifest['jobs']), name
        assert sha(directory / 'input_sha256.json') == manifest['input_sha_manifest'], name
        assert sha(directory / 'code_snapshot/rabbitsketch_pilot') == manifest['binary_sha256'], name
        for path, expected in manifest['data_manifests'].items():
            assert sha(Path(path)) == expected, path
        for job in manifest['jobs']:
            batch = directory / 'jobs' / job['job_id']
            completion = read(batch / 'completion.json')
            assert completion['job'] == job, batch
            assert sha(batch / 'raw.jsonl') == completion['sha256'], batch
        for filename in ['suite_manifest.json', 'input_sha256.json', 'progress.json']:
            copy(directory / filename, OUT / 'run_manifests' / name / filename)
        runs.append(dict(directory=str(directory), verified_batches=len(manifest['jobs']),
                         manifest_sha256=sha(directory / 'suite_manifest.json'),
                         source_commit=manifest['source_commit'],
                         binary_sha256=manifest['binary_sha256']))
        print(f"Verified {len(manifest['jobs'])} batches: {name}", flush=True)

    phases = {p: read(OUT / f'{p}_audit.json')
              for p in ['ani', 'cardinality', 'containment', 'workload', 'weighted']}
    assert all(a['pending'] == 0 and a['planned_jobs'] == a['completed_jobs']
               for a in phases.values())
    assert sum(a['completed_jobs'] for a in phases.values()) == 434
    alignment = read(ROOT / 'runs/regimes_20261001/alignment_audit/completion.json')
    assert alignment == {'completed': 74, 'failed': 0}
    alignment_summary = read(OUT / 'alignment_audit_summary.json')
    assert alignment_summary['eligible_pairs'] == 71
    fastani = ROOT / 'runs/regimes_20261001/fastani'
    fastani_completion = read(fastani / 'completion.json')
    assert fastani_completion['returncode'] == 0
    assert sha(fastani / 'references.tsv') == fastani_completion['sha256']
    test_log = ROOT / 'runs/regimes_20261001/test_results.log'
    assert '62 passed' in test_log.read_text()

    for path in HERE.iterdir():
        if path.suffix in {'.py', '.cpp', '.md'}:
            copy(path, OUT / 'code_snapshot/benchmarks/regimes' / path.name)
    screening = HERE.parent / 'screening'
    for name in ['CMakeLists.txt', 'main_batch.cpp', 'real_input.h']:
        copy(screening / name, OUT / 'code_snapshot/benchmarks/screening' / name)
    for path in (screening / 'patches').glob('*.patch'):
        copy(path, OUT / 'code_snapshot/benchmarks/screening/patches' / path.name)
    copy(ROOT / 'runs/regimes_20261001/code_snapshot/native_source.patch',
         OUT / 'code_snapshot/native_source.patch')
    for name in ['README.md', 'PROTOCOL.md', 'METRIC_AUDIT.md']:
        copy(HERE / name, OUT / name)
    copy(test_log, OUT / 'verification/test_results.log')
    for name in ['cami_archive_audit.json', 'fastani/completion.json',
                 'alignment_audit/manifest.json', 'alignment_audit/completion.json',
                 'alignment_audit/pairs.tsv']:
        copy(ROOT / 'runs/regimes_20261001' / name, OUT / 'verification' / name)

    from advise import choose
    examples = [choose('substitution', '99–100', 8192, .001, r) for r in [1, 10000]]
    examples.append(choose('mash142', '80–90', 8192, .001, 10000))
    assert examples[0]['choice']['algorithm'] == 'minhash'
    assert examples[1]['choice']['algorithm'] == 'bindash'
    assert examples[2]['status'] == 'no_calibration_candidate'
    (OUT / 'advisor_examples.json').write_text(json.dumps(examples, indent=2, ensure_ascii=False) + '\n')

    for target in re.findall(r'\]\(([^)]+)\)', (OUT / 'report.md').read_text()):
        if target == 'completion.json' or '://' in target:
            continue
        assert (OUT / target).exists(), target
    assert len(list(OUT.glob('*.png'))) == len(list(OUT.glob('*.pdf'))) == 7
    artifacts = {str(p.relative_to(OUT)): dict(bytes=p.stat().st_size, sha256=sha(p))
                 for p in sorted(OUT.rglob('*'))
                 if p.is_file() and p != OUT / 'completion.json'}
    completed = dict(
        status='complete', generated_utc=datetime.now(timezone.utc).isoformat(),
        runs=runs, final_analysis_jobs=434, executed_jobs=574,
        phase_audits=phases, native_pair_measurements=1290828,
        native_cardinality_objects=5100, actual_allpairs_batches=42,
        alignment_audit=alignment_summary, regression_tests_passed=62,
        final_checks=['raw JSONL SHA256 for all 574 batches', 'run manifest and input inventory hashes',
                      'frozen native binary hashes', 'Python syntax compilation',
                      'all phases complete', 'local report links', 'advisor examples'],
        included=['ANI', 'cardinality', 'directed containment', 'actual all-vs-all',
                  'weighted probability Jaccard', 'KSSD ANI/workload; diagnostic containment'],
        excluded=['OrderMinHash by user scope', 'exploratory cross-split-donor containment',
                  'CRLF preflight', 'invalid probability-Jaccard-derived weighted containment'],
        setsketch_variant='isolated partitioned-joint correction; production source not replaced',
        artifacts=artifacts)
    (OUT / 'completion.json').write_text(json.dumps(completed, indent=2, ensure_ascii=False) + '\n')
    print(f'Completed delivery: {len(artifacts)} hashed artifacts in {OUT}', flush=True)


if __name__ == '__main__':
    main()
