#!/usr/bin/env python3
"""Sequential full-233 A/B test on identical native sketch-key snapshots."""
import csv
import datetime as dt
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
ROOT = Path('/home/bigssd/zt/rabbitsketch2_benchmark')
TAG = 'index_build_233_20261112'
RUN = ROOT / 'runs' / TAG
BUILD = ROOT / ('build_' + TAG)
PREVIOUS = ROOT / 'runs/thread_scaling_index_233_20261109'
OUTPUT = REPO.parent / '.rabbitsketch_benchmark_outputs' / TAG
PUBLISHED = REPO / 'benchmarks/results' / TAG
CHILD = None


def write(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def environment(threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_THREAD_LIMIT=str(threads),
               OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close', OMP_PLACES='cores',
               OMP_MAX_ACTIVE_LEVELS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', LC_ALL='C')
    for name in ('RABBIT_FORCE_ALLPAIRS', 'RABBIT_FORCE_INDEX', 'RABBITSKETCH_SIMD'):
        env.pop(name, None)
    return env


def status(phase, **values):
    state = dict(phase=phase, updated_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                 runner_pid=os.getpid(), **values)
    write(RUN / 'progress.json', state)
    print(phase, values, flush=True)


def execute(label, command, threads, directory):
    global CHILD
    directory.mkdir(parents=True, exist_ok=True)
    stdout, stderr = directory / (label + '.stdout'), directory / (label + '.stderr')
    env = environment(threads)
    write(directory / (label + '.command.json'), dict(argv=list(map(str, command)),
          cwd=str(REPO / 'examples'), threads=threads,
          environment={k: v for k, v in env.items() if k.startswith(('OMP_', 'OPENBLAS_', 'MKL_')) or k == 'LC_ALL'}))
    started = time.monotonic()
    with stdout.open('w') as out, stderr.open('w') as err:
        CHILD = subprocess.Popen(list(map(str, command)), cwd=REPO / 'examples', env=env,
                                 stdout=out, stderr=err, start_new_session=True)
        while True:
            try:
                code = CHILD.wait(timeout=10)
                break
            except subprocess.TimeoutExpired:
                log = stderr.read_text(errors='replace')
                status('running', job=label, child_pid=CHILD.pid,
                       elapsed_seconds=time.monotonic() - started,
                       log=str(stderr), last_lines=log.splitlines()[-3:])
    CHILD = None
    if code:
        raise RuntimeError(f'{label}: exit {code}; see {stderr}')
    wall = time.monotonic() - started
    status('job_complete', job=label, wall_seconds=wall)
    return wall, stdout


def driver(mode, algorithm, listing, snapshot, threads, output):
    return ['taskset', '--cpu-list', ','.join(map(str, range(threads))), BUILD / 'index_build_233',
            mode, algorithm, listing, snapshot, str(threads), output]


def fingerprint(label, output, directory):
    _, stdout = execute(label, [BUILD / 'output_digest', output, '16'], 16, directory)
    result = json.loads(stdout.read_text())
    write(directory / (label + '.json'), result)
    return result


def validate():
    directory = RUN / 'validation'
    directory.mkdir(exist_ok=True)
    listing = PREVIOUS / 'smoke.list'
    cases = []
    for algorithm in ('fastkmv', 'probminhash'):
        snapshot = directory / (algorithm + '.keys')
        execute(algorithm + '_prepare', driver('prepare', algorithm, listing, snapshot, 8, '-'), 8, directory)
        reference = directory / (algorithm + '_cli.tsv')
        execute(algorithm + '_cli', [REPO / 'examples/rabbitsketch', '--' + algorithm,
                '--index', '-l', '-i', listing, '-t', '8', '-o', reference], 8, directory)
        expected = sorted(reference.read_text().splitlines())
        serial = ROOT / 'build_thread_scaling_index_233_20261109/output_digest'
        serial_digest = json.loads(subprocess.check_output([serial, reference], text=True))
        if fingerprint(algorithm + '_cli_digest', reference, directory) != serial_digest:
            raise ValueError('parallel fingerprint differs from existing serial checker')
        for mode in ('map', 'bucket'):
            output = directory / f'{algorithm}_{mode}.tsv'
            _, data = execute(algorithm + '_' + mode,
                driver(mode, algorithm, listing, snapshot, 8, output), 8, directory)
            if sorted(output.read_text().splitlines()) != expected:
                raise ValueError(f'{algorithm}/{mode} differs from actual CLI on 128 inputs')
            digest = fingerprint(algorithm + '_' + mode + '_digest', output, directory)
            if digest != serial_digest:
                raise ValueError('replay fingerprint mismatch')
            cases.append(dict(algorithm=algorithm, mode=mode, rows=len(expected),
                              statistics=json.loads(data.read_text()), matches_current_cli=True))
    # Exercise chunk boundaries, reversed pairs and variable-length paths.
    edge = directory / 'digest_edges.tsv'
    edge.write_text(''.join(f'z{i}\ta{i}{"x" * (i % 97)}\t0.012345\n' for i in range(10007)))
    serial_digest = json.loads(subprocess.check_output([serial, edge], text=True))
    for threads in (1, 3, 16):
        digest = json.loads(subprocess.check_output([BUILD / 'output_digest', edge, str(threads)], text=True))
        if digest != serial_digest:
            raise ValueError('parallel checker failed chunk-boundary validation')
    write(RUN / 'validation.json', dict(all_passed=True, cases=cases, digest_matches_serial=True))


def register():
    listing = REPO / 'examples/233.list'
    expected = json.loads((PREVIOUS / 'dataset.json').read_text())
    if sha(listing) != expected['original_list_sha256']:
        raise ValueError('233 list differs from the preserved dataset')
    count = total = 0
    with (PREVIOUS / 'input_metadata.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            info = Path(row['path']).stat()
            if (info.st_size, info.st_mtime_ns, info.st_ino, info.st_dev) != tuple(
                    int(row[x]) for x in ('bytes', 'mtime_ns', 'inode', 'device')):
                raise ValueError('input metadata changed: ' + row['path'])
            count += 1
            total += info.st_size
            if count % 32768 == 0:
                status('checking_inputs', checked=count)
    if count != expected['files'] or total != expected['file_bytes']:
        raise ValueError('dataset count or bytes changed')
    shutil.copyfile(listing, RUN / '233.list')
    sources = [p for p in (REPO / 'src').rglob('*') if p.is_file()]
    sources += [REPO / 'examples/sub_command.cpp', REPO / 'examples/sub_command.h', REPO / 'CMakeLists.txt']
    sources += [p for p in HERE.iterdir() if p.is_file() and p.suffix in ('.py', '.cpp', '.txt', '.md')]
    plan = dict(files=count, input_bytes=total, input_list_sha256=sha(listing), threads=64,
        algorithms=['fastkmv', 'probminhash'], modes=['map', 'bucket'], repetitions=1,
        defaults=dict(K=1024, registers=1024, kmer=21, seed=42, distance_cutoff=.05),
        source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
        source_sha256={str(p.relative_to(REPO)): sha(p) for p in sources},
        binary_sha256={str(p): sha(p) for p in (BUILD / 'index_build_233', BUILD / 'output_digest', REPO / 'examples/rabbitsketch')},
        machine=subprocess.check_output(['lscpu'], text=True),
        compiler=subprocess.check_output(['/usr/bin/g++', '--version'], text=True),
        cache_policy='normal filesystem cache; no cache eviction; snapshots prepared once per algorithm',
        timing='index timer covers local-map insertion + shard merge + CSR, or count/scatter/sort + CSR; query/output timed separately; snapshot IO and native restoration excluded from these two timers',
        scope='full real 233 sketch keys, with native similarities and output; not repeated end-to-end FASTX CLI measurements')
    write(RUN / 'plan.json', plan)


def publish(results):
    PUBLISHED.mkdir(parents=True, exist_ok=True)
    write(PUBLISHED / 'results.json', results)
    lines = ['# 233 全量真实 Sketch 的倒排索引构建对照', '',
        '208,753 个输入文件，945,020,243,458 字节（约 880.119 GiB）。64 线程，各任务串行运行。',
        'FastKMV 的 K=1024，ProbMinHash 的寄存器数=1024；k=21、seed=42、距离阈值=0.05。', '',
        '每个算法调用实际 main 的构建函数生成一次全量键快照。旧、新构建器在独立进程中读取同一份快照，',
        '分别构建索引、执行原生相似度计算并输出距离。每种构建器测量一次，不给出方差或置信区间。',
        '索引计时包含旧路径的局部表插入和分片归并。快照读取、FastKMV 对象恢复及结果校验不计入索引计时。',
        '这不是从 FASTA/FASTQ 输入开始的两次完整 CLI 端到端对照；原始 Sketch 准备时间单独保存。', '',
        '| 算法 | 构建方式 | 完整索引构建/秒 | 相似度及输出/秒 | 两阶段合计/秒 | 索引结束前峰值/GiB | 输出行数 |',
        '|---|---|---:|---:|---:|---:|---:|']
    for row in results:
        lines.append(f"| {row['algorithm']} | {row['mode']} | {row['index_seconds']:.6f} | {row['query_seconds']:.6f} | {row['index_query_seconds']:.6f} | {row['index_peak_rss_kib']/2**20:.3f} | {row['output_digest']['rows']:,} |")
    lines += ['', '内存是到索引构建结束为止的进程峰值，包含输入键；FastKMV 还包含原生 Sketch 对象。',
              '输出校验使用与既有实验一致的行数、字节数、规范化样本对的 XXH3-128 累加及异或指纹。',
              '指纹是完整性检查，不是数学上的无碰撞证明。', '',
              f'原始日志、键快照、命令与源代码身份：`{RUN}`。']
    for algorithm in ('fastkmv', 'probminhash'):
        pair = {r['mode']: r for r in results if r['algorithm'] == algorithm}
        if len(pair) == 2:
            old, new = pair['map'], pair['bucket']
            speed = old['index_seconds'] / new['index_seconds']
            reduction = 1 - new['index_peak_rss_kib'] / old['index_peak_rss_kib']
            lines.append(f"\n{algorithm}：索引构建加速比 {speed:.3f}×；索引结束前峰值内存降低 {reduction*100:.2f}%；新旧输出指纹一致。")
    (PUBLISHED / 'report.md').write_text('\n'.join(lines) + '\n')


def measure():
    results = []
    listing = RUN / '233.list'
    directory = RUN / 'full'
    directory.mkdir(exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for algorithm in ('fastkmv', 'probminhash'):
        snapshot = RUN / (algorithm + '.keys')
        label = algorithm + '_prepare'
        command = ['/usr/bin/time', '-v', '-o', directory / (label + '.time'),
                   *driver('prepare', algorithm, listing, snapshot, 64, '-')]
        wall, stdout = execute(label, command, 64, directory)
        preparation = json.loads(stdout.read_text())
        preparation.update(process_seconds=wall, snapshot_sha256=sha(snapshot), snapshot_bytes=snapshot.stat().st_size)
        write(directory / (label + '.json'), preparation)
        expected = json.loads((PREVIOUS / 'jobs' / (algorithm + '_t64_r1') / 'output_digest.json').read_text())
        pair = []
        for mode in ('map', 'bucket'):
            label = algorithm + '_' + mode
            output = OUTPUT / (label + '.tsv')
            if output.exists():
                raise RuntimeError('refusing to overwrite ' + str(output))
            command = ['/usr/bin/time', '-v', '-o', directory / (label + '.time'),
                       *driver(mode, algorithm, listing, snapshot, 64, output)]
            wall, stdout = execute(label, command, 64, directory)
            row = json.loads(stdout.read_text())
            timing = (directory / (label + '.time')).read_text()
            row.update(algorithm=algorithm, mode=mode, process_seconds=wall,
                       process_peak_rss_kib=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', timing)[1]))
            write(directory / (label + '.timing.json'), row)
            digest = fingerprint(label + '_digest', output, directory)
            row.update(output_digest=digest, matches_historical_full_output=digest == expected)
            if pair and (digest != pair[0]['output_digest'] or
                         row['postings'] != pair[0]['postings'] or row['index_keys'] != pair[0]['index_keys']):
                raise ValueError('map/bucket index or full output mismatch; retaining output files')
            pair.append(row)
            results.append(row)
            write(RUN / 'results.json', results)
            publish(results)
            with output.open() as stream:
                (directory / (label + '.sample.tsv')).write_text(''.join(itertools.islice(stream, 101)))
        if not all(row['matches_historical_full_output'] for row in pair):
            raise ValueError('new/old agree but historical full output differs; retain files and investigate')
        # Only remove this run's large, fully verified output; keep snapshots,
        # samples, complete fingerprints and all measured timings.
        for mode in ('map', 'bucket'):
            (OUTPUT / f'{algorithm}_{mode}.tsv').unlink()
        status('algorithm_complete', algorithm=algorithm)
    status('complete', completed_runs=len(results), report=str(PUBLISHED / 'report.md'))


def interrupted(signum, frame):
    if CHILD is not None and CHILD.poll() is None:
        os.killpg(CHILD.pid, signal.SIGTERM)
        try:
            CHILD.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(CHILD.pid, signal.SIGKILL)
            CHILD.wait()
    raise InterruptedError(f'signal {signum}')


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    if (RUN / 'plan.json').exists():
        raise RuntimeError('registered experiment already exists; refusing to repeat measurements')
    lock = (ROOT / 'runs/.active_benchmark.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, interrupted)
    status('validating')
    validate()
    register()
    measure()


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        status('failed', error=str(error))
        raise
