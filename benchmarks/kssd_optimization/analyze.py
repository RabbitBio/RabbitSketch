#!/usr/bin/env python3
from validate import REPO, ROOT, RUN, OLD, NEW
from run import sha, write
import csv
import json
import shutil
import statistics


def main():
    rows = json.loads((RUN/'cli.json').read_text())
    full = [r for r in rows if r['panel'] == 'full233' and r['repeat'] == 0]
    assert len(full) == 2
    groups = {version: [r for r in full if r['version'] == version] for version in ['before','after']}
    assert all(len(group) == 1 for group in groups.values())
    assert all(r['output_digest'] == full[0]['output_digest'] for r in full)
    published = REPO/'benchmarks/results/kssd_optimization_20261111'
    published.mkdir(parents=True, exist_ok=True)
    summary = dict(files=208753, input_bytes=945020243458, threads=64, k=20, drlevel=3,
                   cutoff=.05, repetitions=1, output_digest=full[0]['output_digest'])
    for version, group in groups.items():
        summary[version] = {
            key: statistics.mean(r[key] for r in group)
            for key in ['process_seconds','prepare_seconds','similarity_seconds','peak_rss_kib']}
        summary[version]['process_range'] = [min(r['process_seconds'] for r in group), max(r['process_seconds'] for r in group)]
    before, after = summary['before'], summary['after']
    summary['end_to_end_speedup'] = before['process_seconds']/after['process_seconds']
    summary['prepare_speedup'] = before['prepare_seconds']/after['prepare_seconds']
    summary['time_reduction_percent'] = 100*(1-after['process_seconds']/before['process_seconds'])
    summary['memory_reduction_percent'] = 100*(1-after['peak_rss_kib']/before['peak_rss_kib'])
    write(published/'summary.json', summary)
    fields = ['version','repeat','threads','process_seconds','prepare_seconds','sketch_seconds',
              'legacy_sketch_local_index_seconds','index_seconds','similarity_seconds','peak_rss_kib','binary_sha256']
    with (published/'measurements.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore', lineterminator='\n'); writer.writeheader(); writer.writerows(full)
    lines = [
        '# KSSD 优化：233 全量、64 线程对照', '',
        f"优化后端到端用时 **{after['process_seconds']:.2f} 秒**，原版本为 **{before['process_seconds']:.2f} 秒**；"
        f"加速 **{summary['end_to_end_speedup']:.2f}×**，用时减少 **{summary['time_reduction_percent']:.1f}%**。", '',
        '208,753 个输入文件，945,020,243,458 字节（880.12 GiB）；有效 k=20、drlevel=3、距离阈值 0.05。'
        '两版本均使用倒排索引，并完整计算相似度、写出结果。每版本完成一次；按用户要求停止第二组测量，以下仅报告完成并验证的第一组。', '',
        '| 指标 | 原版本 | 优化后 |', '|---|---:|---:|',
        f"| 构建及索引准备 | {before['prepare_seconds']:.2f} s | {after['prepare_seconds']:.2f} s |",
        f"| 相似度计算及写盘 | {before['similarity_seconds']:.2f} s | {after['similarity_seconds']:.2f} s |",
        f"| 端到端运行时间 | {before['process_seconds']:.2f} s | {after['process_seconds']:.2f} s |",
        f"| 进程峰值 RSS | {before['peak_rss_kib']/1048576:.2f} GiB | {after['peak_rss_kib']/1048576:.2f} GiB |", '',
        '构建及索引准备统一以内部总时间减相似度阶段时间计算。原版本的 `sketch + local index` '
        '混合了输入、Sketch 和局部倒排表，不能与优化后的纯 `sketch time` 直接相除。'
        '端到端时间额外包含参数初始化及进程开销；结果校验时间不计入。', '',
        '运行顺序为优化后→原版本，两次任务串行。CPU 固定为 0–63，'
        'OpenMP 固定绑定；使用普通文件系统缓存，没有清空缓存。由于每版仅一次，'
        '结果是本机、该数据集和默认参数下的实测，不能作为重复测量的统计结论。', '',
        '## 结果一致性', '',
        f"两次均输出 **{full[0]['output_digest']['rows']:,} 行**。行数、基因组对/距离的 128 位 sum/XOR 指纹、"
        '以及包含 common、集合大小、Jaccard、距离的完整行指纹均一致；基因组对/距离指纹也与此前'
        '已验证的原生全量结果一致。无序指纹用于完整性校验，不等同于数学证明。', '',
        '本机优化构建与可移植构建各通过 12 项 CTest；独立窗口枚举 oracle 检查无效字符、'
        '分块、32/64 位状态、投影、读后继续更新及索引；CLI 对照覆盖 drlevel=0–3、'
        'FASTA/FASTQ/gzip 和距离阈值边界。最终默认使用可复用哈希集合，并通过上述全量结果对照。', '',
        '## 已实现的改动', '',
        '- 降维映射表按实际保留的维度数分配；默认用精确位图排除未选维度，再查小映射表。',
        '- CLI 与库共享 SIMD 编码、滚动规范编码和降维内核，保留原排列、采样和估计器。',
        '- 统一 builder 共享不可变参数；带长度的流式更新累积去重，在需要有序结果时整理。旧单参数更新保持公开数组可读。',
        '- gzip 仅预读文件头识别格式；普通 FASTA 直接处理带换行的序列区间，省去序列拼接副本。',
        '- 只分配当前键宽所需的去重工作区并跨文件复用。比较数组批处理后，默认保留实测更快的哈希集合。',
        '- 用连续的 `(降维值, 基因组编号)` 数组分桶、排序后生成 CSR，替代大量小 posting 数组和加锁合并。', '',
        '原 CLI 本来已经共享参数并按基因组整理，因此全量收益不能归因于 API 中的逐 contig 排序修复。', '',
        '默认入口推广后，main（含不带 --index 的列表模式）、独立 test、双文件模式和评测程序均使用共享优化构建路径。'
        '推广后的 main/test 与全量测量二进制进行逐行一致性验证；按用户要求没有再次运行完整性能测量。', '',
        '## 独立构建实验', '',
        '以下为 1 次预热、3 次测量的中位数。输入预先读入内存；RSS 包含输入，不能当作 Sketch 单独内存。', '',
        '| 用例 | 原版本 | 优化后 | 加速 |', '|---|---:|---:|---:|',
    ]
    micro = json.loads((RUN/'micro.json').read_text())
    def median(row): return statistics.median(r['seconds'] for r in row['trials'] if r['repeat'] >= 0)
    for panel, mode, threads, label in [
        ('real128','api',1,'128 基因组，流式 API，1 线程'),
        ('real128','api',8,'128 基因组，流式 API，8 线程'),
        ('one','api',1,'单基因组拆成 1000 bp records'),
        ('builder8','builder',1,'8 个独立文件 builder（含参数初始化）')]:
        values = {r['version']: median(r) for r in micro
                  if (r['panel'],r['mode'],r['threads']) == (panel,mode,threads)}
        lines.append(f"| {label} | {values['before']:.4f} s | {values['after']:.4f} s | {values['before']/values['after']:.2f}× |")
    lines += ['', '同一批 128 个基因组、相同采样结果的查询/去重候选：', '',
              '| 候选 | 单线程中位数 |', '|---|---:|']
    for mode, label in [('bitmap','位图＋小表，数组缓冲'),('map','紧凑哈希表，数组缓冲'),
                        ('direct','直接排列数组，数组缓冲'),('bitmap-set','位图＋小表，哈希集合')]:
        row = next(r for r in micro if r['mode'] == mode)
        lines.append(f'| {label} | {median(row):.4f} s |')
    lines += ['', '这些候选输出的逐基因组有序哈希数组 SHA256 一致。参数共享用例的完整查询矩阵也一致。', '',
              '## 复现与原始记录', '',
              '- 实验流程：[README](../../kssd_optimization/README.md)',
              '- 完整测量：[measurements.csv](measurements.csv)',
              '- 汇总：[summary.json](summary.json)',
              f'- 原始日志、源代码快照、命令及结果指纹：`{RUN}`', '',
              '此前线程扩展性任务按用户要求停在 47/49；本次独立测试没有覆盖、混用或恢复旧任务。', '']
    (published/'report.md').write_text('\n'.join(lines))
    for name in ['micro.json','cli.json','validation.json','full_provenance.json','micro_provenance.json',
                 'native_tests.log','portable_tests.log','limited_team_test.log',
                 'default_validation.json','default_validation.log','python_native_tests.log','stop_second_measurement.json']:
        if (RUN/name).exists(): shutil.copy2(RUN/name, published/name)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__': main()
