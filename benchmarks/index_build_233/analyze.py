#!/usr/bin/env python3
"""Audit finished measurements and report both construction and query costs."""
import json
import shutil
from pathlib import Path
from run import RUN, PUBLISHED, REPO, sha


def main():
    state = json.loads((RUN / 'progress.json').read_text())
    assert state['phase'] == 'complete' and state['completed_runs'] == 4
    plan = json.loads((RUN / 'plan.json').read_text())
    for name, expected in plan['source_sha256'].items():
        assert sha(REPO / name) == expected, 'source changed: ' + name
    for name, expected in plan['binary_sha256'].items():
        assert sha(name) == expected, 'binary changed: ' + name
    results = json.loads((RUN / 'results.json').read_text())
    assert len(results) == 4
    lines = ['# 233 全量真实 Sketch 的倒排索引构建对照', '',
        '**分桶排序显著降低了索引构建耗时和内存，但本轮没有在两个算法上都改善“索引＋相似度及输出”总时间。**', '',
        '输入为完整 233 数据集：208,753 个文件、945,020,243,458 字节（约 880.119 GiB）。',
        '64 线程、任务串行；FastKMV K=1024，ProbMinHash registers=1024，k=21、seed=42，严格距离阈值 0.05。', '',
        '每个算法调用实际 main 的构建函数生成一次全量真实键快照。旧哈希表路径和新分桶路径分别在独立进程中读取同一快照，',
        '执行完整索引构建、原生相似度计算和 TSV 输出。每个版本测量一次，顺序为旧→新，没有方差或置信区间。', '',
        '| 算法 | 方案 | 索引构建/秒 | 相似度及输出/秒 | 两阶段合计/秒 | 索引结束前峰值/GiB |',
        '|---|---|---:|---:|---:|---:|']
    names = {'fastkmv': 'MinHash（FastKMV）', 'probminhash': 'ProbMinHash'}
    modes = {'map': '线程局部表＋分片归并', 'bucket': '分桶排序＋CSR'}
    for row in results:
        lines.append(f"| {names[row['algorithm']]} | {modes[row['mode']]} | {row['index_seconds']:.3f} | {row['query_seconds']:.3f} | {row['index_query_seconds']:.3f} | {row['index_peak_rss_kib']/2**20:.3f} |")
    summary = {}
    for algorithm in names:
        pair = {r['mode']: r for r in results if r['algorithm'] == algorithm}
        old, new = pair['map'], pair['bucket']
        assert old['output_digest'] == new['output_digest']
        assert old['index_keys'] == new['index_keys'] and old['postings'] == new['postings']
        assert old['matches_historical_full_output'] and new['matches_historical_full_output']
        preparation = json.loads((RUN / 'full' / (algorithm + '_prepare.json')).read_text())
        assert sha(RUN / (algorithm + '.keys')) == preparation['snapshot_sha256']
        speed = old['index_seconds'] / new['index_seconds']
        memory = 100 * (1 - new['index_peak_rss_kib'] / old['index_peak_rss_kib'])
        query = 100 * (new['query_seconds'] / old['query_seconds'] - 1)
        total = 100 * (new['index_query_seconds'] / old['index_query_seconds'] - 1)
        summary[algorithm] = dict(index_speedup=speed, index_peak_memory_reduction_percent=memory,
            query_time_change_percent=query, index_query_time_change_percent=total,
            original_sketch_seconds=preparation['sketch_seconds'],
            output_rows=new['output_digest']['rows'], matches_previous_full_output=True)
        lines += ['', f"**{names[algorithm]}**：索引构建加速 **{speed:.2f}×**，索引结束前峰值内存降低 **{memory:.2f}%**。",
                  f"相似度及输出阶段耗时增加 **{query:.2f}%**；两阶段合计{'增加' if total > 0 else '减少'} **{abs(total):.2f}%**。",
                  f"输出 **{new['output_digest']['rows']:,}** 行，新旧指纹及历史全量参考一致。"]
    lines += ['', '## 计时和内存口径', '',
        '- 旧索引计时包含局部哈希表插入、64 分片归并、单例删除与 CSR 展平；没有遗漏插入成本。',
        '- 新索引计时包含分桶计数、前缀和、散射、桶内排序、单例删除与 CSR 构建。',
        '- 快照读取、FastKMV 原生对象恢复及结果校验不计入上述两个阶段。全进程 wall/RSS 另存于 results.json。',
        '- 表中内存为到索引构建结束为止的进程峰值，包含输入键；FastKMV 还包含恢复后的原生 Sketch。',
        '- 原始数据读取和 Sketch 构建独立完成一次：FastKMV '
        f"{summary['fastkmv']['original_sketch_seconds']:.3f} 秒，ProbMinHash {summary['probminhash']['original_sketch_seconds']:.3f} 秒。",
        '- 因此本表不是两次从原始 FASTA/FASTQ 开始的 CLI 端到端测量，也不把公共准备时间相加后冒充端到端结果。',
        '- 使用正常文件缓存，不清空全机缓存。四项计时及校验任务都串行执行。', '',
        '## 结论范围', '',
        '这轮真实数据支持“分桶排序更快地构建索引，并降低索引构建阶段峰值内存”。',
        '两个算法的相似度及输出阶段都观察到变慢；FastKMV 的构建收益未抵消这部分增加，ProbMinHash 则有净收益。',
        '当前不能宣称统一后的流程全面更快，也不能把查询变慢直接归因于某一种缓存或列表顺序效应；这需要单独验证。', '',
        '## 校验与复现', '',
        '128 样本面板逐行对照当前生产 CLI；并行输出校验器与已有串行校验器在真实面板和分块边界用例中一致。',
        '全量比较校验索引键数、posting 数、输出行数、字节数、规范化样本对的 XXH3-128 累加及异或指纹，',
        '并与既有全量实验参考核对。指纹是完整性检查，不是数学上的无碰撞证明。',
        '结束后再次确认源文件、二进制和键快照哈希未改变。', '',
        f'原始日志、命令、键快照、完整指纹及输出样本：`{RUN}`。',
        '大体积 TSV 在完成新旧与历史参考校验后删除；键快照保留。']
    (PUBLISHED / 'report.md').write_text('\n'.join(lines) + '\n')
    (PUBLISHED / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    for name in ('plan.json', 'validation.json'):
        shutil.copyfile(RUN / name, PUBLISHED / name)
    for algorithm in names:
        name = algorithm + '_prepare.json'
        shutil.copyfile(RUN / 'full' / name, PUBLISHED / name)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
