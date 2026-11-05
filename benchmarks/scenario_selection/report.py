#!/usr/bin/env python3
"""Publish auditable tables, frozen choices/model, figures and a factual report."""
from concurrent.futures import ThreadPoolExecutor
import re
import shutil
import subprocess
import sys

from advisor import PILOT_FEATURES, read_scores
from scenario_common import *


def fmt(value):
    if value is None or (isinstance(value,(float,np.floating)) and not np.isfinite(value)):return '—'
    if isinstance(value,(float,np.floating)):return f'{value:.5f}'
    return str(value).replace('|','/')


def table(frame,columns):
    names=list(columns);headers=[columns[n] for n in names]
    lines=['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(names))+' |']
    for values in frame[names].itertuples(index=False,name=None):lines.append('| '+' | '.join(fmt(v) for v in values)+' |')
    return '\n'.join(lines)


def publication_readme():
    return '''# Scenario selection: independent genome validation

Start with the [Chinese selection guide](../../../RabbitSketch2_scenario_selection.md)
or the [complete report](report.md). The [candidate map](figures/scenario_candidates_k24_8k.png)
explicitly marks cells where no measured configuration meets the standard target.

The experiment uses 450 new genomes from 30 genera held out from the 60-genus
development panel, 14 validation conditions, three k values and four query-payload
budgets. These genera occurred in older broad panels; see
[the scope correction](REPORTING_NOTES.md). SetSketch results use the isolated
corrected implementation. No runtime ranking is made.

The frozen rule does not establish a universal selector or reliable downstream
gains. Relative error reduction and practical target attainment are reported
separately. The full tables, model, frozen decisions, truth/pilot data, source and
input identities, and [verification results](publication_verification.json) are
included here. [Workflow and reproduction](../../scenario_selection/README.md).
'''


def input_integrity():
    manifest=RUN/'input_sha256.json';done=OUT/'inputs_verified.json'
    if done.exists():
        prior=json.loads(done.read_text())
        if prior['manifest_sha256']!=sha(manifest):raise ValueError('verified input manifest changed')
        return prior
    inputs=json.loads(manifest.read_text());paths=sorted(inputs)
    with ThreadPoolExecutor(max_workers=16) as pool:
        for i,(p,digest) in enumerate(zip(paths,pool.map(sha,paths)),1):
            if digest!=inputs[p]:raise ValueError('input content changed: '+p)
            if i%1000==0:print('Verified input SHA256',i,'/',len(paths),flush=True)
    result=dict(manifest_sha256=sha(manifest),files=len(paths),all_match=True,
        total_input_bytes=sum(Path(p).stat().st_size for p in paths))
    atomic(done,result);return result


def input_condition_summary():
    definitions=CONFIG['conditions']+json.loads((HERE/'validation_extra_conditions.json').read_text())
    requested={c['name']:c['length'] for c in definitions if c['method']=='fragment'}
    rows=[]
    for split in ['development','validation']:
        d=pd.read_csv(DATA/(split+'_objects.tsv'),sep='\t');d=d[d.role=='query']
        for (condition,kind),g in d.groupby(['condition','kind']):
            rows.append(dict(split=split,condition=condition,kind=kind,queries=len(g),
                source_genomes=g.source.nunique(),genera=g.genus.nunique(),
                min_bases=int(g.bases.min()),median_bases=float(g.bases.median()),max_bases=int(g.bases.max()),
                median_records=float(g.records.median()),median_record_length=float(g.median_record_length.median()),
                requested_fragment_length=requested.get(condition),
                truncated_fragments=int((g.bases<requested[condition]).sum()) if condition in requested else 0))
    save_csv(OUT/'input_condition_summary.csv',pd.DataFrame(rows))


def shared_exact_check():
    rows=[]
    for k in CONFIG['k']:
        old=ROOT/'analysis/comprehensive_accuracy_20261030/real_panels'/f'natural_k{k}'
        a=pd.read_csv(old/'exact_truth.csv')
        b=pd.read_csv(OUT/'development'/f'exact_k{k}_pairs.csv.gz')
        b=b.merge(pd.read_csv(DATA/'development_pairs.tsv',sep='\t'),on='pair_id',validate='one_to_one')
        lookup={tuple(sorted([r.a,r.b])):int(r.intersection) for r in a.itertuples()}
        matched=set();observations=0
        for r in b.itertuples():
            key=tuple(sorted([r.a,r.b]))
            if key not in lookup:continue
            if int(r.intersection)!=lookup[key]:raise ValueError('shared exact intersection changed')
            matched.add(key);observations+=1
        if len(matched)!=6300:raise ValueError('missing original within-genus full-genome pairs')
        ao=pd.read_csv(old/'objects.csv.gz',usecols=['algorithm','object_id','n'])
        ao=ao[ao.algorithm=='exact'].set_index('object_id').n
        bo=pd.read_csv(OUT/'development'/f'exact_k{k}_objects.csv.gz').set_index('object_id').n
        if not np.array_equal(ao.to_numpy(),bo.loc[ao.index].to_numpy()):raise ValueError('shared exact cardinality changed')
        rows.append(dict(k=k,shared_genomes=len(ao),shared_unordered_pairs=len(matched),
            ordered_pair_checks=observations,all_counts_match=True,prior_exact_sha256=sha(old/'exact_truth.csv')))
    save_csv(OUT/'shared_exact_verification.csv',pd.DataFrame(rows))


def native_inventory():
    rows=[]
    for plan_name in ['development','validation-pilot','validation-outcomes']:
        plan=json.loads((RUN/(plan_name+'_plan.json')).read_text())
        for job in plan['jobs']:
            directory=RUN/'native'/job['job_id'];completion=json.loads((directory/'completion.json').read_text())
            rows.append(dict(job_id=job['job_id'],split=job['split'],purpose=job['purpose'],
                algorithm=job['algorithm'],parameter=job['parameter'],half_subk=job['half_subk'],k=job['k'],seed=job['seed'],
                pairs=completion['rows']['pair'],objects=completion['rows']['object'],
                raw_sha256=completion['raw_sha256'],manifest_sha256=completion['manifest_sha256'],
                binary_sha256=completion['binary_sha256']))
    d=pd.DataFrame(rows)
    if d.job_id.duplicated().any():raise ValueError('duplicate published native job')
    save_csv(OUT/'native_archive_manifest.csv',d)
    return d


def report_text(inventory):
    audit=json.loads((DATA/'independence_audit.json').read_text())
    dev=json.loads((DATA/'development_preparation.json').read_text())
    val=json.loads((DATA/'validation_preparation.json').read_text())
    evaluation=json.loads((OUT/'evaluation_complete.json').read_text())
    summary=pd.read_csv(OUT/'policy_summary.csv')
    intervals=pd.read_csv(OUT/'advisor_vs_fixed_intervals.csv')
    profiles=pd.read_csv(OUT/'precision_profile_summary.csv')
    candidates=pd.read_csv(OUT/'scenario_candidates.csv').fillna('')
    downstream=pd.read_csv(OUT/'downstream_summary.csv')
    primary=intervals[(intervals.policy=='pilot')&(intervals.condition=='all')]
    worsened=primary[(primary.simultaneous_low>0)&(primary.se>1e-12)]
    focus=summary[(summary.k==24)&(summary.budget==8192)&(summary.condition=='all')&summary.policy.isin(['fixed','pilot','metadata'])]
    main=primary[(primary.k==24)&(primary.budget==8192)].copy()
    fixed=pd.read_csv(OUT/'development_fixed_choices.csv')
    fixed=fixed[(fixed.k==24)&(fixed.budget==8192)&(fixed.policy=='fixed')]
    candidate_focus=candidates[(candidates.k==24)&(candidates.budget==8192)].copy()
    candidate_focus['display']=candidate_focus.apply(lambda r:
        ('全部已测配置未达标；最低损失：'+NAMES.get(r.best_observed,r.best_observed))
        if r.feasible_alternative_rate==0 else
        (' / '.join(NAMES.get(a,a) for a in str(r.candidates).split(';') if a) or '无可用候选')+
        (' †' if r.feasible_alternative_rate<1 else ''),axis=1)
    cand=candidate_focus.pivot(index='condition',columns='metric',values='display').reset_index()
    c_comparison=main[main.metric=='containment'].iloc[0]
    cp=profiles[(profiles.k==24)&(profiles.budget==8192)&(profiles.condition=='all')&
                (profiles.metric=='containment')&(profiles.profile=='standard')].set_index('policy')
    chosen_profiles=profiles[(profiles.k==24)&(profiles.budget==8192)&(profiles.condition=='all')&
                             (profiles.policy=='pilot')&(profiles.profile=='standard')]
    transfer=intervals[(intervals.policy=='pilot')&(intervals.k==24)&(intervals.budget==8192)&
                       intervals.condition.isin(['original_conditions','unseen_conditions'])]
    open_set=downstream[(downstream.metric=='species')&(downstream.k==24)&(downstream.budget==8192)&
                        downstream.policy.isin(['fixed','pilot'])&downstream.condition.isin(['full','fragment_1k','reads_cov1_error03','unseen_reads_cov05_error02'])]
    zero=summary[(summary.metric=='jaccard')&(summary.k==24)&(summary.budget==8192)&(summary.policy=='pilot')&
                  summary.condition.isin(['full','fragment_1k','fragment_10k','unseen_fragment_5k'])]
    counts=inventory[(inventory.purpose=='outcomes')&(inventory.algorithm!='exact')].pairs.sum()
    exact_pairs=inventory[inventory.algorithm=='exact'].pairs.sum()
    lines=[
        '# RabbitSketch2：可观测场景选型的独立验证（2026-11-05）',
        '',
        '本轮完成了“可观测输入 → 冻结规则和参数 → 按属留出的新基因组、新随机种子、未见输入条件验证”的实验链条。'
        '目标是检验精度与空间约束下的选择，不预设自动选择一定有效。下面的候选图与冻结规则是两个不同产物：'
        '候选图总结验证后的证据，规则的全部决定则在验证结果计算前已经保存。',
        '',
        f'冻结的 pilot 规则在 72 个预定义任务/k/预算汇总对比中，有 **{int(primary.material_improvement.sum())} 项**同时满足'
        '跨 72 项的 95% bootstrap 区间和至少 10%（且绝对 0.0001）的损失下降门槛；'
        f'另有 **{len(worsened)} 项**的同时区间明确落在固定方案更好的方向。'
        '这些计数不能替代各任务的效应大小，也不能解释成通用最优选择器已经成立。',
        '',
        f'尤其是 k=24/8 KiB：containment 平均损失从 {c_comparison.fixed_loss:.6f} 降到 {c_comparison.policy_loss:.6f}，'
        f'但实际达标率仍为固定方案 {cp.loc["fixed","attainment_rate"]:.2%}、规则 {cp.loc["pilot","attainment_rate"]:.2%}。'
        '该汇总收益全部来自 1 kb 片段；这个场景下所有已测配置仍未达到预设 containment 目标。'
        '因此这项统计收益不能作为“选择产生了可用方案或改善下游”的证据。',
        '',
        '## 数据与独立性',
        '',
        f'- 开发：{dev["genomes"]} 个旧基因组、{dev["genera"]} 个属、{dev["workloads"]} 个工作负载；旧验证半区也归入开发，不能再当作新证据。',
        f'- 验证：**{val["genomes"]} 个新基因组、{audit["species"]} 个物种、{val["genera"]} 个留出属**，覆盖 {audit["families"]} 个科、{audit["phyla"]} 个门；accession/配对别名、BioSample、序列指纹与全部历史输入无重叠，属与本轮 60 个开发属无重叠。',
        '- 范围更正：这 30 个属都在更早的大规模面板中出现过，因此不能称为整个项目从未使用过的属。最初登记措辞过宽；`historical_genus_audit.json` 与 `REPORTING_NOTES.md` 明确记录了这一点。没有因此更换验证数据或重新训练规则。',
        f'- 验证查询来自 {val["query_sources"]} 个独立来源基因组，对照库为 {val["reference_genomes"]} 个其他基因组。每个属是 9 个查询来源、6 个参考；每种条件得到一个工作负载。',
        f'- 10 种原条件 + 4 种未见条件，共 {val["workloads"]} 个验证工作负载；k=20/24/32，pair payload=512 B/2 KiB/8 KiB/32 KiB。',
        f'- {len(inventory)} 个原生作业、**{int(counts):,} 条算法/配置/配对测量**，另有 {int(exact_pairs):,} 条精确集合配对结果。重复配置、种子、派生查询不增加生物学独立样本数；区间的独立单位仍为 30 个验证属。',
        '- 所有读段均为真实组装的受控派生：均匀位置采样、均匀替换错误、保留 record 边界。这不是实测测序平台数据，也不覆盖蛋白质或全局未知物种数据库检索。',
        '- 100 kb 是目标片段长度：开发集 540 条中 6 条、验证集 270 条中 2 条因 contig 长度限制截短，最短分别为 59,242/70,301 bp。所有真值及 pilot 特征均使用实际序列；实际输入规模见 `input_condition_summary.csv`。',
        '- 物种标签来自冻结的 NCBI 元数据。验证集无名称为 sp./uncultured/unclassified/bacterium 的分组；开发集含 1 个这样的 TaxID 分组（5 个基因组）。标签检索成绩不等同于经过比对确认的 ANI。',
        '',
        '## 先看可以使用的场景候选',
        '',
        '下表固定 k=24、最终 pair payload≤8 KiB。每个算法的参数由开发数据选择；候选保留了不能以本场景内同时区间和实际效应门槛排除的算法。'
        '不同表格单元之间没有全局多重比较保证，不能把描述性第一名写成每个场景唯一赢家。',
        '“†”表示该场景仍有一部分工作负载在整个候选网格中都没有满足标准目标的配置；灰格/“全部未达标”需要回退或放宽约束，不能推荐最低 MAE 的算法作为可用解。',
        '',table(cand,{'condition':'输入条件','jaccard':'Jaccard','containment':'查询 containment','cardinality':'查询基数','neighbor':'最近邻恢复','species':'物种检索','clustering':'聚类'}),
        '',
        '![场景候选](figures/scenario_candidates_k24_8k.png)',
        '',
        'SetSketch-fixed 始终指隔离实验中的修正版；不能把结果直接归于未修正的生产版本。'
        'KSSD/无权 ProbMinHash 的 containment 和 cardinality 被明确视为不支持；BinDash 的可用性还受其原生基数估计影响。',
        '',
        '## 冻结选择是否真正改善结果',
        '',
        '基线是在开发集按任务、输入 kind、k 和预算选出的最强固定配置，不是随意挑一个较弱算法。'
        '下表为 k=24/8 KiB，全部 14 种验证条件按属等权汇总。损失越小越好：J/C 为绝对误差；基数为截断到 1 的相对误差；'
        '最近邻/物种为 1−期望并列评分；聚类为 (1−ARI)/2。不可用或超预算结果损失为 1。不同任务的损失不能混成一个总分。',
        '',table(focus,{'metric':'任务','policy':'选择策略','loss':'失败惩罚后损失','availability':'可用且预算内比例','mean_regret':'相对事后最小损失的遗憾'}),
        '',table(main,{'metric':'任务','fixed_loss':'固定方案损失','policy_loss':'pilot 规则损失','difference':'规则−固定','simultaneous_low':'同时区间下界','simultaneous_high':'同时区间上界','material_improvement':'达到实际收益门槛'}),
        '',
        '![规则和固定方案](figures/advisor_vs_fixed_simultaneous.png)',
        '',
        '全部 72 项以及各场景的配对区间见 `advisor_vs_fixed_intervals.csv`。场景区间是逐项区间；'
        '主要汇总区间才在每个选择器的 72 项对比内同时校正。随机种子先在属内平均，区间不覆盖所有可能哈希种子的额外不确定性。',
        '',
        '## 预测达标不等于实际达标',
        '',
        '标准精度要求：J MAE≤0.02/P95≤0.05；C≤0.05/0.10；基数相对误差≤0.10/0.25；最近邻恢复≥90%；物种准确率≥95%；聚类 ARI≥0.80。'
        '同时要求至少 99% 可用且在预算内的输出、最多 1% 超预算输出。P95 先在每个工作负载内计算；预测取邻近工作负载 P95 的平均，属于经验预测而非置信上界。',
        '',table(chosen_profiles,{'metric':'任务','recommendation_coverage':'预测达标覆盖率','attainment_rate':'实际达标率','false_feasible_rate':'预测达标中的错误比例',
            'selection_miss_rate':'有可行方案时选错的比例','no_measured_feasible_rate':'所有已测配置均未达标的比例'}),
        '',
        '![实际精度达标](figures/precision_attainment_k24_8k.png)',
        '',
        '所选配置不达标，不等于这个预算下没有可行算法。表中分别检查“存在达标配置却选错”和“全部已测配置都未达标”；后者也只限本次参数网格。'
        '当前规则优化失败惩罚后的平均损失，精度门槛用于推荐后的筛查，还不是直接以 MAE/P95 约束优化的选择器。'
        '不能删除空 Sketch、超预算或选错的工作负载后再报好看的 MAE。',
        '',
        '## 未见输入条件的迁移',
        '',
        '四种未见条件是：5 kb 片段；30% 完整度/10% 污染的 tile 组装；150 bp、0.5×、2% 替换读段；250 bp、0.5×、1% 替换读段。'
        '它们在模型冻结前单独登记，只用于验证，没有重新拟合规则或阈值。下表的原条件与未见条件有不同难度分布，应比较各自内部“规则−固定”，不能直接把两组绝对成绩当作同难度比较。',
        '',table(transfer,{'metric':'任务','condition':'条件组','difference':'规则−固定','pointwise_low':'逐项下界','pointwise_high':'逐项上界'}),
        '',
        '## 短序列的零估计陷阱与 pilot 局限',
        '',table(zero,{'condition':'输入','conditional_mae':'J 条件 MAE','constant0_mae':'恒定输出 0 的 MAE','positive_pair_relative_mae':'正真值配对相对 MAE','zero_fraction':'J 输出 0 的比例'}),
        '',
        '小片段对完整基因组可以具有很高的查询 containment，同时 Jaccard 很低。只满足绝对 J 误差并不证明能保持排序、区分近缘对象或准确估计包含率。'
        'K=32 pilot 在低交集概率下会频繁输出 0，所以不能把 pilot J=0 理解为两个序列完全不相关。',
        '',
        '![pilot 分辨率](figures/pilot_resolution_limits.png)',
        '',
        'pilot 每对象最多 256 B，9+6 对象工作负载最多 3,840 B 临时查询载荷。实验在构建最终 Sketch 前释放 pilot；'
        '同时保留时必须另计空间。没有评测其输入扫描和构建时间，也没有证明选择本身有端到端性能收益。',
        '',
        '## 下游任务与开放集合失败',
        '',
        '最近邻参照是精确 k-mer Jaccard 的最优参考；物种评价按来源 TaxID；聚类使用仅在开发数据拟合的单链接阈值。'
        '三个任务不能用同一个 MAE 排名互相代替。每个属只有三个候选物种；不能外推为任意大参考库的分类能力。',
        '短片段来自不同基因组位置，查询间可能没有同源重叠；本次单链接阈值网格下，精确集合的 1/5/10 kb 片段 ARI 也接近 0。'
        '这属于集合 Jaccard 加当前聚类流程的任务限制，不能把它全归因于 Sketch 误差。',
        '',
        '![下游任务](figures/downstream_tasks_k24_8k.png)',
        '',
        '开放集合测试移除查询物种的参考，再检查是否错误接受其他物种。接受阈值只在开发集拟合，目标是 90% 已知查询召回；'
        '同时列出实际已知接受率和未知错误接受率。全为零的相似度会导致并列，阈值为零时也可能错误接受未知查询。',
        '',table(open_set,{'condition':'输入','policy':'策略','known_acceptance':'已知接受率','unknown_false_acceptance':'未知错误接受率','open_availability':'已知库可用率','unknown_availability':'移除物种后可用率'}),
        '',
        '## 可复现的规则与能力边界',
        '',
        '实验性的 `advisor.py` 返回候选及经验损失预测。它只使用输入 kind、k、碱基/record 统计以及独立 seed=919 的 K=32 pilot。'
        '真实 J/C/基数、物种/属编号、生成覆盖度、错误率和污染比例不参与选择。每个查询工作负载选一套兼容的 backend/参数，不能给已建好的索引临时换算法后跨算法比较。',
        '',table(fixed,{'kind':'输入 kind','metric':'任务','config':'开发集冻结的固定基线配置'}),
        '',
        '加权 ProbMinHash 的 probability Jaccard、min-max weighted Jaccard 和普通集合 Jaccard 不同。'
        '本轮七算法比较只覆盖普通集合语义；加权概率 Jaccard 继续使用上一轮独立加权结果。OrderMinHash 估计顺序相似性，也不应插入这张集合 J 排名表。'
        'Legacy MinHash 仍是传统实现基线。不能从不同估计对象之间的能力差异宣称某算法统计上全面占优。',
        '',
        '## 文件与审计',
        '',
        '- `scenario_candidates.csv` / `algorithm_contrasts.csv`：全部 k/预算/场景的候选与场景内同时区间。',
        '- `validation_decisions.csv.gz`：验证测量之前冻结的 302,400 条决定；`decision_outcomes.csv.gz`：决定与实际结果逐项对应。',
        '- `validation_native_scores.csv.gz` / `validation_native_tasks.csv.gz`：全部候选配置的评分；失败没有删除。',
        '- `precision_profile_summary.csv`：标准/严格门槛下的覆盖率、达标率与错误推荐率。',
        '- `model/`、`example_features.json`、`example_recommendation.json`：可复查的模型和可运行示例。',
        '- `independence_audit.json`、输入/派生来源清单、`native_archive_manifest.csv`：独立性、源文件与全部原生归档的 SHA256。',
        '- `shared_exact_verification.csv`：与旧实验共享的 900 个基因组和每个 k 的 6,300 个无序配对，精确基数/交集全部一致。',
        '- `analysis_tests.log`、`artifact_sha256.json`：分析测试及公开产物校验；完整 FASTA/原生归档仍保存在外部实验目录。',
        '',
        '本轮的结论范围是固定实现、细菌组装及上述受控派生、三种 k、四档查询载荷，以及已知属候选库。'
        '优化实现本身的性能贡献留待统一性能实验；这里没有加入标量对比或速度排名。'
    ]
    return '\n'.join(lines)+'\n'


def main():
    if not (OUT/'evaluation_complete.json').exists() or not (OUT/'figures_complete.json').exists():
        raise ValueError('complete evaluation and figures before publication')
    input_integrity();input_condition_summary();shared_exact_check();inventory=native_inventory()
    if not (OUT/'historical_genus_audit.json').exists():raise ValueError('complete the full historical genus-scope audit')
    # Include every configuration so the public data can audit all alternatives.
    save_csv(OUT/'validation_native_scores.csv.gz',read_scores('validation'))
    tasks=pd.concat([pd.read_csv(p) for p in sorted((OUT/'validation/tasks').glob('*.csv.gz'))],ignore_index=True)
    save_csv(OUT/'validation_native_tasks.csv.gz',tasks)
    (OUT/'report.md').write_text(report_text(inventory))
    (OUT/'README.md').write_text(publication_readme())
    features=pd.read_csv(OUT/'validation_features.csv')
    example=features[(features.k==24)&(features.condition=='full')].sort_values('workload_id').iloc[0]
    atomic(OUT/'example_features.json',{c:float(example[c]) for c in PILOT_FEATURES})
    test=subprocess.run([sys.executable,'-m','pytest','-q',str(HERE/'test_selection.py')],capture_output=True,text=True)
    (OUT/'analysis_tests.log').write_text(test.stdout+test.stderr)
    if test.returncode:raise RuntimeError('analysis verification tests failed')
    PUBLISHED.mkdir(parents=True,exist_ok=True)
    for p in OUT.iterdir():
        if p.is_file() and (p.suffix in ['.csv','.gz','.json','.md'] or p.name=='analysis_tests.log'):
            shutil.copyfile(p,PUBLISHED/p.name)
    for folder in ['figures','model']:
        shutil.copytree(OUT/folder,PUBLISHED/folder,dirs_exist_ok=True)
    truth_destination=PUBLISHED/'truth_and_pilot';truth_destination.mkdir(exist_ok=True)
    for pattern in ['exact_*.csv.gz','pilot_*.csv.gz']:
        for p in (OUT/'validation').glob(pattern):shutil.copyfile(p,truth_destination/p.name)
    for filename in ['PROTOCOL.md','ANALYSIS_ADDENDUM.md','REPORTING_NOTES.md','config.json','validation_extra_conditions.json']:
        shutil.copyfile(HERE/filename,PUBLISHED/filename)
    for filename in ['independence_audit.json','registration.json','development_preparation.json','validation_preparation.json',
        'validation_extension_registration.json','validation_extension_complete.json','development_genomes.tsv','validation_genomes.tsv',
        'development_objects.tsv','validation_objects.tsv','development_pairs.tsv','validation_pairs.tsv','validation_workloads.json']:
        shutil.copyfile(DATA/filename,PUBLISHED/filename)
    for filename in ['model_frozen.json','decisions_frozen.json','input_sha256.json']:
        shutil.copyfile(RUN/filename,PUBLISHED/filename)
    shutil.copyfile(ROOT/'build_comprehensive_20261030/build_manifest.json',PUBLISHED/'native_build_manifest.json')
    shutil.copyfile(ROOT/'runs/comprehensive_accuracy_20261030/identity.json',PUBLISHED/'native_source_identity.json')
    atomic(PUBLISHED/'workflow_sha256.json',{p.name:sha(p) for p in sorted(HERE.iterdir())
        if p.is_file() and p.suffix in ['.py','.md','.json']})
    result=subprocess.run([sys.executable,str(HERE/'advisor.py'),'--stage','recommend','--model',str(PUBLISHED/'model'),
        '--features-json',str(PUBLISHED/'example_features.json'),'--k','24','--kind','assembly','--metric','jaccard','--budget','8192'],
        capture_output=True,text=True,check=True)
    recommendation=json.loads(result.stdout);atomic(PUBLISHED/'example_recommendation.json',recommendation)
    frozen=pd.read_csv(OUT/'validation_decisions.csv.gz')
    chosen=frozen[(frozen.workload_id==example.workload_id)&(frozen.k==24)&(frozen.metric=='jaccard')&
                  (frozen.budget==8192)&(frozen.policy=='pilot')].iloc[0]
    if chosen.config!=recommendation['config']:raise ValueError('portable model differs from frozen choice')
    manifest={str(p.relative_to(PUBLISHED)):sha(p) for p in sorted(PUBLISHED.rglob('*'))
              if p.is_file() and p.name not in ['artifact_sha256.json','publication_verification.json']}
    atomic(PUBLISHED/'artifact_sha256.json',manifest)
    print('Published artifacts',len(manifest),'at',PUBLISHED,flush=True)


if __name__=='__main__':main()
