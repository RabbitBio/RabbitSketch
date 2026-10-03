#!/usr/bin/env python3
"""Data-derived tables for the final robustness addendum; retain all failures."""
import pandas as pd
import numpy as np
from common import *
from analyze import NAMES
REPO=HERE.parent.parent
def load(name):return pd.read_csv(ANALYSIS/(name+'.csv'),low_memory=False)
def number(x,d=5):
    if pd.isna(x):return '—'
    return f'{x:.{d}f}'
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(v).replace('|','/') for v in r)+' |' for r in rows])
def link(name):return f'[{name}]({ANALYSIS/name})'
def main():
    prep=json.loads((DATA/'natural_preparation.json').read_text());sel=json.loads((DATA/'natural_selection.json').read_text());prov=json.loads((ANALYSIS/'data_provenance.json').read_text());audit=json.loads((ANALYSIS/'raw_audit.json').read_text())
    s=load('budget_accuracy');ci=load('paired_cluster_intervals');cost=load('build_summary');model=load('mixed_model_summary');advisor=load('frozen_advisor_external');align=load('alignment_reference_comparison');alignraw=load('alignment_results');exact=load('exact_truth_and_models')
    def result(panel,alg,stratum,metric='j_error',budget=8192,k=24):
        g=s[(s.panel==panel)&(s.algorithm==alg)&(s.stratum==stratum)&(s.metric==metric)&(s.budget==budget)&(s.k==k)]
        return None if g.empty else g.iloc[0]
    def value(panel,alg,stratum,metric='j_error',budget=8192):
        r=result(panel,alg,stratum,metric,budget);return np.nan if r is None else r.mae
    mh=value('natural','minhash','all');bd=value('natural','bindash','all');sshigh=value('natural','setsketch','J>0.5');bdhigh=value('natural','bindash','J>0.5')
    bc=cost[(cost.panel=='natural')&(cost.k==24)].set_index(['algorithm','parameter'])
    buildratio=bc.loc[('bindash',2048),'median_ready']/bc.loc[('minhash',512),'median_ready']
    naturalad=advisor[(advisor.panel=='natural')&advisor.status.isin(['pass','fail'])]
    text=[]
    def add(x):text.append(x.strip()+'\n')
    add(f'''# RabbitSketch2：独立多物种与混合突变验证，及论文结论收口

本报告衔接 [逐算法详细分析]({REPO/'RabbitSketch2_algorithm_selection_analysis.md'})，只补两项已授权缺口：独立真实 genomes，以及 substitution＋indel。未寻找第二台 CPU，未新增大规模性能矩阵，未扩展 HLL、containment、weighted 参数网格，未恢复 OrderMinHash。当前 MinHash 始终指 FastKMV；KSSD 保留。

**这一轮的目的，是检验旧结论能够迁移到什么范围，而不是再次从新数据中选出最有利参数。** 配置、抽样规则、指标及对照在新 sketch 结果生成前固定，见 [PROTOCOL.md]({HERE/'PROTOCOL.md'})。SetSketch 继续使用上一轮隔离副本中的 partitioned-joint 修正版，生产源码未被这个补丁替换。

完成规模：独立真实数据 **{prep['genomes']} genomes / {prep['species']} species / {prep['genera']} genera / {prep['families']} families / {prep['phyla']} phyla，{sel['pairs']:,} pairs**；混合突变 **20 references、320 pairs**。共完成 **{audit['jobs']} 个 native/exact 批次，{audit['rows']['pair']:,} 条 pair 记录**，另有双向 FastANI 与 {len(alignraw)} 对预先选定的 NUCmer 审核。sketch 与 exact 测量沿用当前 main 批次入口，单个 benchmark 进程、64 OpenMP threads 顺序运行；参考工具另外顺序执行。

## 1. 先读这四个验证结果

在新真实数据、k=24、8 KiB/pair 下，BinDash 的整体 exact-Jaccard MAE / MinHash MAE 为 **{bd/mh:.3f}**；SetSketch 在 J>0.5 区间的 MAE / BinDash MAE 为 **{sshigh/bdhigh:.3f}**。这两个比值描述不同条件，不能直接合成一个全局优胜顺序。

900-genome 构建的对应比较中，BinDash2048 / MinHash512 的 build＋query-preparation 批次时间比为 **{buildratio:.3f}**。这是同一 CPU 上对新数据的构建复核；没有重新估计或推广旧的 all-vs-all crossover。

被冻结的 Mash-family advisor 在这批真实数据上有 **{len(naturalad)} 个可评估的既定选择，其中 {int((naturalad.status=='pass').sum())} 个通过、{int((naturalad.status=='fail').sum())} 个失败**。请求可能共享配置及 pairs，这不是独立 Bernoulli 实验。完整原因及 exact-J 转换残差见第 7 节。

这些观察分别支持或限制“统计精度方向”“实现成本方向”“具体误差阈值迁移”。尤其不能用第一项的成立，自动宣称第三项也成立。''')
    add('## 2. 独立数据是怎样构成的')
    add(f'''每属取 3 个 species taxid，每 species 取 5 个不同 BioSample。先按 phylum/family 平衡并用固定哈希决定候选顺序，再做 QC：1–10 Mb、完整 genome representation、≤300 contigs、≥99% ACGT、N50≥10 kb。这里“完整 representation”不等于所有 assembly 都是闭合染色体。

旧实验中识别到的 {prep['prior_sequence_fingerprints']:,} 个不同 sequence fingerprints 被用于排重；与旧输入的 accession roots、BioSample、相同序列 fingerprint 重叠均为 0。指纹忽略 contig 次序和正反链、保留 contig 边界，非 ACGT 标准化为 N。仍有 **{prov['overlap_with_earlier_inputs']['species']} 个 species、{prov['overlap_with_earlier_inputs']['genera']} 个 genera** 与旧数据的分类学范围重叠；本实验验证的是独立 assemblies/BioSamples，不宣称所有分类群均未见过，也不排除不同 BioSample 的近克隆。

{prov['unnamed_or_placeholder_species']} 个 species taxid 的 scientific name 含 sp./bacterium 等未命名或占位形式。因此“{prep['species']} species”是按 NCBI species taxid 统计，不等于 {prep['species']} 个已正式命名物种。来源 snapshot 的 SHA256 和文件时间见 {link('data_provenance.json')}；文件时间不冒充下载或发布日期。

首先计算同属候选及远缘 controls 的 FastANI，再在不读取任何 sketch 误差的条件下选择最终 pairs：每物种 4 对覆盖样本，各主要 ANI bin 最多补入 600 对，保留可用跨属/科 controls。它是有意覆盖场景的 benchmark，整体 MAE 不是 RefSeq 自然分布下的总体误差估计。''')
    rows=[[b,sel['candidate_bin_counts'].get(b,'—'),sel['final_bin_counts'].get(b,0)] for b in ['90-95','95-99','99-100','below90','missing']]
    add(f'''如果只与先前 **accuracy panel** 比较，而不计入 10,000-genome 性能集，本轮只有 {prov['overlap_with_earlier_accuracy_panel']['species']} 个 species、{prov['overlap_with_earlier_accuracy_panel']['genera']} 个属曾出现过；即新增 {prep['species']-prov['overlap_with_earlier_accuracy_panel']['species']} 个 species taxid 和 {prep['genera']-prov['overlap_with_earlier_accuracy_panel']['genera']} 个属。因而它扩大了 accuracy 的分类学覆盖，但不是严格只由未见属构成的 leave-genus-out 测试。CSV 中另外保留“两个端点均来自未见 accuracy 属”的描述性分层，不能把这个补充分层的显著性算作预先指定的独立主结果。

最终远缘 controls 是 {prep['controls']['cross_genus_same_family']} 对同科跨属、{prep['controls']['cross_family']} 对跨科。前者低于计划最多 120 对，因为本轮 60 属分散在 55 科，仅取每属一个 control representative；这个缺口没有通过改变抽样策略补齐。三个 ANI 分层与物种覆盖面板取并集，所以最终某个 bin 可能超过额外抽入的 600 对。''')
    add(table(['双向 FastANI 区间 (%)','候选 qualified pairs','最终 pairs'],rows))
    add(f'''FastANI 要求两个方向都输出且 aligned fraction 均 ≥0.2，使用两方向均值。没有输出的 pair 保持缺失；仍计算其 exact Jaccard，不将 ANI 设成 0。候选不足的 bin 不另挑数据去填出期望效果。完整分类组成见 {link('natural_composition.csv')}，每个 genome 的来源和 QC 见 {link('natural_taxonomy.csv')}。

## 3. 固定配置与公平口径

主分析 k=24：MinHash K=128/512/2048；BinDash bins=128/512/2048/8192（16 bit）；SetSketch p=12；FracMinHash scaled=5000/10000；KSSD dr3/half_subk6 和 dr4/half_subk7。次分析 k=20：MH512、BD2048、SS12、Frac10000、KSSD3。除 KSSD42 外，种子为 42/101/202。两种 k 均计算完整 canonical 2-bit k-mer 集合。

2/8/32 KiB 是 A+B 两个 sketch 的 payload 总和。固定大小的 MH/BD 等字节对应为 128/512、512/2048、2048/8192。Frac/KSSD 只在全部测试 seeds 的实测 pair bytes 均不超预算时纳入，并按字节选可容纳的最大已测档位，完全不看该 pair 的误差。这个预算 envelope 不等价于部署一个固定参数；advisor 复核则严格固定参数。SetSketch 只验证 p=12，因此它在 32 KiB 预算下仍只使用 8 KiB/pair。

ANI 使用当前接口对应的转换 $g(J)=(2J/(1+J))^{{1/k}}$。所有 ANI 误差表以百分点（pp）表示；Jaccard 误差为 [0,1] 上绝对误差。P95 是误差分位数，95% CI 是另一个量。将误差分为 $|\\hat J-J|$、$|g(\\hat J)-g(J)|$ 和 $|g(\\hat J)-ANI_{{ref}}|$，另列 $|g(J)-ANI_{{ref}}|$。后者是转换模型与参考定义之间的差异，不是 sketch 采样误差。

误差均值/P95 包括全部 seeds；配对 CI 先在 pair 内平均 seeds，再按真实数据的属、或受控数据的 reference 重采样 2,000 次。跨属边采用两端 cluster 权重，同属边只取一次 cluster 权重。区间是 conditional on 测过的 seeds、逐项且未校正多重比较；不能据此声称整个配置空间被统计淘汰。''')
    add('## 4. 新真实 genomes：逐算法的数据优势')
    rows=[]
    for bin_ in ['90-95','95-99','99-100']:
        for alg in NAMES:
            j=result('natural',alg,bin_);p=result('natural',alg,bin_,'sampling_error');a=result('natural',alg,bin_,'fastani_error')
            if j is not None:rows.append([bin_,NAMES[alg],int(j.pairs),number(j.mae),number(j.within_cluster_macro_mae),number(100*p.mae,4),number(100*a.mae,4),number(100*a.p95,4)])
    add('k=24，8 KiB/pair。变长 sketch 的有效 pair 覆盖必须与精度一起读：')
    add(table(['FastANI (%)','算法','pairs','J MAE','属内宏平均 J MAE','ANI sampling MAE (pp)','总 ANI MAE vs FastANI (pp)','总 ANI P95 (pp)'],rows))
    add(f'原始分层统计：{link("budget_accuracy.csv")}。属内宏平均先求每属均值再平均，减少拥有更多 selected pairs 的属对结论的影响；远缘跨属 pairs 不纳入这个宏平均，但保留在 pair-weighted 结果中。')
    rows=[]
    for budget in [2048,8192,32768]:
        for stratum in ['all','J>0.5']:
            mhv=value('natural','minhash',stratum,budget=budget);bdv=value('natural','bindash',stratum,budget=budget);ssv=value('natural','setsketch',stratum,budget=budget)
            rows.append([budget//1024,stratum,number(mhv),number(bdv),number(bdv/mhv,3),number(ssv),number(ssv/bdv,3)])
    add(table(['KiB/pair','范围','MH J MAE','BD J MAE','BD/MH','SS J MAE','SS/BD'],rows))
    add('“all” 的 Jaccard 统计包含全部 1,967 对，包括 FastANI 未输出的远缘 controls。CSV 中的 overall sampling error 也能在这些 pair 上计算，但 overall 总 ANI 误差只对有 reference 的 pairs 定义；两列不能在不同有效样本集合下直接相减。主要 ANI 比较以上面三个 qualified ANI bins 为准。')
    rows=[]
    wanted=ci[(ci.panel=='natural')&(ci.k==24)&(ci.budget==8192)&((ci.left=='bindash')&(ci.right=='minhash')&ci.stratum.isin(['all','90-95','95-99','99-100'])|(ci.left=='setsketch')&(ci.right=='bindash')&ci.stratum.isin(['J>0.5','99-100']))]
    for r in wanted.itertuples():
        if r.metric not in ['j_error','sampling_error','fastani_error']:continue
        scale=1 if r.metric=='j_error' else 100
        verdict='左侧误差更低，CI 不跨 0' if r.ci_high<0 else ('右侧误差更低，CI 不跨 0' if r.ci_low>0 else 'CI 跨 0，方向未分离')
        rows.append([r.stratum,f'{NAMES[r.left]} − {NAMES[r.right]}',r.metric,int(r.pairs),int(r.clusters),number(scale*r.difference),f'[{number(scale*r.ci_low)}, {number(scale*r.ci_high)}]',verdict])
    add(table(['范围','配对差','误差指标','pairs','clusters','平均差','95% CI','解释'],rows))
    add(f'''误差差值为左减右，负值有利于左侧。J 指标仍为 fraction，ANI 指标在表内换算为 pp。这个表比“总平均第一名”更适合决定论文应写到什么强度。完整 k20/24、2/8/32 KiB 区间见 {link('paired_cluster_intervals.csv')}。

MinHash 的存在理由仍应由满足精度后的构建成本支撑；不能把它称为每字节精度的默认最优。BinDash 的推荐需要同时得到等字节 J 和 sampling error 的支持，不能用总 ANI 误差中的偶然模型抵消代替统计精度。SetSketch 的 high-overlap 优势须在相应分层读取，不能升级成所有真实 ANI 最准。FracMinHash 与 KSSD 的变长 payload 需要先检查覆盖率：缺少可行配置不是一个应被删掉的失败样本。''')
    cov=load('budget_coverage');rows=[]
    for r in cov[(cov.panel=='natural')&(cov.k==24)&(cov.budget==8192)].itertuples():rows.append([NAMES[r.algorithm],int(r.budget_feasible_pairs),int(r.total_pairs),f'{100*r.budget_coverage:.2f}%',int(r.empty_observations)])
    add(table(['算法','空间可行 pairs','所有 pairs','预算覆盖率','可行输出中的 empty 次数'],rows))
    add('这一轮 FracMinHash 只保留 scaled=5000/10000，因此 8 KiB 下的约 80% 覆盖率是**冻结配置集合的覆盖率**，不是 FracMinHash 本身无法支持其余 genomes。KSSD 的 reduction 档位也较粗，8 KiB envelope 中常只使用远少于 8 KiB 的 payload；不能根据它在这个表中较大的误差宣称严格等字节全局劣势。逐配置空间和误差必须一起读。')
    add('### 4.1 构建优势在新输入上是否成立')
    qc=load('pair_query_cost');qc=qc[(qc.panel=='natural')&(qc.k==24)].set_index(['algorithm','parameter'])
    rows=[]
    for alg,p in [('minhash',512),('bindash',2048),('setsketch',12),('fracminhash',5000),('fracminhash',10000),('kssd',3),('kssd',4)]:
        r=bc.loc[alg,p];rows.append([NAMES[alg],p,int(r.replicates),number(r.median_build,4),number(r.median_ready,4),f'{r.min_build:.4f}–{r.max_build:.4f}',number(r.median_payload/prep['genomes'],1),number(qc.loc[(alg,p),'median_j_query_seconds']*1e6,3)])
    add(table(['算法','native 参数','seed 批次','构建含 setup (s)','构建＋准备 (s)','构建 min–max (s)','平均 B/object','pair elapsed 中位数 (μs)'],rows))
    add(f'''这是整个 {prep['genomes']}-genome batch 的 wall time，包含 KSSD 的共享初始化；输入解析和 pair 查询不混入构建。三个 seed 批次不是大量独立性能重复，KSSD 只有一批。系统运行时波动与种子差异不能由这张表完全分离。完整原始批次和 input time 在 {link('build_batches.csv')}；本轮每对查询计时在 {link('pair_query_cost.csv')}，后者仍是 64 线程负载中的 elapsed time，不能冒充独占 pair latency 或新的 all-vs-all crossover。

旧的 N*≈8593/2177 和实际 all-vs-all 证据保留原有机器及 workload 限定。本轮没有第二 CPU，论文不能增加 across-architectures 的结论。''')
    add('### 4.2 k=20 的预定敏感性复核')
    rows=[]
    for k in [20,24]:
        m=result('natural','minhash','all',k=k);b=result('natural','bindash','all',k=k);hb=result('natural','bindash','J>0.5',k=k);hs=result('natural','setsketch','J>0.5',k=k)
        rows.append([k,number(m.mae),number(b.mae),number(b.mae/m.mae,3),int(hb.pairs),number(hb.mae),number(hs.mae),number(hs.mae/hb.mae,3)])
    add(table(['k','all MH J MAE','all BD J MAE','BD/MH','J>0.5 pairs','高 J BD MAE','高 J SS MAE','SS/BD'],rows))
    add('两种 k 下 BinDash 的等字节优势和 SetSketch 的高 J 改善方向均保留，幅度不同。J>0.5 是各自 k 下的 truth 分层，所以两个 k 的高 J 行不是完全相同的 pair 集合；不能把差额直接解释成 k 的因果收益。')
    add('## 5. 混合 substitution/indel：什么保持稳定，什么会改变')
    add('20 个真实 5 Mb references 来自上一轮，只用于配对的 mutation-model sensitivity，不计为第二份独立生物验证。每个 reference 有 4 个目标 alignment-column identities × 4 种 affected-base composition，共 16 个变体。indel 的 insertion/deletion 碱基数相等，事件长度为截断 geometric(p=0.5)，上限 10；事件不重叠，不修改模糊碱基。composition 以受影响碱基数统计，不是 event 数。\n\n每个变体保存可重放的 compressed edit script，并精确记录 M、S、I、D。目标定义为 M/(M+S+I+D)，同时保留 M/(M+S) 的 ungapped homologous identity。脚本按独立的字符串实现测试过重放，真实 320 个输出也逐一通过保存后重放、长度及计数验证。它们是已知构造真值；FastANI 是另外一列参考。')
    rows=[]
    for r in model[model.k==24].itertuples():rows.append([f'{r.target_identity:.0%}',f'{r.indel_share:.0%}',number(r.mean_j),number(100*r.homologous_identity,4),number(100*r.fastani,4),number(100*r.column_model_mae,4),number(100*r.homologous_model_mae,4),number(100*r.fastani_model_mae,4)])
    add(table(['column target','indel bases share','mean exact J','ungapped identity (%)','FastANI (%)','exact-J→column MAE (pp)','exact-J→ungapped MAE (pp)','exact-J→FastANI MAE (pp)'],rows))
    add('这个表中的差异不随 sketch 分辨率自动消失。对于 substitution、indel 事件和同源/非同源内容，Jaccard→ANI 的映射前提不同。增加 sketch bytes 主要减少 sampling error，不能把它当作修复任意 ANI 参考偏差的方法。')
    rows=[]
    for target in [90,95,97,99]:
        for q in [0,10,25,50]:
            st=f'id{target}_indel{q}';vals=[100*value('mixed',a,st,'sampling_error') for a in NAMES]
            ratio=value('mixed','bindash',st)/value('mixed','minhash',st)
            rows.append([target,q,*[number(v,4) for v in vals],number(ratio,3)])
    add('k=24、8 KiB/pair；下表按同一 reference/target 直接比较。除最后一列为 Jaccard MAE 比值外，算法列均为 ANI sampling MAE（pp）：')
    add(table(['column target (%)','indel (%)',*NAMES.values(),'BD/MH J MAE'],rows))
    mci=ci[(ci.panel=='mixed')&(ci.k==24)&(ci.budget==8192)&(ci.left=='bindash')&(ci.right=='minhash')]
    rows=[]
    for metric in ['j_error','sampling_error','column_error','homologous_error','fastani_error']:
        x=mci[mci.metric==metric];rows.append([metric,len(x),int((x.difference<0).sum()),int((x.ci_high<0).sum()),int((x.ci_low>0).sum()),int(((x.ci_low<=0)&(x.ci_high>=0)).sum())])
    add(table(['BD − MH 指标','已测 target×composition cells','平均差有利 BD','逐项 CI 有利 BD','逐项 CI 有利 MH','CI 跨 0'],rows))
    add(f'''这些是冻结的 16 个 cells 的方向核查，不能把 16 个相关 cells 当 16 个独立证据倍增信心。支持“robust”的正文应指明是 Jaccard / sampling 精度，还是对某个 ANI 参考的总误差；若总误差方向逆转，保留该结果。补充按 observed exact J 分层的结果在 {link('mixed_observed_j_accuracy.csv')}，用于分辨 mutation composition 改变了 J 分布之后的效应，不作为匹配了所有变量的因果实验。''')
    add('## 6. NUCmer 参考检查：避免由单一 ANI reference 决定胜负')
    rows=[]
    for panel in ['natural','mixed']:
        raw=alignraw[alignraw.panel==panel];g=align[(align.panel==panel)&(align.k==24)]
        rows.append([panel,len(raw),int((raw.status=='ok').sum()),len(g),number(100*abs(g.fastani_nucmer_difference).mean(),4),number(100*abs(g.fastani_nucmer_difference).max(),4)])
    add(table(['面板','预选 pairs','命令成功','1-to-1 AF 双向≥0.2','FastANI/NUCmer 平均绝对差 (pp)','最大绝对差 (pp)'],rows))
    astats=load('alignment_sketch_accuracy');rows=[]
    for bin_ in ['90-95','95-99','99-100']:
        for alg,p in [('minhash',512),('bindash',2048),('setsketch',12)]:
            g=astats[(astats.panel=='natural')&(astats.k==24)&(astats.ani_bin==bin_)&(astats.algorithm==alg)&(astats.parameter==p)].set_index('metric')
            rows.append([bin_,NAMES[alg],int(g.loc['nucmer_error','pairs']),number(100*g.loc['fastani_error','mae'],4),number(100*g.loc['nucmer_error','mae'],4)])
    add('同一审核子集、k=24、8 KiB/pair，分别更换 ANI reference（均为 MAE，pp）：')
    add(table(['FastANI 分层 (%)','算法','审核 pairs','对 FastANI','对 NUCmer'],rows))
    add(f'''自然数据每 ANI bin 20 对，混合数据为预先按哈希选出的 5 个 references ×16 条件；不按 sketch 误差挑审核 pairs。NUCmer/dnadiff 使用 1-to-1 identity 及两个方向的对齐覆盖率。NUCmer 与 FastANI 衡量/筛选的对齐内容不完全相同，彼此都不是适用于所有场景的绝对真理。

逐算法在两个参考下的 MAE/P95 见 {link('alignment_sketch_accuracy.csv')}，exact-J 模型与各参考的差异见 {link('alignment_reference_comparison.csv')}。正文要把参考敏感性作为限制保留；不能只引用对偏好算法有利的参考。''')
    add('## 7. 冻结 advisor 的外部检验')
    rows=[]
    for group,g in advisor.groupby(['panel','indel_share'],dropna=False):
        ev=g[g.status.isin(['pass','fail'])];f=ev[ev.status=='fail']
        rows.append([group[0],'—' if pd.isna(group[1]) else f'{group[1]:.0%}',len(g),int((g.status=='no_frozen_candidate').sum()),int((g.status=='no_external_pairs').sum()),int((ev.status=='pass').sum()),len(f),int((f.exact_j_model_satisfies==0).sum())])
    add(table(['外部面板','indel share','固定 requests','原 calibration 无候选','本面板无对应 pairs','通过','失败','失败且 exact-J model 也超阈值'],rows))
    rows=[]
    examples=advisor[(advisor.panel=='natural')&advisor.status.isin(['pass','fail'])].drop_duplicates(['ani_bin','budget','target_mae','algorithm','parameter'])
    for r in examples.itertuples():
        rows.append([r.ani_bin,int(r.budget)//1024,f'{100*r.target_mae:.1f} / {200*r.target_mae:.1f}',f'{NAMES[r.algorithm]} {int(r.parameter)}',number(100*r.mae,4),number(100*r.p95,4),number(100*r.exact_j_model_p95,4),r.status])
    add('合并重复 workload 请求后，以下是真实数据的具体通过/失败证据（均为 pp）：')
    add(table(['ANI (%)','KiB/pair','MAE/P95 限值','冻结配置','新 MAE','新 P95','exact-J model P95','结果'],rows))
    add(f'''直接读取上一轮 `heldout_selection.csv`，不重新做 calibration、不训练 ML、不用本轮错误选择备用配置。Natural 使用原 Mash-family 请求；controlled 使用原 substitution 请求，并对四种 mutation composition 分别复核。Mixed 的区间归属用 nominal column target，防止整数碱基舍入让刚好 95% 的 pair 漂进相邻 bin；误差仍对 achieved column identity 计算。Natural 对 FastANI 的真实 bin 归属计算。

通过条件仍是 MAE≤ε、P95≤2ε、全部观测 payload≤预算。只有同时有冻结候选和对应外部 pairs 的请求才进入通过/失败分母。原 calibration 无候选的请求保留 `no_frozen_candidate`，不算成功；若已有候选但外部缺少对应 pairs，则标记 `no_external_pairs`。成本模型仍来自上一轮 10,000-genome CPU 测量，本轮不更新其 crossover。

最后一列的 exact-J residual 不是一个严格的误差下界：有限 sketch 的 bias 可能碰巧抵消 reference/model residual。但若完整集合的转换都超过约束，它已足以说明“只增加 sketch size”没有解决问题的保证。新数据失败必须限制 advisor 的适用域。实用接口应在超出已验证数据域或参考定义时显示证据不足，而非承诺 ε。所有逐请求选择、原 holdout 状态及新失败见 {link('frozen_advisor_external.csv')}。''')
    add('## 8. 每个算法的最终推荐应怎样写')
    add(table(['算法','保留的使用理由','本轮补充证据如何限制表述'],[
        ['MinHash/FastKMV','accuracy 已满足后的构建型 workload 候选','以第 4 节新输入的同预算 build time 复核；不能等同于最佳 accuracy/byte；ANI reference/model 残差需要单独达标'],
        ['BinDash','普通集合的紧凑统计估计；旧实测中 query reuse 高时有成本优势','用新真实与 mixed 的等字节 J、sampling error 和 paired CI 支撑；性能仍限同机；真实 ANI 总误差另验'],
        ['SetSketch-fixed','愿意支付查询成本以换取部分 high-overlap joint-estimation 精度','只在第 4 节支持的分层推荐；不能从 J>0.5 自动推断任意 99–100% ANI 都更好；实现版本必须标明'],
        ['FracMinHash','directed containment 且 query 有足够保留样本','本轮普通 similarity 是外部检验，不替代旧 containment 证据；变长空间与缺失可行档位必须报告'],
        ['HLL','cardinality-only 的轻量起点','本轮未重测，继续依赖上一阶段的 cardinality evidence，不将新 Jaccard 结论套给 HLL'],
        ['weighted ProbMinHash','probability Jaccard 的 abundance comparison','本轮未重测；CAMI profile 层验证与完整 reads→k-mer workflow 的边界保持不变；不是 min/max weighted Jaccard'],
        ['KSSD','提供 reduction backend；允许局部低空间/精度竞争','不因必须分配 niche 而强行设默认入口；本轮变长 payload、单 reduction seed 和含 setup 成本共同解读；未证明全局无用'],
        ['unweighted ProbMinHash / OrderMinHash','前者降低普通集合默认优先级；后者排除当前 scope','本轮没有新增这两者的证据，不宣称本轮完成它们的全面淘汰']]))
    add('论文的三条主贡献可以继续围绕 unified characterization、条件化 selection regimes、empirical constrained advisor。Advisor 的贡献应包括可明确返回无合格候选或外部验证失败，而不是保证任意 genome 都能达到输入的误差阈值。算法精度的方向性复现，与固定阈值和 crossover 的普适性，是需要分开报告的结果。')
    add('## 9. 实验冻结与交付清单')
    add(f'''本轮不再增加普通 Jaccard 参数、size ratio、weighted registers 或性能平台。当前证据仍不覆盖跨架构、reads 全流程、长 indel、rearrangement 或普通 similarity 的大型真核 genomes；这些是论文适用边界，不是本轮继续扩实验的理由。

原始输入与编辑记录：`{DATA}`。原始命令、stdout/stderr、逐批 SHA256、代码与二进制快照：`{RUN}`。分析表、PNG/PDF 与审计：`{ANALYSIS}`。所有旧阶段文件保留。

- {link('raw_audit.json')}：94 批次完整性、raw SHA256、exact/native 输入摘要一致性。
- {link('data_provenance.json')}：排重、分类组成、来源 snapshot 及限制。
- {link('configuration_accuracy.csv')}、{link('budget_accuracy.csv')}、{link('budget_coverage.csv')}：固定配置与预算比较。
- {link('paired_cluster_intervals.csv')}：全部配对差与区间。
- {link('mixed_model_summary.csv')}：编辑真值、FastANI 和 exact-J 模型差异。
- {link('frozen_advisor_external.csv')}：成功、失败及无证据请求。
- {link('independent_genome_accuracy.png')} / {link('independent_genome_accuracy.pdf')}。
- {link('mixed_mutation_sampling.png')} / {link('mixed_mutation_sampling.pdf')}。
- {link('mixed_mutation_model_residual.png')} / {link('mixed_mutation_model_residual.pdf')}。

本报告由 [write_report.py]({HERE/'write_report.py'}) 从冻结 CSV 生成表格；主要分析为 [analyze.py]({HERE/'analyze.py'})。时间值与统计误差分别按以上口径读取。最终 artifact hashes 见 {link('completion.json')}。''')
    report='\n'.join(text)
    (REPO/'RabbitSketch2_robustness_validation.md').write_text(report)
    (ANALYSIS/'report.md').write_text(report)
    print('Wrote report',len(report),'characters',flush=True)

if __name__=='__main__':main()
