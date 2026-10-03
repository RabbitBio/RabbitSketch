# RabbitSketch2：独立多物种与混合突变验证，及论文结论收口

本报告衔接 [逐算法详细分析](RabbitSketch2_algorithm_selection_analysis.md)，只补两项已授权缺口：独立真实 genomes，以及 substitution＋indel。未寻找第二台 CPU，未新增大规模性能矩阵，未扩展 HLL、containment、weighted 参数网格，未恢复 OrderMinHash。当前 MinHash 始终指 FastKMV；KSSD 保留。

**这一轮的目的，是检验旧结论能够迁移到什么范围，而不是再次从新数据中选出最有利参数。** 配置、抽样规则、指标及对照在新 sketch 结果生成前固定，见 [PROTOCOL.md](benchmarks/robustness/PROTOCOL.md)。SetSketch 继续使用上一轮隔离副本中的 partitioned-joint 修正版，生产源码未被这个补丁替换。

完成规模：独立真实数据 **900 genomes / 180 species / 60 genera / 55 families / 10 phyla，1,967 pairs**；混合突变 **20 references、320 pairs**。共完成 **94 个 native/exact 批次，107,489 条 pair 记录**，另有双向 FastANI 与 140 对预先选定的 NUCmer 审核。sketch 与 exact 测量沿用当前 main 批次入口，单个 benchmark 进程、64 OpenMP threads 顺序运行；参考工具另外顺序执行。

## 1. 先读这四个验证结果

在新真实数据、k=24、8 KiB/pair 下，BinDash 的整体 exact-Jaccard MAE / MinHash MAE 为 **0.484**；SetSketch 在 J>0.5 区间的 MAE / BinDash MAE 为 **0.816**。这两个比值描述不同条件，不能直接合成一个全局优胜顺序。

900-genome 构建的对应比较中，BinDash2048 / MinHash512 的 build＋query-preparation 批次时间比为 **1.482**。这是同一 CPU 上对新数据的构建复核；没有重新估计或推广旧的 all-vs-all crossover。

被冻结的 Mash-family advisor 在这批真实数据上有 **18 个可评估的既定选择，其中 6 个通过、12 个失败**。请求可能共享配置及 pairs，这不是独立 Bernoulli 实验。完整原因及 exact-J 转换残差见第 7 节。

这些观察分别支持或限制“统计精度方向”“实现成本方向”“具体误差阈值迁移”。尤其不能用第一项的成立，自动宣称第三项也成立。

## 2. 独立数据是怎样构成的

每属取 3 个 species taxid，每 species 取 5 个不同 BioSample。先按 phylum/family 平衡并用固定哈希决定候选顺序，再做 QC：1–10 Mb、完整 genome representation、≤300 contigs、≥99% ACGT、N50≥10 kb。这里“完整 representation”不等于所有 assembly 都是闭合染色体。

旧实验中识别到的 10,504 个不同 sequence fingerprints 被用于排重；与旧输入的 accession roots、BioSample、相同序列 fingerprint 重叠均为 0。指纹忽略 contig 次序和正反链、保留 contig 边界，非 ACGT 标准化为 N。仍有 **104 个 species、60 个 genera** 与旧数据的分类学范围重叠；本实验验证的是独立 assemblies/BioSamples，不宣称所有分类群均未见过，也不排除不同 BioSample 的近克隆。

44 个 species taxid 的 scientific name 含 sp./bacterium 等未命名或占位形式。因此“180 species”是按 NCBI species taxid 统计，不等于 180 个已正式命名物种。来源 snapshot 的 SHA256 和文件时间见 [data_provenance.json](benchmarks/results/robustness_20261001/data_provenance.json)；文件时间不冒充下载或发布日期。

首先计算同属候选及远缘 controls 的 FastANI，再在不读取任何 sketch 误差的条件下选择最终 pairs：每物种 4 对覆盖样本，各主要 ANI bin 最多补入 600 对，保留可用跨属/科 controls。它是有意覆盖场景的 benchmark，整体 MAE 不是 RefSeq 自然分布下的总体误差估计。

如果只与先前 **accuracy panel** 比较，而不计入 10,000-genome 性能集，本轮只有 3 个 species、18 个属曾出现过；即新增 177 个 species taxid 和 42 个属。因而它扩大了 accuracy 的分类学覆盖，但不是严格只由未见属构成的 leave-genus-out 测试。CSV 中另外保留“两个端点均来自未见 accuracy 属”的描述性分层，不能把这个补充分层的显著性算作预先指定的独立主结果。

最终远缘 controls 是 6 对同科跨属、120 对跨科。前者低于计划最多 120 对，因为本轮 60 属分散在 55 科，仅取每属一个 control representative；这个缺口没有通过改变抽样策略补齐。三个 ANI 分层与物种覆盖面板取并集，所以最终某个 bin 可能超过额外抽入的 600 对。

| 双向 FastANI 区间 (%) | 候选 qualified pairs | 最终 pairs |
| --- | --- | --- |
| 90-95 | 629 | 600 |
| 95-99 | 1032 | 624 |
| 99-100 | 574 | 574 |
| below90 | — | 44 |
| missing | — | 125 |

FastANI 要求两个方向都输出且 aligned fraction 均 ≥0.2，使用两方向均值。没有输出的 pair 保持缺失；仍计算其 exact Jaccard，不将 ANI 设成 0。候选不足的 bin 不另挑数据去填出期望效果。完整分类组成见 [natural_composition.csv](benchmarks/results/robustness_20261001/natural_composition.csv)，每个 genome 的来源和 QC 见 [natural_taxonomy.csv](benchmarks/results/robustness_20261001/natural_taxonomy.csv)。

## 3. 固定配置与公平口径

主分析 k=24：MinHash K=128/512/2048；BinDash bins=128/512/2048/8192（16 bit）；SetSketch p=12；FracMinHash scaled=5000/10000；KSSD dr3/half_subk6 和 dr4/half_subk7。次分析 k=20：MH512、BD2048、SS12、Frac10000、KSSD3。除 KSSD42 外，种子为 42/101/202。两种 k 均计算完整 canonical 2-bit k-mer 集合。

2/8/32 KiB 是 A+B 两个 sketch 的 payload 总和。固定大小的 MH/BD 等字节对应为 128/512、512/2048、2048/8192。Frac/KSSD 只在全部测试 seeds 的实测 pair bytes 均不超预算时纳入，并按字节选可容纳的最大已测档位，完全不看该 pair 的误差。这个预算 envelope 不等价于部署一个固定参数；advisor 复核则严格固定参数。SetSketch 只验证 p=12，因此它在 32 KiB 预算下仍只使用 8 KiB/pair。

ANI 使用当前接口对应的转换 $g(J)=(2J/(1+J))^{1/k}$。所有 ANI 误差表以百分点（pp）表示；Jaccard 误差为 [0,1] 上绝对误差。P95 是误差分位数，95% CI 是另一个量。将误差分为 $|\hat J-J|$、$|g(\hat J)-g(J)|$ 和 $|g(\hat J)-ANI_{ref}|$，另列 $|g(J)-ANI_{ref}|$。后者是转换模型与参考定义之间的差异，不是 sketch 采样误差。

误差均值/P95 包括全部 seeds；配对 CI 先在 pair 内平均 seeds，再按真实数据的属、或受控数据的 reference 重采样 2,000 次。跨属边采用两端 cluster 权重，同属边只取一次 cluster 权重。区间是 conditional on 测过的 seeds、逐项且未校正多重比较；不能据此声称整个配置空间被统计淘汰。

## 4. 新真实 genomes：逐算法的数据优势

k=24，8 KiB/pair。变长 sketch 的有效 pair 覆盖必须与精度一起读：

| FastANI (%) | 算法 | pairs | J MAE | 属内宏平均 J MAE | ANI sampling MAE (pp) | 总 ANI MAE vs FastANI (pp) | 总 ANI P95 (pp) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 90-95 | MinHash | 600 | 0.01211 | 0.01238 | 0.3420 | 1.4163 | 2.6635 |
| 90-95 | BinDash | 600 | 0.00555 | 0.00575 | 0.1519 | 1.4728 | 2.5826 |
| 90-95 | SetSketch-fixed | 600 | 0.00553 | 0.00577 | 0.1569 | 1.4486 | 2.5895 |
| 90-95 | FracMinHash | 520 | 0.01080 | 0.01072 | 0.3044 | 1.5062 | 2.6872 |
| 90-95 | KSSD | 600 | 0.02309 | 0.02084 | 0.8400 | 1.5626 | 3.0474 |
| 95-99 | MinHash | 624 | 0.01684 | 0.01709 | 0.1150 | 0.3957 | 0.9030 |
| 95-99 | BinDash | 624 | 0.00847 | 0.00852 | 0.0580 | 0.3801 | 0.8473 |
| 95-99 | SetSketch-fixed | 624 | 0.00694 | 0.00694 | 0.0475 | 0.3844 | 0.8551 |
| 95-99 | FracMinHash | 481 | 0.01734 | 0.01773 | 0.1105 | 0.4039 | 0.9134 |
| 95-99 | KSSD | 624 | 0.04154 | 0.04098 | 0.2858 | 0.4440 | 1.0404 |
| 99-100 | MinHash | 574 | 0.00981 | 0.00949 | 0.0294 | 0.0831 | 0.2781 |
| 99-100 | BinDash | 574 | 0.00470 | 0.00442 | 0.0142 | 0.0796 | 0.2546 |
| 99-100 | SetSketch-fixed | 574 | 0.00389 | 0.00388 | 0.0117 | 0.0798 | 0.2541 |
| 99-100 | FracMinHash | 449 | 0.01033 | 0.00929 | 0.0297 | 0.0831 | 0.3005 |
| 99-100 | KSSD | 574 | 0.02381 | 0.02336 | 0.0729 | 0.1098 | 0.3709 |

原始分层统计：[budget_accuracy.csv](benchmarks/results/robustness_20261001/budget_accuracy.csv)。属内宏平均先求每属均值再平均，减少拥有更多 selected pairs 的属对结论的影响；远缘跨属 pairs 不纳入这个宏平均，但保留在 pair-weighted 结果中。

| KiB/pair | 范围 | MH J MAE | BD J MAE | BD/MH | SS J MAE | SS/BD |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | all | 0.02388 | 0.01207 | 0.506 | — | — |
| 2 | J>0.5 | 0.02282 | 0.01181 | 0.518 | — | — |
| 8 | all | 0.01203 | 0.00582 | 0.484 | 0.00519 | 0.892 |
| 8 | J>0.5 | 0.01169 | 0.00576 | 0.492 | 0.00470 | 0.816 |
| 32 | all | 0.00599 | 0.00303 | 0.506 | 0.00519 | 1.712 |
| 32 | J>0.5 | 0.00591 | 0.00293 | 0.496 | 0.00470 | 1.604 |

“all” 的 Jaccard 统计包含全部 1,967 对，包括 FastANI 未输出的远缘 controls。CSV 中的 overall sampling error 也能在这些 pair 上计算，但 overall 总 ANI 误差只对有 reference 的 pairs 定义；两列不能在不同有效样本集合下直接相减。主要 ANI 比较以上面三个 qualified ANI bins 为准。

| 范围 | 配对差 | 误差指标 | pairs | clusters | 平均差 | 95% CI | 解释 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 90-95 | BinDash − MinHash | j_error | 600 | 30 | -0.00656 | [-0.00835, -0.00511] | 左侧误差更低，CI 不跨 0 |
| 90-95 | BinDash − MinHash | sampling_error | 600 | 30 | -0.19009 | [-0.24968, -0.14631] | 左侧误差更低，CI 不跨 0 |
| 90-95 | BinDash − MinHash | fastani_error | 600 | 30 | 0.05646 | [-0.02122, 0.13434] | CI 跨 0，方向未分离 |
| 95-99 | BinDash − MinHash | j_error | 624 | 56 | -0.00837 | [-0.00941, -0.00738] | 左侧误差更低，CI 不跨 0 |
| 95-99 | BinDash − MinHash | sampling_error | 624 | 56 | -0.05692 | [-0.06511, -0.04931] | 左侧误差更低，CI 不跨 0 |
| 95-99 | BinDash − MinHash | fastani_error | 624 | 56 | -0.01560 | [-0.02935, -0.00097] | 左侧误差更低，CI 不跨 0 |
| 99-100 | BinDash − MinHash | j_error | 574 | 57 | -0.00511 | [-0.00589, -0.00432] | 左侧误差更低，CI 不跨 0 |
| 99-100 | BinDash − MinHash | sampling_error | 574 | 57 | -0.01521 | [-0.01788, -0.01263] | 左侧误差更低，CI 不跨 0 |
| 99-100 | BinDash − MinHash | fastani_error | 574 | 57 | -0.00352 | [-0.00810, 0.00151] | CI 跨 0，方向未分离 |
| 99-100 | SetSketch-fixed − BinDash | j_error | 574 | 57 | -0.00081 | [-0.00121, -0.00032] | 左侧误差更低，CI 不跨 0 |
| 99-100 | SetSketch-fixed − BinDash | sampling_error | 574 | 57 | -0.00248 | [-0.00380, -0.00091] | 左侧误差更低，CI 不跨 0 |
| 99-100 | SetSketch-fixed − BinDash | fastani_error | 574 | 57 | 0.00013 | [-0.00175, 0.00199] | CI 跨 0，方向未分离 |
| J>0.5 | SetSketch-fixed − BinDash | j_error | 809 | 60 | -0.00106 | [-0.00141, -0.00071] | 左侧误差更低，CI 不跨 0 |
| J>0.5 | SetSketch-fixed − BinDash | sampling_error | 809 | 60 | -0.00391 | [-0.00524, -0.00261] | 左侧误差更低，CI 不跨 0 |
| J>0.5 | SetSketch-fixed − BinDash | fastani_error | 809 | 60 | -0.00054 | [-0.00294, 0.00173] | CI 跨 0，方向未分离 |
| all | BinDash − MinHash | j_error | 1967 | 60 | -0.00621 | [-0.00677, -0.00573] | 左侧误差更低，CI 不跨 0 |
| all | BinDash − MinHash | sampling_error | 1967 | 60 | -0.48845 | [-0.79863, -0.24415] | 左侧误差更低，CI 不跨 0 |
| all | BinDash − MinHash | fastani_error | 1842 | 60 | 0.01482 | [-0.01403, 0.04460] | CI 跨 0，方向未分离 |

误差差值为左减右，负值有利于左侧。J 指标仍为 fraction，ANI 指标在表内换算为 pp。这个表比“总平均第一名”更适合决定论文应写到什么强度。完整 k20/24、2/8/32 KiB 区间见 [paired_cluster_intervals.csv](benchmarks/results/robustness_20261001/paired_cluster_intervals.csv)。

MinHash 的存在理由仍应由满足精度后的构建成本支撑；不能把它称为每字节精度的默认最优。BinDash 的推荐需要同时得到等字节 J 和 sampling error 的支持，不能用总 ANI 误差中的偶然模型抵消代替统计精度。SetSketch 的 high-overlap 优势须在相应分层读取，不能升级成所有真实 ANI 最准。FracMinHash 与 KSSD 的变长 payload 需要先检查覆盖率：缺少可行配置不是一个应被删掉的失败样本。

| 算法 | 空间可行 pairs | 所有 pairs | 预算覆盖率 | 可行输出中的 empty 次数 |
| --- | --- | --- | --- | --- |
| MinHash | 1967 | 1967 | 100.00% | 0 |
| BinDash | 1967 | 1967 | 100.00% | 0 |
| SetSketch-fixed | 1967 | 1967 | 100.00% | 0 |
| FracMinHash | 1581 | 1967 | 80.38% | 0 |
| KSSD | 1967 | 1967 | 100.00% | 0 |

这一轮 FracMinHash 只保留 scaled=5000/10000，因此 8 KiB 下的约 80% 覆盖率是**冻结配置集合的覆盖率**，不是 FracMinHash 本身无法支持其余 genomes。KSSD 的 reduction 档位也较粗，8 KiB envelope 中常只使用远少于 8 KiB 的 payload；不能根据它在这个表中较大的误差宣称严格等字节全局劣势。逐配置空间和误差必须一起读。

### 4.1 构建优势在新输入上是否成立

| 算法 | native 参数 | seed 批次 | 构建含 setup (s) | 构建＋准备 (s) | 构建 min–max (s) | 平均 B/object | pair elapsed 中位数 (μs) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MinHash | 512 | 3 | 0.1927 | 0.1930 | 0.1926–0.1931 | 4096.0 | 2.888 |
| BinDash | 2048 | 3 | 0.2859 | 0.2860 | 0.2859–0.2868 | 4096.0 | 1.173 |
| SetSketch-fixed | 12 | 3 | 0.2548 | 0.2553 | 0.2537–0.2549 | 4096.0 | 63.195 |
| FracMinHash | 5000 | 3 | 0.2257 | 0.2263 | 0.2256–0.2260 | 5992.3 | 3.299 |
| FracMinHash | 10000 | 3 | 0.2263 | 0.2275 | 0.2256–0.2271 | 2993.1 | 2.347 |
| KSSD | 3 | 1 | 1.2043 | 1.2045 | 1.2043–1.2043 | 7202.7 | 5.800 |
| KSSD | 4 | 1 | 21.8718 | 21.8719 | 21.8718–21.8718 | 230.2 | 0.611 |

这是整个 900-genome batch 的 wall time，包含 KSSD 的共享初始化；输入解析和 pair 查询不混入构建。三个 seed 批次不是大量独立性能重复，KSSD 只有一批。系统运行时波动与种子差异不能由这张表完全分离。完整原始批次和 input time 在 [build_batches.csv](benchmarks/results/robustness_20261001/build_batches.csv)；本轮每对查询计时在 [pair_query_cost.csv](benchmarks/results/robustness_20261001/pair_query_cost.csv)，后者仍是 64 线程负载中的 elapsed time，不能冒充独占 pair latency 或新的 all-vs-all crossover。

旧的 N*≈8593/2177 和实际 all-vs-all 证据保留原有机器及 workload 限定。本轮没有第二 CPU，论文不能增加 across-architectures 的结论。

### 4.2 k=20 的预定敏感性复核

| k | all MH J MAE | all BD J MAE | BD/MH | J>0.5 pairs | 高 J BD MAE | 高 J SS MAE | SS/BD |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | 0.01232 | 0.00606 | 0.492 | 863 | 0.00583 | 0.00512 | 0.879 |
| 24 | 0.01203 | 0.00582 | 0.484 | 809 | 0.00576 | 0.00470 | 0.816 |

两种 k 下 BinDash 的等字节优势和 SetSketch 的高 J 改善方向均保留，幅度不同。J>0.5 是各自 k 下的 truth 分层，所以两个 k 的高 J 行不是完全相同的 pair 集合；不能把差额直接解释成 k 的因果收益。

## 5. 混合 substitution/indel：什么保持稳定，什么会改变

20 个真实 5 Mb references 来自上一轮，只用于配对的 mutation-model sensitivity，不计为第二份独立生物验证。每个 reference 有 4 个目标 alignment-column identities × 4 种 affected-base composition，共 16 个变体。indel 的 insertion/deletion 碱基数相等，事件长度为截断 geometric(p=0.5)，上限 10；事件不重叠，不修改模糊碱基。composition 以受影响碱基数统计，不是 event 数。

每个变体保存可重放的 compressed edit script，并精确记录 M、S、I、D。目标定义为 M/(M+S+I+D)，同时保留 M/(M+S) 的 ungapped homologous identity。脚本按独立的字符串实现测试过重放，真实 320 个输出也逐一通过保存后重放、长度及计数验证。它们是已知构造真值；FastANI 是另外一列参考。

| column target | indel bases share | mean exact J | ungapped identity (%) | FastANI (%) | exact-J→column MAE (pp) | exact-J→ungapped MAE (pp) | exact-J→FastANI MAE (pp) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 90% | 0% | 0.04176 | 90.0000 | 88.8546 | 0.0199 | 0.0199 | 1.1641 |
| 90% | 10% | 0.04722 | 90.9091 | 89.3615 | 0.4611 | 0.4480 | 1.0995 |
| 90% | 25% | 0.05704 | 92.3077 | 90.1411 | 1.1406 | 1.1671 | 0.9995 |
| 90% | 50% | 0.07889 | 94.7368 | 91.4634 | 2.3022 | 2.4346 | 0.8389 |
| 95% | 0% | 0.17153 | 95.0000 | 94.4042 | 0.0126 | 0.0126 | 0.6072 |
| 95% | 10% | 0.18392 | 95.4774 | 94.6682 | 0.2461 | 0.2313 | 0.5779 |
| 95% | 25% | 0.20502 | 96.2025 | 95.0640 | 0.6077 | 0.5949 | 0.5436 |
| 95% | 50% | 0.24733 | 97.4359 | 95.7405 | 1.2196 | 1.2163 | 0.4791 |
| 97% | 0% | 0.31714 | 97.0000 | 96.5842 | 0.0045 | 0.0045 | 0.4169 |
| 97% | 10% | 0.33285 | 97.2919 | 96.7411 | 0.1488 | 0.1431 | 0.4077 |
| 97% | 25% | 0.35825 | 97.7330 | 96.9871 | 0.3703 | 0.3627 | 0.3832 |
| 97% | 50% | 0.40722 | 98.4772 | 97.4053 | 0.7471 | 0.7301 | 0.3418 |
| 99% | 0% | 0.64646 | 99.0000 | 98.8020 | 0.0023 | 0.0023 | 0.1959 |
| 99% | 10% | 0.65981 | 99.0991 | 98.8605 | 0.0489 | 0.0502 | 0.1884 |
| 99% | 25% | 0.68047 | 99.2481 | 98.9490 | 0.1251 | 0.1230 | 0.1761 |
| 99% | 50% | 0.71698 | 99.4975 | 99.0963 | 0.2523 | 0.2452 | 0.1560 |

这个表中的差异不随 sketch 分辨率自动消失。对于 substitution、indel 事件和同源/非同源内容，Jaccard→ANI 的映射前提不同。增加 sketch bytes 主要减少 sampling error，不能把它当作修复任意 ANI 参考偏差的方法。

k=24、8 KiB/pair；下表按同一 reference/target 直接比较。除最后一列为 Jaccard MAE 比值外，算法列均为 ANI sampling MAE（pp）：

| column target (%) | indel (%) | MinHash | BinDash | SetSketch-fixed | FracMinHash | KSSD | BD/MH J MAE |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 90 | 0 | 0.5723 | 0.3450 | 0.3823 | 0.4402 | 0.8339 | 0.590 |
| 90 | 10 | 0.5530 | 0.2837 | 0.3113 | 0.3881 | 1.2156 | 0.519 |
| 90 | 25 | 0.5377 | 0.2789 | 0.3008 | 0.3831 | 1.2057 | 0.528 |
| 90 | 50 | 0.4577 | 0.2093 | 0.2191 | 0.3424 | 1.0536 | 0.463 |
| 95 | 0 | 0.2691 | 0.1347 | 0.1170 | 0.2090 | 0.4187 | 0.504 |
| 95 | 10 | 0.2177 | 0.1124 | 0.1193 | 0.1967 | 0.4981 | 0.513 |
| 95 | 25 | 0.2453 | 0.1104 | 0.0893 | 0.1720 | 0.3755 | 0.448 |
| 95 | 50 | 0.2305 | 0.0949 | 0.0893 | 0.1728 | 0.3157 | 0.419 |
| 97 | 0 | 0.1575 | 0.0769 | 0.0645 | 0.1316 | 0.3358 | 0.487 |
| 97 | 10 | 0.1496 | 0.0699 | 0.0647 | 0.1247 | 0.2751 | 0.476 |
| 97 | 25 | 0.1325 | 0.0744 | 0.0592 | 0.1157 | 0.3074 | 0.560 |
| 97 | 50 | 0.1367 | 0.0554 | 0.0485 | 0.1027 | 0.2368 | 0.402 |
| 99 | 0 | 0.0555 | 0.0337 | 0.0247 | 0.0520 | 0.2237 | 0.613 |
| 99 | 10 | 0.0569 | 0.0320 | 0.0254 | 0.0578 | 0.1486 | 0.564 |
| 99 | 25 | 0.0645 | 0.0315 | 0.0254 | 0.0695 | 0.1320 | 0.494 |
| 99 | 50 | 0.0617 | 0.0306 | 0.0223 | 0.0604 | 0.1114 | 0.499 |

| BD − MH 指标 | 已测 target×composition cells | 平均差有利 BD | 逐项 CI 有利 BD | 逐项 CI 有利 MH | CI 跨 0 |
| --- | --- | --- | --- | --- | --- |
| j_error | 16 | 16 | 16 | 0 | 0 |
| sampling_error | 16 | 16 | 16 | 0 | 0 |
| column_error | 16 | 11 | 6 | 0 | 10 |
| homologous_error | 16 | 13 | 6 | 0 | 10 |
| fastani_error | 16 | 6 | 1 | 0 | 15 |

这些是冻结的 16 个 cells 的方向核查，不能把 16 个相关 cells 当 16 个独立证据倍增信心。支持“robust”的正文应指明是 Jaccard / sampling 精度，还是对某个 ANI 参考的总误差；若总误差方向逆转，保留该结果。补充按 observed exact J 分层的结果在 [mixed_observed_j_accuracy.csv](benchmarks/results/robustness_20261001/mixed_observed_j_accuracy.csv)，用于分辨 mutation composition 改变了 J 分布之后的效应，不作为匹配了所有变量的因果实验。

## 6. NUCmer 参考检查：避免由单一 ANI reference 决定胜负

| 面板 | 预选 pairs | 命令成功 | 1-to-1 AF 双向≥0.2 | FastANI/NUCmer 平均绝对差 (pp) | 最大绝对差 (pp) |
| --- | --- | --- | --- | --- | --- |
| natural | 60 | 60 | 60 | 0.2280 | 1.0623 |
| mixed | 80 | 80 | 80 | 0.3950 | 1.2062 |

同一审核子集、k=24、8 KiB/pair，分别更换 ANI reference（均为 MAE，pp）：

| FastANI 分层 (%) | 算法 | 审核 pairs | 对 FastANI | 对 NUCmer |
| --- | --- | --- | --- | --- |
| 90-95 | MinHash | 20 | 1.4477 | 1.2466 |
| 90-95 | BinDash | 20 | 1.5493 | 1.3415 |
| 90-95 | SetSketch-fixed | 20 | 1.5184 | 1.3099 |
| 95-99 | MinHash | 20 | 0.3257 | 0.2599 |
| 95-99 | BinDash | 20 | 0.3064 | 0.2418 |
| 95-99 | SetSketch-fixed | 20 | 0.3139 | 0.2481 |
| 99-100 | MinHash | 20 | 0.0699 | 0.0765 |
| 99-100 | BinDash | 20 | 0.0630 | 0.0737 |
| 99-100 | SetSketch-fixed | 20 | 0.0628 | 0.0689 |

自然数据每 ANI bin 20 对，混合数据为预先按哈希选出的 5 个 references ×16 条件；不按 sketch 误差挑审核 pairs。NUCmer/dnadiff 使用 1-to-1 identity 及两个方向的对齐覆盖率。NUCmer 与 FastANI 衡量/筛选的对齐内容不完全相同，彼此都不是适用于所有场景的绝对真理。

逐算法在两个参考下的 MAE/P95 见 [alignment_sketch_accuracy.csv](benchmarks/results/robustness_20261001/alignment_sketch_accuracy.csv)，exact-J 模型与各参考的差异见 [alignment_reference_comparison.csv](benchmarks/results/robustness_20261001/alignment_reference_comparison.csv)。正文要把参考敏感性作为限制保留；不能只引用对偏好算法有利的参考。

## 7. 冻结 advisor 的外部检验

| 外部面板 | indel share | 固定 requests | 原 calibration 无候选 | 本面板无对应 pairs | 通过 | 失败 | 失败且 exact-J model 也超阈值 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| mixed | 0% | 48 | 21 | 0 | 27 | 0 | 0 |
| mixed | 10% | 48 | 21 | 0 | 22 | 5 | 3 |
| mixed | 25% | 48 | 21 | 0 | 10 | 17 | 15 |
| mixed | 50% | 48 | 21 | 0 | 6 | 21 | 21 |
| natural | — | 48 | 30 | 0 | 6 | 12 | 6 |

合并重复 workload 请求后，以下是真实数据的具体通过/失败证据（均为 pp）：

| ANI (%) | KiB/pair | MAE/P95 限值 | 冻结配置 | 新 MAE | 新 P95 | exact-J model P95 | 结果 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 95-99 | 8 | 0.5 / 1.0 | MinHash 128 | 0.4290 | 1.0612 | 0.8122 | fail |
| 95-99 | 8 | 0.5 / 1.0 | BinDash 128 | 0.4288 | 1.0767 | 0.8122 | fail |
| 95-99 | 32 | 0.5 / 1.0 | MinHash 128 | 0.4290 | 1.0612 | 0.8122 | fail |
| 95-99 | 32 | 0.5 / 1.0 | BinDash 128 | 0.4288 | 1.0767 | 0.8122 | fail |
| 99-100 | 8 | 0.1 / 0.2 | SetSketch-fixed 12 | 0.0798 | 0.2541 | 0.2503 | fail |
| 99-100 | 8 | 0.1 / 0.2 | BinDash 512 | 0.0833 | 0.2724 | 0.2503 | fail |
| 99-100 | 8 | 0.5 / 1.0 | MinHash 128 | 0.0983 | 0.3449 | 0.2503 | pass |
| 99-100 | 8 | 0.5 / 1.0 | BinDash 128 | 0.0998 | 0.3286 | 0.2503 | pass |
| 99-100 | 32 | 0.1 / 0.2 | MinHash 2048 | 0.0793 | 0.2517 | 0.2503 | fail |
| 99-100 | 32 | 0.1 / 0.2 | BinDash 512 | 0.0833 | 0.2724 | 0.2503 | fail |
| 99-100 | 32 | 0.5 / 1.0 | MinHash 128 | 0.0983 | 0.3449 | 0.2503 | pass |
| 99-100 | 32 | 0.5 / 1.0 | BinDash 128 | 0.0998 | 0.3286 | 0.2503 | pass |

直接读取上一轮 `heldout_selection.csv`，不重新做 calibration、不训练 ML、不用本轮错误选择备用配置。Natural 使用原 Mash-family 请求；controlled 使用原 substitution 请求，并对四种 mutation composition 分别复核。Mixed 的区间归属用 nominal column target，防止整数碱基舍入让刚好 95% 的 pair 漂进相邻 bin；误差仍对 achieved column identity 计算。Natural 对 FastANI 的真实 bin 归属计算。

通过条件仍是 MAE≤ε、P95≤2ε、全部观测 payload≤预算。只有同时有冻结候选和对应外部 pairs 的请求才进入通过/失败分母。原 calibration 无候选的请求保留 `no_frozen_candidate`，不算成功；若已有候选但外部缺少对应 pairs，则标记 `no_external_pairs`。成本模型仍来自上一轮 10,000-genome CPU 测量，本轮不更新其 crossover。

最后一列的 exact-J residual 不是一个严格的误差下界：有限 sketch 的 bias 可能碰巧抵消 reference/model residual。但若完整集合的转换都超过约束，它已足以说明“只增加 sketch size”没有解决问题的保证。新数据失败必须限制 advisor 的适用域。实用接口应在超出已验证数据域或参考定义时显示证据不足，而非承诺 ε。所有逐请求选择、原 holdout 状态及新失败见 [frozen_advisor_external.csv](benchmarks/results/robustness_20261001/frozen_advisor_external.csv)。

## 8. 每个算法的最终推荐应怎样写

| 算法 | 保留的使用理由 | 本轮补充证据如何限制表述 |
| --- | --- | --- |
| MinHash/FastKMV | accuracy 已满足后的构建型 workload 候选 | 以第 4 节新输入的同预算 build time 复核；不能等同于最佳 accuracy/byte；ANI reference/model 残差需要单独达标 |
| BinDash | 普通集合的紧凑统计估计；旧实测中 query reuse 高时有成本优势 | 用新真实与 mixed 的等字节 J、sampling error 和 paired CI 支撑；性能仍限同机；真实 ANI 总误差另验 |
| SetSketch-fixed | 愿意支付查询成本以换取部分 high-overlap joint-estimation 精度 | 只在第 4 节支持的分层推荐；不能从 J>0.5 自动推断任意 99–100% ANI 都更好；实现版本必须标明 |
| FracMinHash | directed containment 且 query 有足够保留样本 | 本轮普通 similarity 是外部检验，不替代旧 containment 证据；变长空间与缺失可行档位必须报告 |
| HLL | cardinality-only 的轻量起点 | 本轮未重测，继续依赖上一阶段的 cardinality evidence，不将新 Jaccard 结论套给 HLL |
| weighted ProbMinHash | probability Jaccard 的 abundance comparison | 本轮未重测；CAMI profile 层验证与完整 reads→k-mer workflow 的边界保持不变；不是 min/max weighted Jaccard |
| KSSD | 提供 reduction backend；允许局部低空间/精度竞争 | 不因必须分配 niche 而强行设默认入口；本轮变长 payload、单 reduction seed 和含 setup 成本共同解读；未证明全局无用 |
| unweighted ProbMinHash / OrderMinHash | 前者降低普通集合默认优先级；后者排除当前 scope | 本轮没有新增这两者的证据，不宣称本轮完成它们的全面淘汰 |

论文的三条主贡献可以继续围绕 unified characterization、条件化 selection regimes、empirical constrained advisor。Advisor 的贡献应包括可明确返回无合格候选或外部验证失败，而不是保证任意 genome 都能达到输入的误差阈值。算法精度的方向性复现，与固定阈值和 crossover 的普适性，是需要分开报告的结果。

## 9. 实验冻结与交付清单

本轮不再增加普通 Jaccard 参数、size ratio、weighted registers 或性能平台。当前证据仍不覆盖跨架构、reads 全流程、长 indel、rearrangement 或普通 similarity 的大型真核 genomes；这些是论文适用边界，不是本轮继续扩实验的理由。

原始输入与编辑记录：`/home/bigssd/zt/rabbitsketch2_benchmark/data/robustness_20261001`。原始命令、stdout/stderr、逐批 SHA256、代码与二进制快照：`/home/bigssd/zt/rabbitsketch2_benchmark/runs/robustness_20261001`。分析表、PNG/PDF 与审计：`/home/bigssd/zt/rabbitsketch2_benchmark/analysis/robustness_20261001`。所有旧阶段文件保留。

- [raw_audit.json](benchmarks/results/robustness_20261001/raw_audit.json)：94 批次完整性、raw SHA256、exact/native 输入摘要一致性。
- [data_provenance.json](benchmarks/results/robustness_20261001/data_provenance.json)：排重、分类组成、来源 snapshot 及限制。
- [configuration_accuracy.csv](benchmarks/results/robustness_20261001/configuration_accuracy.csv)、[budget_accuracy.csv](benchmarks/results/robustness_20261001/budget_accuracy.csv)、[budget_coverage.csv](benchmarks/results/robustness_20261001/budget_coverage.csv)：固定配置与预算比较。
- [paired_cluster_intervals.csv](benchmarks/results/robustness_20261001/paired_cluster_intervals.csv)：全部配对差与区间。
- [mixed_model_summary.csv](benchmarks/results/robustness_20261001/mixed_model_summary.csv)：编辑真值、FastANI 和 exact-J 模型差异。
- [frozen_advisor_external.csv](benchmarks/results/robustness_20261001/frozen_advisor_external.csv)：成功、失败及无证据请求。
- [independent_genome_accuracy.png](benchmarks/results/robustness_20261001/independent_genome_accuracy.png) / [independent_genome_accuracy.pdf](benchmarks/results/robustness_20261001/independent_genome_accuracy.pdf)。
- [mixed_mutation_sampling.png](benchmarks/results/robustness_20261001/mixed_mutation_sampling.png) / [mixed_mutation_sampling.pdf](benchmarks/results/robustness_20261001/mixed_mutation_sampling.pdf)。
- [mixed_mutation_model_residual.png](benchmarks/results/robustness_20261001/mixed_mutation_model_residual.png) / [mixed_mutation_model_residual.pdf](benchmarks/results/robustness_20261001/mixed_mutation_model_residual.pdf)。

本报告由 [write_report.py](benchmarks/robustness/write_report.py) 从冻结 CSV 生成表格；主要分析为 [analyze.py](benchmarks/robustness/analyze.py)。时间值与统计误差分别按以上口径读取。最终 artifact hashes 见 [completion.json](benchmarks/results/robustness_20261001/completion.json)。
