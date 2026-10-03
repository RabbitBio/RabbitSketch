**RabbitSketch2 算法选择的详细数据分析：每个算法为什么值得推荐，以及推荐到什么范围。**

分析日期：2026-10-01。本文主体保留普通 Jaccard 初筛和第二阶段任务评测的冻结结果，面向 RabbitSketch2 的 selection guide 与论文 Results/Discussion；表格直接从源 CSV 提取，数值按展示精度四舍五入。

**后续验证已完成，最新结论见同目录的 [独立多物种与混合突变验证](RabbitSketch2_robustness_validation.md)。** 新增 900 genomes、180 species taxid、60 属的 1,967 对真实比较，以及 20 references 的 320 对 substitution/indel 扰动；94 个 exact/native 批次和 140 对 NUCmer 审核全部完成。未增加第二 CPU 平台，OrderMinHash 仍排除。

新数据复现了主要选择方向：k=24、8 KiB/pair 下，BinDash 的 Jaccard MAE 为 MinHash 的 0.484 倍；SetSketch-fixed 在 J>0.5 的 MAE 为 BinDash 的 0.816 倍；900-genome 构建加准备时间为 MinHash512 0.1930 s、BinDash2048 0.2860 s。混合突变的 16 种 target×composition 条件中，BinDash 的等字节 Jaccard 与 ANI sampling error 均优于 MinHash，逐项 reference-cluster 区间均不跨零。

**更需要收紧的是 advisor 的阈值迁移，而不是这些采样优势。** 原 Mash-family 的 18 个冻结选择只有 6 个通过新真实数据的误差约束；短 indel affected-base 占比为 0/10/25/50% 时，原 substitution-family 的 27 个冻结选择分别通过 27/22/10/6 个。失败全部保留，没有用外部数据换候选。99–100% ANI 的新真实数据中，完整 exact J 转换本身的 P95 约 0.2503 pp，已超过原 0.2 pp 约束。论文应将“优势方向复现”与“特定数据域内的误差承诺”分开。下面主体中的旧样本量、41/45 holdout 结果及 substitution-only 限制均指第一、二阶段；新增证据和修订后的适用边界以链接报告为准。

当前 MinHash 指 FastKMV，只计一个候选。KSSD 已恢复，OrderMinHash 按用户要求排除。SetSketch 的结果均指隔离工作副本中的 partitioned-joint 修正版；生产目录尚未被该实验补丁替换。第一阶段旧报告中“KSSD 不再参与”的范围声明已被第二阶段恢复实验取代。

本文遵循一个判定顺序：先确认算法估计的是用户需要的指标，再检查误差、有效输出比例和空间限制，最后比较满足要求的配置的实际成本。某种算法在一个指标上有优势，不自动意味着它适合所有 genome tasks；某种算法没有成为默认推荐，也不等于已证明它在所有条件下被支配。

**1. 当前推荐由下面这些数据支撑。** 表中的“强证据”只指本次已测数据和实现中有直接、反复的观测支持，不代表跨物种、跨机器或所有参数的数学保证。

| 算法或模式 | 推荐场景 | 最直接的数据优势 | 推荐边界 |
|---|---|---|---|
| MinHash / FastKMV | 精度达标、构建为主、少量查询 | 每对象 4 KiB 时，10,000 genomes 构建及准备 2.390 s，低于 BinDash 的 3.533 s；较小规模的实际全矩阵总成本也更低 | 等字节精度通常弱于 BinDash；低 J、小 K 的零命中会损害 ANI |
| BinDash | 普通集合相似度、有限字节预算、重复查询 | 合成 Jaccard 的等预算 MAE 比值 BD/MH 中位数约 0.48–0.51；4 KiB/object 的 10,000-genome 查询约比 MinHash 快 4.6 倍 | 16-bit 配置下的实测；极低 J 仍需检查相对误差；原生 occupancy cardinality 饱和后不能支持可靠 containment |
| SetSketch-fixed | 联合估计精度重要，查询较少，或精度约束筛掉较快配置 | 真实 J>0.5、8 KiB/pair 时，Jaccard MAE 比 BinDash 低约 24%–31%；部分独立验证的约束请求选择它 | 查询昂贵；结果依赖模型修正；高重叠采样优势不能直接升级成真实 ANI 最优性 |
| HLL | Cardinality-only 的优先起点，尤其小输入和中等精度要求 | p=12、4 KiB/object 时，10 kb–50 Mb 的相对 MAE 约 1.0%–1.2%；小输入构造与计数更轻 | 高精度、大 register 配置中 SetSketch 也可更准；不是所有 cardinality 约束下的唯一赢家 |
| FracMinHash | Directed containment，query 有足够保留样本，重视查询成本 | C≈0.5、ratio=50、8 KiB/pair 时，MAE 0.0968，低于 MinHash 的 0.1561；非空的完全包含 query 可返回 1 | 同条件 SetSketch 更准；ratio=500 时空样本与大误差明显；不能仅按大小比例推荐 |
| ProbMinHash-weighted | 相对 abundance/composition 的 probability Jaccard | CAMI profile 的 m=2048：MAE 0.00774、P95 0.01871；增大 m 后误差持续降低 | 不估计 min/max weighted Jaccard；不表示绝对总量差异；真实验证是 profile 层面 |
| ProbMinHash-unweighted | 普通集合任务中降低默认优先级 | 未观察到稳定的等字节精度收益；其专门价值来自 weighted 模式 | 降低优先级不是正式全局淘汰；不影响 weighted 模式的独立定位 |
| KSSD | 保留 backend 和完整结果，暂不设默认入口 | dr=4 平均约 264 B/object；部分中等 ANI 区间 MAE 优于相近字节 BinDash | 初始化、构建、查询成本较高，且没有形成独立默认选择区；不能宣称所有 KSSD 配置无价值 |

支持这些判断的两个主报告是[第一阶段系统分析][stage1_report]和[第二阶段任务分析][stage2_report]。下面逐项展开证据、机制和限制。

**2. 公平比较首先需要统一预算、误差和时间的定义。** 若忽略这些定义，表格中的数字很容易被解读成错误的算法结论。

本报告中，2/8/32 KiB 的 accuracy budget 均为 A+B 两个 sketch 的查询 payload 总和。对于对称固定大小 sketch，分别等于每对象 1/4/16 KiB；cardinality 表单独标明每对象空间。KiB=1024 bytes。

预算精度表按实际字节选择可容纳的已测档位，选择时不使用该 pair 的误差；变长的 FracMinHash/KSSD 可能在不同 pair 上使用不同档位，因此预算表不能被解释为一个固定 native 参数在所有输入上的表现。后面的 advisor 则选择固定配置，以 calibration 中观察到的最大 pair bytes 检查预算，再在 holdout 检查是否仍符合。

| 8 KiB/pair 下的典型配置 | 每对象保存的信息 | 每对象 payload |
|---|---|---:|
| MinHash K=512 | 512 个 64-bit hashes | 4 KiB |
| BinDash bins=2048 | 2048 个 16-bit signatures | 4 KiB |
| HLL p=12 | 4096 个 8-bit registers | 4 KiB |
| SetSketch p=12 | 4096 个 8-bit registers | 4 KiB |
| ProbMinHash m=512 | 512 个 64-bit winners，用于相似度查询 | 4 KiB |
| FracMinHash | 由集合大小与 scaled 决定的 hashes | 按实测字节限制选档，不固定 |
| KSSD | reduction 后的 32/64-bit tuples | 按实测字节限制选档，不固定 |

因此，“K=512、bins=512、m=512”不属于等资源比较。这里的 payload 也不是完整序列化文件大小、整个对象的 resident memory 或构建峰值内存。ProbMinHash 的 winner-only 空间针对查询；若还需要继续更新，不能忽略所需的额外状态。大型运行的 RSS 包含全部驻留 FASTA、共享字典和临时分配，不能拿来代替 sketch payload。

Jaccard、containment、probability Jaccard 的误差默认以 [0,1] 的 fraction 表示。ANI 表使用“百分点”：0.001 fraction=0.1 个百分点。Cardinality 表使用相对误差百分比。P95 是绝对误差分布的第 95 百分位，既不是 95% 置信区间，也不是对未来输入的保证。

性能表区分两类测量。第二阶段 all-vs-all 是真实不同对象间的完整矩阵，查询时间取三遍的中位数；构建时间包括需要的共享 setup、构建和查询缓存准备，排除 FASTA/gzip 解析与完整结果写盘。第一阶段 per-pair 计时是在 64 线程负载中的每对 elapsed time，不是独占单核延迟；本文不会把它直接外推成大矩阵性能。

本轮有独立 calibration/validation 拆分：受控数据按 reference 分成 10/10，Mash subset 按 genome 分成 71/71。不同 sketch seeds 是同一 biological pair 的重复测量，不是额外独立物种。第二阶段 ANI 区间按 genome/reference 重采样，第一阶段 synthetic 区间以独立随机数据重复为单位；均为逐项区间，未作多重比较校正。KSSD 当前 seed 字段不改变 reduction sampling，因此只测一个 seed，不把重复运行当独立统计样本。

**3. 这批数据覆盖了算法推荐需要的不同任务。** 不同任务使用各自的 truth，避免拿普通 Jaccard 代替所有需求。

| 实验 | 数据规模与设计 | Truth 或比较参考 | 主要回答的问题 |
|---|---|---|---|
| 普通 Jaccard 初筛 | 84 个可行 synthetic regimes，每场景 20 次数据重复；修正及 FracMinHash 补档后 86,240 条测量 | 已知集合大小与交集；保持 J≤1/size ratio 的可行性 | 等字节采样精度、相似度依赖 |
| 真实集合 Jaccard | 200 对、396 个 microbial assemblies，k=20/24、三个 seeds，共 32,400 条 native 测量 | 完整 canonical 2-bit k-mer 集合 | 合成观察是否延续到真实序列 |
| 自然 genome ANI | Mash 原始 accession 的版本匹配 subset：142 genomes、10,011 pairs；另有 200 对多样性样本 | 双向 FastANI；缺失保持缺失 | Jaccard 排序是否能直接转为 ANI 排序 |
| 受控 ANI | 20 个不同属的真实 5 Mb references，14 档 substitution identity，共 280 pairs | 已知 positional identity，并计算 exact Jaccard | 分离采样误差与转换误差 |
| 独立 alignment 复核 | 固定分层抽取 74 对，71 对达到双侧覆盖条件 | NUCmer/dnadiff 1-to-1 identity，与 FastANI 分开报告 | 推荐是否对参考方法敏感 |
| Cardinality | 85 个对象：10 kb、100 kb、1 Mb、5 Mb，以及五个 50 Mb 窗口 | Exact distinct canonical k-mers | HLL 是否有独立保留价值 |
| Directed containment | 20 references ×7 query lengths ×5 overlap levels，共 700 pairs；最终 donor/reference 保持 split 隔离 | Exact intersection/query cardinality，双方向 | FracMinHash 与 register methods 的使用边界 |
| 实际 all-vs-all | N=256/1024/4096/10000；补充紧凑配置在 N=1024/10000 | 每个不同 accession 构建一次，实际遍历 i<j | 构建与查询优势何时反转 |
| Weighted similarity | 720 组 synthetic weights；CAMI Marine 10 profiles、45 pairs；20 sketch seeds | Exact probability Jaccard，与 min/max Jaccard 分开 | abundance 专用指标是否被正确估计 |

第二阶段最终分析纳入 434 个批次、1,290,828 条 native pair measurements，另有 5,100 条 native cardinality 对象测量和 42 个全矩阵配置/规模批次。共执行过 574 个批次；其中旧 containment 生成方式造成 donor 跨 split 的 140 个探索性批次不进入最终结论。原始输出保留，没有覆盖或隐去这些差异。62 项回归检查通过，运行及交付哈希见[完成记录][completion]。

Mash 这里不是完整重现原始 500-genome benchmark。归档中同一 pair ID 存在冲突记录，因此本实验没有任取一条作为真值，而是对匹配版本的 142-genome subset 重算参考。50 Mb 的五个对象是同一 human reference 的不同窗口，不是五个独立个体。受控 ANI 仅覆盖 substitutions；reads、indels、rearrangements 和跨架构表现未被本轮结果证明。

**4. MinHash/FastKMV 的推荐依据是：满足精度之后，它的构建优势可以降低总成本。** 它不是普通 Jaccard 等字节精度最好的算法，但“不是最精确”不等于“没有最佳使用场景”。

先看实际计算成本。以下均为 k=24、每对象 4 KiB，MinHash K=512 与 BinDash bins=2048。时间为秒，已排除输入解析与结果写盘。

| genomes | MH 构建/准备 | BD 构建/准备 | MH 查询 | BD 查询 | MH 合计 | BD 合计 |
| --- | --- | --- | --- | --- | --- | --- |
| 256 | 0.0731 | 0.1086 | 0.0011 | 0.0003 | 0.0742 | 0.1088 |
| 1,024 | 0.2564 | 0.3888 | 0.0174 | 0.0034 | 0.2738 | 0.3922 |
| 4,096 | 0.9887 | 1.4560 | 0.2802 | 0.0556 | 1.2689 | 1.5116 |
| 10,000 | 2.3897 | 3.5329 | 1.6960 | 0.3654 | 4.0856 | 3.8984 |

在 N=256、1024、4096 时，MinHash 的构建优势足以覆盖查询劣势，实测总成本更低；到 N=10000，BinDash 的总成本反超。固定这些配置，反转确实发生在测试规模范围内，而不只是来自同一 pair 重复缓存查询的外推。

这与当前实现吻合：FastKMV 保留最小 K 个完整 hash，第 K 小值作为阈值；大量后续 hash 只需经过过滤，不进入保留集合更新。当前实现对 rolling/hash/filter 路径做了优化。这个机制解释实测方向，但构建领先的证据来自真实 FASTA 测量，不能只凭算法名称推出。

第一阶段 k=24、8 KiB/pair、200 个真实 pairs 的平均绝对误差和构建中位数进一步说明“快”和“准”分别在哪里：

| 算法 | Jaccard MAE | 每对构建中位数（ms） | warm query 中位数（μs） |
| --- | --- | --- | --- |
| MinHash/FastKMV | 0.004473 | 26.469 | 0.645 |
| BinDash | 0.002538 | 39.176 | 0.368 |
| SetSketch-fixed | 0.002933 | 34.916 | 60.128 |
| HLL | 0.005383 | 32.853 | 12.022 |
| FracMinHash | 0.004388 | 30.407 | 0.521 |
| ProbMinHash-unweighted | 0.004479 | 51.042 | 0.479 |

MinHash 一对构建约 26.47 ms，BinDash 约 39.18 ms；但 Jaccard MAE 为 0.004473 对 0.002538。实际选择要先检查任务是否接受 MinHash 的误差。如果不能接受，构建较快也不能使该配置成为合格推荐。

低相似度是它最需要检查的边界。对于远大于 K 的集合，bottom-k union sample 中共同元素的期望数量近似为 KJ，零共同样本概率近似为 (1−J)^K。以 J=0.001、K=512 为例，期望共同样本约 0.512，零命中概率约 60%。这是用于解释有限样本的近似推导，不是把基因组真实 J 固定为某个理论值。

受控 substitution、k=24、8 KiB/pair、validation 的 80–90% ANI 区间中，结果如下。ANI 误差单位为百分点。

| 算法 | ANI MAE（pp） | ANI P95（pp） | 零 J 比例 |
| --- | --- | --- | --- |
| MinHash/FastKMV | 19.624 | 80.000 | 23.33% |
| BinDash | 1.007 | 2.457 | 0% |
| SetSketch-fixed | 12.040 | 80.000 | 13.33% |

MinHash 的零 J 比例为 23.3%，而当前 ANI 转换将零 J 映射为 ANI=0，导致非常大的尾部误差。此时“MinHash 构建快”已经不足以支持推荐；应增大 K、重新评估 k，或选择该预算下误差更低的候选，而不是忽略失败率。

因此，MinHash 的可执行建议是：普通 set/Jaccard 或经过验证的 ANI 任务，先在 calibration 上确定足够的 K，再用独立 validation 检查；若构建频繁、每个 sketch 参与比较次数较少，优先考察 FastKMV。若 sketch 主要用于长期反复查询，或者字节预算很紧，应同时比较 BinDash。不能据本结果证明所有 MinHash 变体都可被替代，也不能把当前 FastKMV 与 Legacy MinHash 重复计为两个独立推荐。

源数据：[实际矩阵成本][workload]、[真实 Jaccard 汇总][real_summary]、[ANI 分预算误差][ani_budget]、[零 J 比例][ani_zero]。

**5. BinDash 的推荐依据有两项独立证据：等字节精度更高，以及实际查询更快。** 这两项收益可以同时出现，所以它是普通集合任务里最值得优先考虑的候选之一。

在合成实验中，对每个预算、每个相同 regime，配对比较 BinDash 与 MinHash：

| pair budget（KiB） | BD 观测 MAE 更低 | BD/MH MAE 比值中位数 | 配对 95% 区间支持 BD 更低 |
| --- | --- | --- | --- |
| 2 | 84/84 | 0.509 | 59/84 |
| 8 | 84/84 | 0.484 | 64/84 |
| 32 | 84/84 | 0.490 | 71/84 |

每个预算下 84/84 个场景的观测 MAE 都较低；但只有表中部分场景的逐项 95% 区间完全支持差值方向。应写“所有已测场景的点估计更低”，而不是写“每个场景都统计显著”。区间也没有做多重比较校正。

真实 genome 上，六组 k/预算的总体结果如下：

| k | pair KiB | pairs | MinHash MAE | BinDash MAE | BD/MH |
| --- | --- | --- | --- | --- | --- |
| 20 | 2 | 200 | 0.010029 | 0.005177 | 0.516 |
| 24 | 2 | 200 | 0.009102 | 0.004513 | 0.496 |
| 20 | 8 | 200 | 0.005070 | 0.002434 | 0.480 |
| 24 | 8 | 200 | 0.004473 | 0.002538 | 0.567 |
| 20 | 32 | 200 | 0.002461 | 0.001239 | 0.503 |
| 24 | 32 | 200 | 0.002488 | 0.001199 | 0.482 |

这里每个 biological pair 先汇总三个 seeds，再对 200 pairs 统计。总体 MAE 约为 MinHash 的 48%–57%，支持合成结果在真实数据中延续。但它并非在每一对 genomes 上都更准；少数大误差或难点 pair 也可能影响总体均值。因此仍需查看分层与 P95。

机制上，当前 16-bit BinDash 在同样 4 KiB/object 中保存 2048 个 bin signatures，而 MinHash 只能保存 512 个 64-bit hashes。短签名会有随机碰撞，估计器作相应修正；densification 也会导致相关性，因此“四倍槽位”不是“四倍独立样本”的保证。机制能解释精度/字节优势的来源，真正支持推荐的是上述成对实验。

查询端则将 signatures 位打包，用位运算与 popcount 统计匹配。这比在每次比较中遍历完整 hash 集合或优化联合似然更容易实现高吞吐。N=10000、k=24、4 KiB/object 下，BinDash 查询 0.365 s，MinHash 1.696 s，SetSketch-fixed 75.336 s。前述原理与作者的[BinDash 算法说明](https://github.com/zhaoxiaofei/bindash#how-bindash-works)一致；性能倍数均来自本次 RabbitSketch 实现实测。

但推荐需要保留三个边界。第一，极低 J 时，即便 BinDash 的绝对 MAE 较低，相对误差仍可能很大；在真实 J≤0.001 的分组中，部分预算仍存在大量零估计，不能用小绝对误差声称能够精细分辨远缘关系。第二，低 J 的采样优势未必原样继承到真实 ANI，转换模型会引入另一层误差。第三，BinDash 当前 occupancy cardinality 在 genome 大集合上饱和，不能直接把普通相似度优势扩展为 directed containment 优势。

本次另测的 BinDash+HLL 已计入两份 payload 和两次构建；它可补充 cardinality，但没有稳定成为 containment 的最优方案。反复查询的收益也只在确实需要这些 pairs 时成立；本实验没有证明某种索引方法能减少候选 pair 数。

可写入 selection guide 的建议是：对普通 set similarity，在明确误差容忍度之后，优先评估 BinDash 的最小合格 bins；对已构建数据库的重复查询和大规模 all-vs-all，它的查询优势具有直接数据支持。一次性、构建为主的场景仍需与 MinHash 比较总成本。

源数据：[合成配对比较][synthetic_paired]、[真实 Jaccard 汇总][real_summary]、[按 J 分层结果][real_j_ranges]、[全矩阵成本][workload]。

**6. SetSketch-fixed 的推荐依据是高重叠联合估计中的精度收益，但其查询成本使推荐具有条件。** “SetSketch 很准”太宽泛；应说明在哪些 J 区间、哪些预算下更准，以及额外精度是否影响任务达标。

第一阶段 8 KiB/pair、相同大小集合的 synthetic 比较，下面每个 J 有 3 档 cardinality ×2 档 k，共 6 个场景。SS/BD<1 表示 SetSketch MAE 较低。

| true J | SS/BD MAE 比值中位数 | SS MAE 更低的场景 |
| --- | --- | --- |
| 0.001 | 3.955 | 0/6 |
| 0.01 | 1.929 | 0/6 |
| 0.05 | 1.287 | 1/6 |
| 0.1 | 0.946 | 4/6 |
| 0.25 | 0.999 | 3/6 |
| 0.5 | 0.911 | 4/6 |
| 0.75 | 0.869 | 4/6 |
| 0.9 | 0.730 | 5/6 |

这组数据不能支持“SetSketch 在所有相似度上更好”：低 J 时 BinDash 明显占优，而高 J 时 SetSketch 更有竞争力。真实 k=24、8 KiB/pair 的分层结果也表现出相近趋势：

| exact J 区间 | pairs | BinDash MAE | SetSketch MAE | SS/BD |
| --- | --- | --- | --- | --- |
| (0,.001] | 82 | 0.000101 | 0.001222 | 12.149 |
| (.001,.01] | 33 | 0.000963 | 0.002358 | 2.448 |
| (.01,.05] | 17 | 0.002486 | 0.003554 | 1.430 |
| (.05,.1] | 4 | 0.004897 | 0.003816 | 0.779 |
| (.1,.5] | 26 | 0.009501 | 0.008138 | 0.857 |
| (.5,1] | 38 | 0.004174 | 0.003191 | 0.764 |

J>0.5 的 k=24 分组有 38 对，SetSketch MAE 比 BinDash 低约 24%；k=20 对应 41 对，降幅约 31%。J∈(0.05,0.1] 的分组只有 4 对，不能据此设定一个通用切换阈值。分界应被描述为已测区间中的趋势，而不是“J>0.1 必须换 SetSketch”。

机制上，成对 register 的联合状态包含关于独有和共享元素的信息。当前修正版通过与实际 partitioned 构建相匹配的条件联合似然来估计共享部分，而不是只计算 register 相同个数。相似集合更强的相关结构可为估计提供信息；上述数据支持这种方法在部分高重叠条件下有效。这里的机制解释不等于证明它在所有高 ANI 数据中最优。

这轮的 SetSketch 必须明确标为修正版。先前实现把“一元素进入一个 register”的状态套入了不匹配的 dense rank-order 联合模型；修正使用包含空 register 和截断尾部的 partitioned 模型，marginal cardinalities 单独拟合。修正代码和验证说明见[模型验证][ss_validation]及[补丁][ss_patch]。报告不会把未经修正版本的旧异常当作算法固有特征，也不宣称这个补丁已经验证了所有索引功能。

精度收益的成本很高。4 KiB/object 的 10,000-genome 全矩阵查询中，SetSketch 为 75.336 s，BinDash 为 0.365 s，约相差 206 倍。当前 SetSketch 查询含数值似然优化，BinDash 主要执行位运算。这个倍率是本实现、此 workload 的观测，不是所有 SetSketch 实现必须承担的固定比例。

独立验证的选择结果提供了两个具体保留场景。下表的误差门槛和验证误差均为 ANI 百分点；两个请求都只有每 sketch 一次 incident comparison，成本使用实测矩阵吞吐模型。

| 数据 | ANI 区间 | pair KiB | 配置 | MAE 门槛（pp） | P95 门槛（pp） | 验证 MAE（pp） | 验证 P95（pp） |
| --- | --- | --- | --- | --- | --- | --- | --- |
| substitution | 90–95 | 8 | p=12 | 0.50 | 1.00 | 0.2503 | 0.5729 |
| mash142 | 99–100 | 8 | p=12 | 0.10 | 0.20 | 0.0588 | 0.1687 |

在受控 90–95% ANI 的请求中，MinHash K=512 的验证 P95 约 1.504 pp，超过 1 pp；BinDash 与 SetSketch 均能通过，而 SetSketch 的构建及准备比对应 BinDash 配置略低，所以在极少查询条件下被选中。查询次数增多后，数值优化的代价会改变选择。这里不是把 SetSketch 的速度说成普遍更快，而是说明误差约束和查询复用共同影响结果。

对于“真实 ANI 更准确”，结论要更保守。Mash held-out 99–100%、k=24、8 KiB/pair，SetSketch−BinDash 的配对差值如下。负值表示 SetSketch 误差更低，单位是百分点。

| 比较量 | SS−BD（pp） | 逐项 95% 区间（pp） | pairs | 涉及 genomes |
| --- | --- | --- | --- | --- |
| 采样误差 MAE | -0.0126 | [-0.0223, -0.0046] | 73 | 54 |
| 总 ANI 误差 MAE，对 FastANI | -0.0125 | [-0.0305, +0.0040] | 73 | 54 |

采样误差的区间支持 SetSketch 较低；相对 FastANI 的总误差区间仍包含零。把“采样改善”改写成“最佳 strain-level ANI 方法”会超出数据。alignment 参考更换造成的排序变化在第 11 段单独说明。

因此，推荐 SetSketch 的条件是：用户需要较准确的集合联合量，内存约束下其他配置可能不能达标，或当前满足精度的候选中它在少量查询下成本更低。对查询主导的大矩阵，通常需要非常明确的精度理由才能接受它的额外成本。它是应保留的精度候选，尚不是高 ANI 任务的无条件默认项。

源数据：[SS/BD 合成配对][synthetic_paired]、[真实 J 分层][real_j_ranges]、[held-out 选择记录][heldout]、[genome 层配对区间][paired_ani]。

**7. HLL 的推荐依据是 cardinality 任务中的简洁成本与足够精度，而不是普通 Jaccard 的排名。** distinct k-mer 数、共享 k-mer 比例和 abundance 是不同的目标量，不能用其中一个任务淘汰另一个任务的候选。

HLL 保存 hash 极值的压缩统计量，用固定数量 registers 估计 distinct count。它不保存每个元素的完整标识，也不会因同一 k-mer 重复出现而增加 distinct count。估计单个集合有多大，通常比准确推断两个集合共享多少元素更容易；本次 HLL 的 Jaccard 已使用其联合估计路径，不能把结果简单归因于只做了粗糙的容斥运算。

先看 k=24、p=12、4 KiB/object 的真实序列 cardinality。下表为相对误差百分比。

| 长度 | 每算法测量数 | HLL MAE（%） | HLL P95（%） | SS MAE（%） | SS P95（%） |
| --- | --- | --- | --- | --- | --- |
| 10 kb | 60 | 1.090 | 2.649 | 0.976 | 2.065 |
| 100 kb | 60 | 1.143 | 3.273 | 1.160 | 2.673 |
| 1 Mb | 60 | 1.106 | 3.160 | 1.235 | 3.310 |
| 5 Mb | 60 | 1.158 | 3.206 | 1.331 | 3.080 |
| 50 Mb | 15 | 1.014 | 1.937 | 1.041 | 1.730 |

HLL 的 MAE 在这批长度中约 1.0%–1.2%，与 SetSketch 接近。10 kb 上 SetSketch 的 MAE 略低，5 Mb 上 HLL 略低；这支持保留 HLL 的独立角色，不支持 SetSketch 全面替代 HLL。

同一配置的每对象构建与计数耗时如下。两部分分别取中位数再相加，单位微秒；这是 64 线程批次负载中的 elapsed time，不是独占单核 latency。

| 长度 | HLL build+card（μs） | SS build+card（μs） | SS/HLL |
| --- | --- | --- | --- |
| 10 kb | 123.75 | 341.45 | 2.76 |
| 100 kb | 523.85 | 906.07 | 1.73 |
| 1 Mb | 4502.43 | 4960.53 | 1.10 |
| 5 Mb | 21894.57 | 22215.17 | 1.01 |
| 50 Mb | 205252.08 | 212310.13 | 1.03 |

10 kb 时 HLL 约 124 μs、SetSketch 约 341 μs，前者的轻量构造更有吸引力；5 Mb 时约 21.9 ms 对 22.2 ms，差距明显缩小。长序列中扫描与 hash 处理占主要成本，不能把小对象优势直接放大成所有输入长度上的大幅加速。

更严格的精度要求下，还必须看 precision 参数。5 Mb、k=24 的完整 p 扫描如下；仍为相对误差百分比。

| p | bytes/object | HLL MAE（%） | HLL P95（%） | SS MAE（%） | SS P95（%） |
| --- | --- | --- | --- | --- | --- |
| 8 | 256 | 3.925 | 8.601 | 5.488 | 13.139 |
| 10 | 1024 | 2.511 | 6.282 | 2.851 | 6.472 |
| 12 | 4096 | 1.158 | 3.206 | 1.331 | 3.080 |
| 14 | 16384 | 0.661 | 1.901 | 0.590 | 1.360 |
| 16 | 65536 | 0.330 | 0.801 | 0.257 | 0.591 |

这个表包含一个必须保留的例外：p=14/16 时 SetSketch 的 MAE 和 P95 在这批数据中均更低。p=16 时，HLL MAE 约 0.330%，SetSketch 约 0.257%。因此，正确的推荐不是“cardinality 永远选 HLL”，而是“cardinality-only 优先从 HLL 开始；若误差门槛更严格，检查 SetSketch 是否在可用空间和成本下更合适”。

若用户已经为了 joint similarity 构建了 SetSketch，则应评估直接复用它的 cardinality，而不是无理由再构建 HLL；这一点是成本逻辑，不是本轮单独测过的组合工作流。若只需要 abundance/coverage，总量或频率结构本来就不属于 distinct count 的定义，应转到合适的 weighted metric。

源数据：[cardinality 全部精度结果][card_accuracy]、[构建与计数成本][card_build]、[HLL/SetSketch 第一阶段对照][hll_ss]。

**8. FracMinHash 的推荐依据是直接的 containment 采样与较低查询成本，限制则是小 query 的样本数量。** 大小比例是影响样本量的一个因素，不能单独充当选择开关。

定义 A 为 query，B 为 reference：

\[
C(A\rightarrow B)=\frac{|A\cap B|}{|A|}.
\]

FracMinHash 在兼容 hash 配置下使用同一全局阈值。若采样率 s=1/scaled，则较大集合自然保留更多 hashes，小集合保留更少 hashes。由同一阈值采样得到的交集与 query sample 数可直接构成 directed estimate。这与固定 K 的两个集合具有不同采样比例的情况不同。可比较的 scaled sketches 在需要时投影到共同、较稀疏的分辨率；该设计见 [sourmash 的 FracMinHash 说明](https://sourmash.readthedocs.io/en/latest/sourmash-internals.html)。

在完全包含 A⊆B 的情况下，必有 S(A)⊆S(B)；query sample 非空时，sample containment 等于 1。下面是最终 validation、k=24、8 KiB/pair、长度比 500、nominal C=1 的实际结果：

| 算法 | 有效/测量数 | 有效输出率 | 条件 MAE | 条件 P95 |
| --- | --- | --- | --- | --- |
| MinHash/FastKMV | 30/30 | 100.0% | 0.318978 | 1.000000 |
| BinDash | 0/30 | 0% | 不可用 | 不可用 |
| BinDash+HLL | 30/30 | 100.0% | 0.202113 | 0.686767 |
| HLL | 30/30 | 100.0% | 0.036521 | 0.135709 |
| SetSketch-fixed | 30/30 | 100.0% | ≈0 | ≈0 |
| FracMinHash | 25/30 | 83.3% | ≈0 | ≈0 |
| KSSD（实验诊断） | 1/10 | 10.0% | 0 | 0 |

FracMinHash 的有效样本上误差接近数值零，但有效输出比例只有 83.3%；SetSketch 同样接近零误差，且有效比例为 100%。KSSD 的条件 MAE 为零，却只在 10% 的测量中有有效输出。因此，绝不能只挑“零误差”这一列证明某算法最优。

完整包含不是唯一重要场景。nominal C=0.5、相同 k 与预算时：

| 长度比 | MH MAE | Frac MAE | Frac P95 | Frac 有效率 | SS MAE | SS P95 | BD+HLL MAE |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10 | 0.0572 | 0.0371 | 0.0876 | 100.0% | 0.0274 | 0.0653 | 0.0552 |
| 50 | 0.1561 | 0.0968 | 0.2261 | 100.0% | 0.0655 | 0.1650 | 0.1365 |
| 100 | 0.1835 | 0.1175 | 0.3084 | 100.0% | 0.0932 | 0.2335 | 0.1882 |
| 500 | 0.4811 | 0.2282 | 0.5002 | 86.7% | 0.1958 | 0.3929 | 0.3000 |

FracMinHash 相对 MinHash 的改进在 ratio=10/50/100/500 均可见；ratio=50 时 MAE 从约 0.156 降到 0.097。但 SetSketch 在这些行中更准，因此不能把这组结果描述为 FracMinHash 精度全面领先。ratio=500 时，FracMinHash 有效输出率为 86.7%，表中的 MAE 仍只针对可用测量。

查询成本解释了它为什么仍值得保留。以下取 ratio=50、C≈0.5、k=24、8 KiB/pair：

| 算法 | MAE | P95 | joint query 中位数（μs） |
| --- | --- | --- | --- |
| MinHash/FastKMV | 0.1561 | 0.3305 | 1.019 |
| FracMinHash | 0.0968 | 0.2261 | 1.012 |
| BinDash+HLL | 0.1365 | 0.2684 | 1.526 |
| HLL | 0.0747 | 0.1574 | 11.160 |
| SetSketch-fixed | 0.0655 | 0.1650 | 44.953 |

这里的时间是已准备 cardinality 缓存后的 joint Jaccard 查询中位数，单位 μs，不是含输入解析、构建以及所有后续派生操作的独立请求延迟。FracMinHash 约 1.01 μs，SetSketch 约 44.95 μs。二者的选择应取决于误差门槛是否允许 FracMinHash，而不是只看最低 MAE。

其限制可以用样本预算解释。令 r=|B|/|A|≥1，M 是两个 sketch 的 hash payload 总和；忽略随机大小波动，则：

\[
m_A\approx\frac{M}{8(1+r)}.
\]

8 KiB、r≈500 只对应约两个 query hashes。保持总空间不变时，再大的 reference 不能给 query 创造更多观测信息。固定大小 sketches 也受这一问题影响；FracMinHash 合理分配采样空间，不意味着它能绕开有限样本限制。

在简单随机抽样近似下，conditional containment 标准误差约为 sqrt(C(1−C)/m_A)。C≈0.5、希望 P95 绝对误差约≤0.1 时，用正态近似需要约 96 个 query samples；若 r≈500，对应约 376 KiB/pair。这只是预算量级推导，未在本次 8 KiB 网格上实测保证，也未计入所有有限样本、相关性和参数档位效应。

另一个容易误读的情况是 C=1 与 C=0.99 的区分：只有两个 query samples 时，C=0.99 下两个样本都来自共享部分的概率近似为 0.99²≈98%。因此，“完全包含时能给 1”并不证明能够精细识别 1% 的缺失。需要检查错误判满率或针对目标容差增大样本数。

本次主表用 joint Jaccard 和 native cardinalities 推导 containment，并按两侧 cardinality 共同约束交集；没有提供免费的 exact cardinality。对于同 scaled 的 FracMinHash，这与直接 sample ratio 一致。SetSketch/HLL 使用与 UnifiedQuery 一致的 common-intersection cap；原生低级方法的误差另外保存。BinDash+HLL 的两份 payload 和两次构建均计费。KSSD 的结果只是实验性 sample diagnostic，当前 public query contract 不支持 directed containment。

因此，推荐顺序是先确定 query 的 distinct cardinality 与误差要求，保证足够非空样本，再选 scaled 并检查 reference 的空间代价。对近完整包含、查询成本重要且 query 足够大的任务，FracMinHash 有明确位置；对极小 fragment、极小预算、部分 overlap 的精细估计，不能只凭 ratio 大就推荐它。对只有相近大小普通 Jaccard 的任务，本轮证据仍更支持优先比较 MinHash/BinDash。

![最终 containment：上排为有效输出条件下的 MAE，下排为有效输出比例][fig_containment]

源数据：[containment 精度与可用率][contain_accuracy]、[joint 查询成本][contain_cost]、[原生方法单独结果][contain_native]、[metric contract 审计][metric_audit]。

**9. ProbMinHash 的独立推荐理由是它能比较相对 abundance，而普通集合方法无法从相同 support 中恢复频率差异。** weighted 和 unweighted 模式必须分开评价。

例如两个样本只含相同的两个元素，weights 为 A=(90,10)、B=(10,90)。普通 set Jaccard=1，但相对组成明显不同。当前 ProbMinHash winner comparison 估计 probability Jaccard，而不是 min/max weighted Jaccard。对于正权重的共同坐标，定义为：

\[
J_P(a,b)=\sum_{i:a_i>0,b_i>0}
\left[\sum_{j\in\operatorname{supp}(a)\cup\operatorname{supp}(b)}
\max\left(\frac{a_j}{a_i},\frac{b_j}{b_i}\right)\right]^{-1}.
\]

上例 J_P=0.2。该例是依据 exact 定义计算的说明，不是声称来自某一对 CAMI 样本。若 B=10A，J_P=1，而 min/max weighted Jaccard=sum(min)/sum(max)=0.1。因此，整体乘一个系数不会被 J_P 当作组成差异；它并不表达绝对总 abundance 的不同。[ProbMinHash 作者说明](https://github.com/oertl/probminhash)也明确区分这两个指标。

CAMI II Marine 的 45 对官方 abundance profiles，每个配置 20 seeds，实际误差如下。误差均为 fraction，空间仅为 winner 查询 payload。

| m | KiB/object | 45 pairs ×20 seeds | J_P MAE | J_P P95 |
| --- | --- | --- | --- | --- |
| 128 | 1 | 900 | 0.032467 | 0.079247 |
| 512 | 4 | 900 | 0.015243 | 0.037465 |
| 2048 | 16 | 900 | 0.007745 | 0.018708 |
| 8192 | 64 | 900 | 0.004004 | 0.009788 |

register 数每增加四倍，MAE 大致减半。m=2048 时，MAE≈0.00774、P95≈0.01871；m=8192 时下降到约 0.00400 和 0.00979。这是在已知 exact J_P 下的直接准确性证据，而不只是“API 支持 weights”。720 组 uniform/log-normal/Zipf 合成 abundance 对照与缩放不变性测试提供了额外 sanity checks。

错误 metric 会造成另一种完全不同的现象。若把同一个 ProbMinHash 输出拿去对照 min/max Jaccard，得到：

| m | 正确 J_P MAE | 错用 min/max truth 的 MAE | 错用 min/max truth 的 P95 |
| --- | --- | --- | --- |
| 128 | 0.032467 | 0.142610 | 0.302160 |
| 512 | 0.015243 | 0.133680 | 0.276769 |
| 2048 | 0.007745 | 0.132663 | 0.258823 |
| 8192 | 0.004004 | 0.133367 | 0.250178 |

增大 registers 后，正确 J_P 的误差持续下降，而“对错 truth”后的 MAE 仍停在约 0.13。这不是通过继续加大 m 能解决的采样问题，而是估计目标不同。这个对照应保留在 Methods 或 Supplementary，防止论文把 probability Jaccard 泛称为 weighted Jaccard。

当前源代码中部分 weighted 注释及 weighted-containment 派生公式混用了 min/max 定义。一个具体反例是 a=(1,3)、b=(3,1)：exact J_P=1/2，若代入相应 min/max 推导会得到 containment=2/3，但 min-mass containment 实际为 1/2。本实验没有用这些无效派生量评分；weighted ANI 和这种 containment 不在推荐范围。详见[metric 审计][metric_audit]。

真实数据验证的范围也要写清楚：CAMI 输入是官方 genome abundance profiles，10 个 profile 的原文件与 setup archive 逐字节核对。它支持 profile-level composition comparison，不意味着已经验证从 noisy reads 估计 k-mer abundance 的整个流程。整体深度缩放不影响 J_P，也不等于能够自动消除物种间不均匀覆盖偏差。

对于 unweighted 模式，当所有元素 weight=1 时，J_P 退化为普通 set Jaccard，专门的 abundance 能力不再增加任务信息。第一阶段 k=24、8 KiB/pair、200 对真实 genomes：

| 算法 | Jaccard MAE | 每对构建（ms） | warm query（μs） |
| --- | --- | --- | --- |
| MinHash/FastKMV | 0.004473 | 26.469 | 0.645 |
| BinDash | 0.002538 | 39.176 | 0.368 |
| ProbMinHash-unweighted | 0.004479 | 51.042 | 0.479 |

ProbMinHash-unweighted 的 Jaccard MAE 与 MinHash 接近，构建中位数却约为其 1.93 倍；BinDash 的等字节 MAE 更低，构建也更快。ProbMinHash-unweighted 的 warm query 比 MinHash 快一些，因此不能随意宣称三项成本下被严格支配；但没有足够稳定的普通集合收益让它成为默认入口。

推荐应写成：需要 probability-Jaccard 意义上的相对 abundance 时保留 weighted ProbMinHash；只需 presence/absence 时优先比较 MinHash/BinDash/SetSketch。若任务要求绝对总量、min/max similarity 或 min-mass containment，应换成匹配的指标实现，不能直接复用 J_P 的名称和公式。

源数据：[weighted 精度][weighted_accuracy]、[错误 metric 对照][weighted_wrong]、[逐 pair 的 exact 指标][weighted_pairs]、[普通集合对照][real_summary]。

**10. KSSD 有部分精度收益，但本次没有建立独立默认推荐区域。** 这个结论需要同时看压缩程度、误差尾部和实际计算成本，不能只展示 dr=4 的小文件或某一项最小 MAE。

KSSD 通过 k-mer 子串空间采样与重编码进行 reduction。原方法允许兼容 sketches 上的集合运算，且不同 reduction levels 改变保留信息量；其原理见 [KSSD 原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7962209/)。这些设计价值不能自动视作本轮已经证明的性能优势。

本次恢复了 dr=2/3、half_subk=6，以及 dr=4、half_subk=7。后者需要更大的共享字典；初始化在每个批次只计一次，没有错误地向每个 genome 重复收取 setup。下表为 k=24、N=10000，全部时间单位为秒。

| 算法 | 参数 | 平均 bytes/object | 共享 setup（s） | 构建+cache（s） | 查询（s） | 合计（s） |
| --- | --- | --- | --- | --- | --- | --- |
| KSSD | dr=2 | 136871.61 | 0.928 | 6.962 | 98.0507 | 105.941 |
| KSSD | dr=3 | 8347.99 | 0.791 | 6.241 | 6.2852 | 13.317 |
| KSSD | dr=4 | 264.44 | 25.020 | 28.924 | 0.4106 | 54.355 |
| BinDash | bins=128 | 256.00 | 0 | 3.445 | 0.0362 | 3.482 |

dr=4 约 264 B/object，确实很紧凑，但对照的 BinDash 128 bins 只有 256 B/object。KSSD 的共享 setup 约 25.02 s，构建约 28.92 s，合计准备成本约 53.94 s；BinDash 约 3.45 s。KSSD 查询约 0.411 s，BinDash 约 0.0362 s。

因此，成本差距不能只归因于一次性字典初始化。即使讨论预先缓存共享 setup 的情形，仍需处理该批测量中的构建与查询差距；缓存后的完整工作流没有在本轮单独重新验证。也不能把当前 dr=4/half_subk=7 的成本外推为所有 KSSD 参数的固有开销。

精度方面，相近字节的紧凑配置在受控 substitution validation 中如下。单位为 ANI 百分点，k=24。

| ANI 区间 | KSSD MAE（pp） | KSSD P95（pp） | BD128 MAE（pp） | BD128 P95（pp） |
| --- | --- | --- | --- | --- |
| 90–95 | 0.8363 | 2.7668 | 0.9056 | 2.3418 |
| 95–99 | 0.3185 | 1.3125 | 0.4149 | 1.3001 |
| 99–100 | 0.0940 | 0.3444 | 0.0556 | 0.2185 |

KSSD 在 90–95% 和 95–99% 的 MAE 较小，但 P95 没有相应领先；99–100% 时 BinDash 的 MAE 和 P95 都较低。这说明“每个精度格子都被 BinDash 支配”不成立，也说明只看平均误差容易掩盖尾部取舍。

参数粒度也是一个实际因素：整数 dr 使 sampling rate 成倍变化，且 k/dr 会影响 tuple 的字节宽度。它不像固定-size sketches 那样总能在指定空间上限附近细调。这是当前 native grid 和表示的约束，不应写成所有 reduction-based sketches 的理论局限。

在当前 96 个约束选择请求中，KSSD 没有成为默认推荐配置。对于 containment，它没有 public query contract 支持；实验性 sample diagnostic 在最小 query、高 reduction 时经常为空，不能据少数有效样本的零误差建立推荐。索引、sketch union/intersection/subtraction 等其他 workflow 没有在本轮单独建立成本或精度优势，因此也不能凭原论文能力给它补一个已证实的“独占场景”。

现在适合写入论文的是：KSSD 在部分压缩/identity 配置上仍有精度竞争力，但在本次任务与约束范围内，尚未形成足以抵消其计算成本的默认推荐区域。保留 backend 和 supplementary 结果，降低一般用户的选择优先级。没有执行或证明跨全部参数与数据的全局淘汰定理。

源数据：[矩阵成本及 setup 分解][workload]、[按固定配置的 ANI 精度][ani_config]、[选择结果][heldout]、[containment 可用率][contain_accuracy]。

**11. ANI 推荐还必须处理采样之外的误差，以及阈值附近的实际判断。** 这一层决定“Jaccard 更准”能在多大程度上转换成 genome 用户的收益。

当前 native ANI 使用：

\[
g(J)=\left(\frac{2J}{1+J}\right)^{1/k}.
\]

另行报告的 Mash 型近似是 1+log(2J/(1+J))/k，二者不相同。对 J>0，native 变换的导数为 g(J)/[kJ(1+J)]；低 J 下对误差尤其敏感。更重要的是，exact J 本身也未必能够精确描述 alignment ANI。

带符号误差可以分解为：

\[
g(\hat J)-ANI_{ref}
=\{g(\hat J)-g(J_{exact})\}
+\{g(J_{exact})-ANI_{ref}\}.
\]

左侧是总误差；右侧分别是 sketch 采样误差与转换模型/参考定义的残差。绝对值之后不能把两个 MAE 相加，偏差抵消也可能让某个 sketch 的总 MAE 低于 exact-J baseline，这不表示它比 exact sets 含有更多信息。

Mash142 held-out pairs、k=24、8 KiB/pair，相对双向 FastANI 的结果如下，误差单位为百分点。

| 算法 | 95–99 MAE（pp） | 95–99 P95（pp） | 99–100 MAE（pp） | 99–100 P95（pp） |
| --- | --- | --- | --- | --- |
| MinHash/FastKMV | 0.2526 | 0.5736 | 0.0756 | 0.2148 |
| BinDash | 0.2634 | 0.4858 | 0.0712 | 0.1991 |
| SetSketch-fixed | 0.2131 | 0.4216 | 0.0588 | 0.1687 |
| FracMinHash | 0.2661 | 0.6165 | 0.0748 | 0.2506 |
| KSSD | 0.2684 | 0.6255 | 0.1012 | 0.2347 |

SetSketch 的点估计较低，但还要结合第 6 段的配对区间。相同 held-out 数据中，exact J 转换后的 MAE 在 95–99% 区间约 0.250 pp、99–100% 约 0.0596 pp；80–90% 和 90–95% 分别约 2.265 和 1.660 pp。后两个区间的明显残差不会单靠加大 sketch 消失。这些是观察到的模型/reference 差异，不是所有 genome pairs 不可突破的统一误差下界。

74 对固定分层样本的 NUCmer/dnadiff 审计中，71 对通过双侧 1-to-1 aligned-length fraction≥0.2。全部 74 对上，FastANI 与 NUCmer 的平均绝对差约 0.548 pp，最大约 4.102 pp。两者的对齐/匹配区域和定义不同，不能把其中一个强行称为另一个的绝对真值。

在 NUCmer 95–99% 区间的同一组 11 对 Mash audit pairs，固定同预算配置，排序就发生变化：

| 算法 | 共同 audit pairs | 对 FastANI MAE（pp） | 对 NUCmer MAE（pp） |
| --- | --- | --- | --- |
| BinDash | 11 | 0.2427 | 0.1590 |
| SetSketch-fixed | 11 | 0.1980 | 0.2245 |

这组数据支持的直接结论是“参考方法会影响排名”。它不支持用 11 对样本宣布 BinDash 普遍胜过 SetSketch，也不支持忽视 alignment 复核而只选对某算法有利的 FastANI 结果。

阈值任务也要直接评分。受控 substitution validation、k=24、8 KiB/pair 的分类误判率如下；95% 窗口为±1 pp，99% 窗口为±0.5 pp。

| 阈值 | 评价子集 | 每算法测量数 | MinHash 误判率 | BinDash 误判率 | SetSketch 误判率 |
| --- | --- | --- | --- | --- | --- |
| 95% | 含精确阈值点 | 150 | 16.00% | 8.67% | 10.67% |
| 95% | 排除精确阈值点 | 120 | 3.33% | 0% | 0% |
| 99% | 含精确阈值点 | 90 | 18.89% | 17.78% | 22.22% |
| 99% | 排除精确阈值点 | 60 | 0% | 0% | 0% |

同一真实 reference 的多个 mutants 和三个 sketch seeds 有相关性，n 是测量条数，不是独立 genome 数。排除恰好落在阈值的样本后，95% 窗口中 BinDash/SetSketch 在已测离散点上没有误判，而 MinHash 为 3.33%；99% 窗口中三者均没有误判。含阈值点时差异明显，因为真值恰好位于判定边界的小幅随机偏差就会改变标签。这里不能把零观测误判率解释为未来样本零风险，也不能从有限离散 mutation levels 声称已经证明任意小 ANI 差异都可分辨。

因此，selection guide 应明确用户要的是集合相似度、ANI 数值误差还是阈值分类。不同目标可能允许不同预算和参数，不能仅凭整体 MAE 决定所有任务。

源数据：[ANI 精度][ani_budget]、[exact-J 转换残差][ani_model]、[alignment 审计][alignment_accuracy]、[阈值分类][ani_thresholds]、[配对区间][paired_ani]。

**12. 查询复用把上述优势转成了可计算的选择条件。** 应先给定指标、误差和空间限制，再最小化合格配置的总成本。

对于需要构建 N 个对象、执行 P 个不重复 pair queries 的 workload：

\[
T=T_{shared}+N\,t_{build}+P\,t_{pair}.
\]

若把共享准备摊入该规模下的每对象成本 B，则 T/N=B+Q t_pair，其中 Q=P/N 是“按对象分摊的比较数”。每个 sketch 实际参与的比较次数 R=2P/N，因此也可写成 B+(R/2)t_pair。all-vs-all 时 Q=(N−1)/2、R=N−1。不能把 R 直接塞入 B+Q t_pair，造成两倍记账错误。已有 reference sketches 的单向搜索应按真实新建对象数和 pair 数计费，不机械套用所有对象都新建的模型。

采用 N=10000 实测吞吐得到的固定配置 crossover 如下：

| pair KiB | MH K | BD bins | Q*=P/N | R*=2P/N | all-vs-all N* |
| --- | --- | --- | --- | --- | --- |
| 2 | 128 | 512 | 15998.7 | 31997.4 | 31998.4 |
| 8 | 512 | 2048 | 4295.8 | 8591.6 | 8592.6 |
| 32 | 2048 | 8192 | 1087.8 | 2175.7 | 2176.7 |

8 KiB/pair 的 MinHash K=512 与 BinDash bins=2048，模型 N*≈8593；32 KiB/pair 的 K=2048 与 bins=8192，模型 N*≈2177。前者有 N=4096/10000 的实际反转检查，后者有 N=1024/4096 的检查。2 KiB/pair 的模型交点约 32000，超过实际测试的最大 N=10000，必须标为外推，不能说已经实测完成 32000 全矩阵。

这些交点依赖本机、线程数、数据构成和参数。用不同已测 N 校准，8 KiB 配置的交点约在 8.5k–10.6k 之间变化，这不是置信区间。query 三遍计时给出短期重复性，构建没有等量独立重跑，不能把接近的总时间点估计宣传为精确到个位 genome 的边界。all-vs-all 还是混合细菌数据；近缘 strain-only 查询组成可能改变成本。

性能内核实际计量的是 backend Jaccard comparison。ANI advisor 使用它作为比较成本模型，未包含逐 pair ANI 转换、输入解析、完整矩阵序列化及后续分析。它可以用于当前内核取舍，不能当作端到端 ANI 工具的绝对耗时承诺。有效 64-thread throughput 也不等于独立单请求 latency。

当前 advisor 使用 calibration 中 MAE≤epsilon、P95≤2epsilon、observed max pair payload≤预算来筛选固定配置，再最小化模型成本；选定后检查独立 validation。96 个请求中，45 个在 calibration 上选到候选，4 个所选候选未通过 validation，系统保留失败而没有利用 test 集再挑第二名。剩余无解仅表示已测配置中无解，不证明更大预算、不同参数或另一类方法也无解。

上述 41/45 是第二阶段内部 holdout 的结果。后续的独立多物种验证中，原 Mash-family 18 个冻结选择仅 6 个通过；这限制了固定误差阈值在新物种分布上的直接迁移。参见[新增报告的 advisor 外部检验](RabbitSketch2_robustness_validation.md#7-冻结-advisor-的外部检验)，其中保留了逐请求失败、P95 及 exact-J model residual。原型和旧 calibration 没有因新结果而重拟合。

一个实际通过 validation 的例子是 substitution、99–100% ANI、k=24、pair 上限 8 KiB、MAE≤0.1 pp、P95≤0.2 pp：

| 每 sketch incident reuse | 算法 | native 参数 | 实际最大 pair KiB | 验证 MAE（pp） | 验证 P95（pp） |
| --- | --- | --- | --- | --- | --- |
| 1 | MinHash/FastKMV | 512 | 8.0 | 0.0340 | 0.1107 |
| 10000 | BinDash | 512 | 2.0 | 0.0296 | 0.1362 |

高复用时选择的 BinDash bins=512 只需要 2 KiB/pair，说明预算是上限，算法不必总用满预算。这个例子也说明“固定 8 KiB 配置的 crossover”不等于“在误差约束下重新优化所有参数后的 crossover”。

![受控 substitution 下的选择图；×表示 calibration 选择未通过 holdout][fig_selection]

原型调用方式：

```bash
/home/bigssd/zt/rabbitsketch2_benchmark/venv/bin/python \
  benchmarks/regimes/advise.py \
  --family substitution --ani-bin '99–100' \
  --target-error 0.001 --pair-budget 8192 --queries-per-sketch 10000
```

这里 --queries-per-sketch 使用 incident reuse R，代码内部按 R/2 计费；--target-error 为 ANI fraction。该原型只在已测范围查证据，尚未接入生产库成为通用自动 advisor。

**13. 论文中的推荐可以按证据强弱表达，而不必强行让每个算法都有独占优势。** 每条判断都应同时给出比较对象、预算、任务和已测范围。

| 算法或问题 | 数据允许写出的结论 | 当前不应写出的结论 |
|---|---|---|
| MinHash | 在达到误差要求后，当前 FastKMV 在构建占主导的真实 FASTA workload 中具有成本优势 | 所有小规模任务必选 MinHash；低 ANI、小 K 也足够准确 |
| BinDash | 已测普通集合任务中有稳定的总体等字节精度收益，并在反复比较时获得查询成本优势 | 所有 pair、所有 metric、所有输入上都最优 |
| SetSketch-fixed | 部分高重叠场景采样误差更低；少量查询和精度约束下可成为合格候选 | 无条件最适合高 ANI；旧版与修正版结果可以混用 |
| HLL | Cardinality-only 的轻量方案有保留价值，SetSketch 没有全面覆盖它 | Cardinality 任何精度要求都必须选 HLL |
| FracMinHash | 有足够 query samples 时适合定向 containment，且查询成本低于当前数值联合估计器 | 只要 size ratio 大就一定更准；条件 MAE=0 就没有失败 |
| Weighted ProbMinHash | 对 probability Jaccard 的 synthetic 与真实 profile 估计有效 | 它估计任意 weighted similarity，或已经验证完整 reads abundance 工作流 |
| Unweighted ProbMinHash | 普通集合默认推荐优先级低，weighted 模式有另行验证的价值 | 普通集合结果足以否定整个 ProbMinHash backend |
| KSSD | 在当前参数、任务和成本约束下未形成默认入口，部分精度单元仍有竞争力 | 所有 KSSD variants 都被全局支配，或未测索引已经证明其优势 |
| OrderMinHash | 本轮按范围排除，未形成评测结论 | 它已被实验淘汰 |
| 算法筛选 | 用任务约束与实际成本建立选择边界，比统计宽泛 Pareto 前沿次数更有操作性 | 已满足“90% regimes 被支配”的正式全局淘汰判据 |

可用于论文的主张是：不同 sketch 的优势来自保留信息、空间分配和查询计算方式的不同；RabbitSketch2 的实测把这些差异转换成带指标、误差、空间和复用条件的选择建议。普通集合中 MinHash/BinDash 是主要成本竞争者；SetSketch 提供有代价的联合估计精度；HLL 和 weighted ProbMinHash 分别由 cardinality 与 probability-Jaccard 任务证明保留价值；FracMinHash 的 containment 推荐受 query sample 数限制；KSSD 的默认优先级由本次实测降低。无需把所有 backend 都包装成普遍不可替代。

**14. 本文所有数据均可从以下产物追溯。** 表格由现有 CSV 提取，未重新运行算法；解释性的概率公式与预算近似已在正文标明，不作为新的实验测量。

| 内容 | 入口 |
|---|---|
| 第一阶段精度与真实 FASTA 验证 | [第一阶段报告][stage1_report] |
| 第二阶段完整结果 | [第二阶段报告][stage2_report] |
| 指标、数据划分和补测约定 | [PROTOCOL.md][protocol] |
| API 指标与不适用派生量 | [METRIC_AUDIT.md][metric_audit] |
| 运行完成、SHA256 与代码快照 | [completion.json][completion] |
| 当前 workspace 的实验入口 | [main_regimes.cpp](benchmarks/regimes/main_regimes.cpp) |
| 当前 workspace 的分析脚本 | [analyze_suite.py](benchmarks/regimes/analyze_suite.py) |
| 可运行的选择原型 | [advise.py](benchmarks/regimes/advise.py) |
| SetSketch 修正和测试 | [VALIDATION.md][ss_validation]、[补丁][ss_patch] |

软件基线 commit：`8db783ea75f597068b72705272ae207c94194e83`。机器为 Intel Xeon Platinum 8375C，64 个可用 logical CPU；实验通过已有 main 的批次入口、单进程 64 OpenMP threads 顺序运行。原始结果根目录为 `/home/bigssd/zt/rabbitsketch2_benchmark/runs/`，包含 `regimes_20261001`、`regimes_supplement_20261001`、`containment_holdout_20261001`、`compressed_validation_20261001` 四个第二阶段目录。原始数据、版本与命令、逐批结果均保留。

独立图片与 PDF 可从原分析目录导出：[Jaccard 真实精度][fig_real]、[ANI 误差分解][fig_ani]、[cardinality][fig_card]、[containment][fig_containment]、[all-vs-all][fig_workload]、[weighted][fig_weighted]、[受控选择图][fig_selection]、[真实 genome 选择图][fig_selection_natural]。

为便于追踪本文使用的结果版本，下面记录关键源 CSV 的 SHA256。原始批次的完整校验值仍以冻结的 completion/manifest 为准。

| 源文件（相对于 analysis 根目录） | SHA256 |
| --- | --- |
| regimes_20261001/alignment_sketch_accuracy.csv | `81322d68cf39d87009fcddc04ea248a8fd714b2ddf1917e7598678574de23c51` |
| regimes_20261001/ani_budget_accuracy.csv | `7b6d35b66bf80340e2420da3007e39e5838b79ce59e17798bd70bdea14aa6749` |
| regimes_20261001/ani_config_accuracy.csv | `c9a983aac9183a6ccc90a9874cae6a9631cac73414830b9d17cc4c38831e5c06` |
| regimes_20261001/ani_threshold_errors.csv | `57c2139d449d6adaa0563a839fc33d9b81d5901756f96bac1453b3e14790c92b` |
| regimes_20261001/ani_zero_rate.csv | `8b890fd04d71aacbb248296b8cb6fd8f4dd715068cb831b46563fde0d551bb1c` |
| regimes_20261001/cardinality_accuracy.csv | `d31621ff117368773474ae3ed39f57d99d2468872c9718bdefd07fd87977cdd8` |
| regimes_20261001/cardinality_build.csv | `819573bd60f36b1809455ea97620f21eb103df3da63f80f2e4345dac71c7d34c` |
| regimes_20261001/containment_budget_accuracy.csv | `64445ca4d734b2f86740684db9fa79fa591969e87d62b2ba652da203e6901843` |
| regimes_20261001/containment_query_cost.csv | `26b4651983015113c96404646ec39bcb0fc8726a8e6c7ea27d86096cbdc62362` |
| regimes_20261001/heldout_selection.csv | `9a96c2b0380b7e9ca6fcd938fe5f8dfdd63449b91c044ea63b870dcc5609e50d` |
| regimes_20261001/paired_ani_intervals.csv | `ad5284c5a09b173208b9819af0cb6b61379163b0fa339522ed451b315274eb6f` |
| regimes_20261001/weighted_accuracy.csv | `6b9566c94db27f78a81c8d20c9af54414a3bb2ee79b8dbf20f78ddc681447987` |
| regimes_20261001/weighted_wrong_metric.csv | `81d190773c2746083034413a29621786deff6c9ff7679da890bcfc2e3d6d5eee` |
| regimes_20261001/workload_crossovers.csv | `e35069ad7149a5f1f776ecd4a81ae74cf757d25d16d9cf7c8205b8928f69d16d` |
| regimes_20261001/workload_measured.csv | `2e678b014fddcf90fc7b8b49000f079283ecab0746034408aa40a7b4dadd60cc` |
| systematic_20261001/real_jaccard_ranges.csv | `8d8b9c723833e83aba5edb492b88b3c6af0f8e3c40f47f55f9c1401aa73139e1` |
| systematic_20261001/real_summary.csv | `92e7ca2b1e8af7de9e282a3f1c7dcc9a8e7c08c418f46dc409c0b848c478ccda` |
| systematic_20261001/synthetic_paired_comparisons.csv | `95e814d23ba5e83382692131a15267c63286bca247d2218666fd36e67f3020ee` |

[stage1_report]: benchmarks/results/systematic_20261001/report.md
[stage2_report]: benchmarks/results/regimes_20261001/report.md
[protocol]: benchmarks/results/regimes_20261001/PROTOCOL.md
[metric_audit]: benchmarks/results/regimes_20261001/METRIC_AUDIT.md
[completion]: benchmarks/results/regimes_20261001/completion.json
[workload]: benchmarks/results/regimes_20261001/workload_measured.csv
[real_summary]: benchmarks/results/systematic_20261001/real_summary.csv
[synthetic_paired]: benchmarks/results/systematic_20261001/synthetic_paired_comparisons.csv
[real_j_ranges]: benchmarks/results/systematic_20261001/real_jaccard_ranges.csv
[ani_budget]: benchmarks/results/regimes_20261001/ani_budget_accuracy.csv
[ani_zero]: benchmarks/results/regimes_20261001/ani_zero_rate.csv
[heldout]: benchmarks/results/regimes_20261001/heldout_selection.csv
[paired_ani]: benchmarks/results/regimes_20261001/paired_ani_intervals.csv
[ss_validation]: benchmarks/results/systematic_20261001/VALIDATION.md
[ss_patch]: benchmarks/results/systematic_20261001/setsketch_partitioned_joint.patch
[card_accuracy]: benchmarks/results/regimes_20261001/cardinality_accuracy.csv
[card_build]: benchmarks/results/regimes_20261001/cardinality_build.csv
[hll_ss]: benchmarks/results/systematic_20261001/synthetic_hll_setsketch.csv
[contain_accuracy]: benchmarks/results/regimes_20261001/containment_budget_accuracy.csv
[contain_cost]: benchmarks/results/regimes_20261001/containment_query_cost.csv
[contain_native]: benchmarks/results/regimes_20261001/containment_native_accuracy.csv
[weighted_accuracy]: benchmarks/results/regimes_20261001/weighted_accuracy.csv
[weighted_wrong]: benchmarks/results/regimes_20261001/weighted_wrong_metric.csv
[weighted_pairs]: benchmarks/results/regimes_20261001/weighted_pair_accuracy.csv
[ani_config]: benchmarks/results/regimes_20261001/ani_config_accuracy.csv
[ani_model]: benchmarks/results/regimes_20261001/ani_model_error_summary.csv
[alignment_accuracy]: benchmarks/results/regimes_20261001/alignment_sketch_accuracy.csv
[ani_thresholds]: benchmarks/results/regimes_20261001/ani_threshold_errors.csv
[fig_containment]: benchmarks/results/regimes_20261001/containment_boundaries.png
[fig_selection]: benchmarks/results/regimes_20261001/selection_substitution.png
[fig_selection_natural]: benchmarks/results/regimes_20261001/selection_mash142.png
[fig_real]: benchmarks/results/systematic_20261001/real_accuracy.png
[fig_ani]: benchmarks/results/regimes_20261001/ani_error_decomposition.png
[fig_card]: benchmarks/results/regimes_20261001/cardinality_tradeoff.png
[fig_workload]: benchmarks/results/regimes_20261001/workload_scaling.png
[fig_weighted]: benchmarks/results/regimes_20261001/weighted_validation.png