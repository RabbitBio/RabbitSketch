### SetSketch 小集合异常的独立复现

使用原来的 main 入口、64 线程、相同输入种子及原生库，给冻结 worker 的副本增加诊断输出。p=12/14/16 各复现 20 个原始 pair，全部默认估计值、输入摘要和 payload 与原记录完全一致。诊断调用改变了查询顺序，因此这些新耗时不参与性能统计。

固定 k=20、|A|=10,000、|B|=100,000，真实 J=0.001001001。

| p | 默认 J 估计均值 | 默认估计器 MAE | 实例 inclusion–exclusion MAE | 双方均空的寄存器比例 |
|---|---:|---:|---:|---:|
| 12 | 0.0020056878 | 0.0019425406 | 0.0055287719 | 0.00% |
| 14 | 1.0013561e-09 | 0.001001 | 0.0024851169 | 0.13% |
| 16 | 0.055715509 | 0.054714508 | 0.0006740214 | 18.90% |

p=16 时默认估计器的 MAE 为 0.054715，现有实例 inclusion–exclusion 路径为 0.000674；说明该异常与默认查询路径有关，不能简单归咎于寄存器容量不足。p=12/14 下 inclusion–exclusion 又不占优，因而不能未经验证就全局替换估计器。

源码核查：当前基数估计已使用 partitioned register 的边际似然；默认 joint 内核仍基于寄存器相等/大于/小于计数，并在双方均空或饱和的寄存器超过 30% 时回退。本例双方均空比例约 18.9%，20 次均未触发回退。这支持“低负载区间的 joint 模型或回退策略需要修正”的假设；尚未完成修复验证。

不能仅把 30% 阈值调低：单独诊断现有静态 `jaccardFromCores` 路径时，p=16 的 MAE 为 2.037，输出还可能超出 [0,1]。该静态路径用 `factor/sum` 计算 union，与实例路径的基数处理不同。它在这组原始测量中并未实际触发；原始 Jaccard 估计仍全部位于 [0,1]。

定位：[joint 内核](../../../results/README.md#external-archive)，[实例 inclusion–exclusion](../../../results/README.md#external-archive)，[静态 fallback](../../../results/README.md#external-archive)。

复现入口：`benchmarks/screening/diagnose_setsketch.py`；完整输出与编译/链接摘要在 `setsketch_diagnostic/`。原生库、原始记录和预注册结论均保留。
