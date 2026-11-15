# RabbitSketch 统一操作 API

更新日期：2026-11-15。

七种算法已经接入 `merge`、`project`、`store_hashes` 和 `index_dict`：MinHash
（源码及 Python 类名 `FastKMV`）、FracMinHash、HLL、SetSketch、ProbMinHash、
BinDash、KSSD。Python 的历史 `MinHash` 类保持原义，不在这七种实现之列。

## 1. 已实现的入口

| 操作 | Python | C++ | 契约 |
|---|---|---|---|
| 合并 | `sketch.merge(other)` | 原生类或 `BuiltSketch::merge` | 返回新 Sketch，验证算法、种子、哈希与算法特定参数；不修改输入 |
| 投影 | `sketch.project(target)` | 原生类或 `BuiltSketch::project` | 返回指定信息量的状态；只允许下表定义的方向 |
| 导出 | `sketch.store_hashes()` | 原生 `storeHashes`，或 `BuiltSketch::storeHashes` | 返回原生摘要数据的副本，具体含义见下文 |
| 索引 | `index_dict(sketches, num_threads=1)` | `API::indexDict` | 返回 CSR `SketchIndex`，支持七种原生类集合及 `BuiltSketch` 集合 |
| 类级索引入口 | `FastKMV.index_dict(...)` 等七个类及 `BuiltSketch.index_dict(...)` | `BuiltSketch::indexDict` 或泛型 `API::indexDict` | 与模块入口共用实现 |
| 基数扩展 | `cardinality()` | `BuiltSketch::cardinality` | 统一 API 只对 HLL、SetSketch 开放；其他算法调用会明确报错 |

原生 Python 类中的 `cardinality` 也仅由 HLL、SetSketch 提供。原有 C++ 基数
辅助计算仍用于相似度和 containment 推导；FracMinHash 的历史
`cardinality_at` 保留兼容，不属于统一 `cardinality` 接口。

高层配置与生命周期仍使用 `SketchConfig`、`MultiSketchBuilder.update/finish`、
`build_fastx/build_fastx_files` 和 `BuiltSketch.query`。C++ 可包含
`api/SketchOperations.h` 使用原生集合到统一索引的适配。

## 2. 投影和合并语义

| 算法 | `project` 的目标参数与约束 | 合并方式 |
|---|---|---|
| MinHash（FastKMV） | 较小 K | 保留合并后的 Bottom-K，允许投影到共同的较小 K |
| FracMinHash | 较大 scaled | 在共同采样阈值下合并、去重 |
| HLL | 较小 p | 精确 folding；相同寄存器取最大值 |
| SetSketch | 较小 p，需要保留完整分区最大哈希值 | 相同精度逐寄存器取最大值，同时合并分区最大哈希值 |
| ProbMinHash | 较少寄存器，截取**原始采样族的前缀** | 相同采样族、相同规模下逐寄存器取最小竞争值，并保留对应 winner |
| BinDash | 较少槽位；目标须为 64 的倍数且整除源槽位数 | 对未填充的完整桶最小值求最小值，再重新填充与压缩 |
| KSSD | 较大 drlevel | 投影到共同降维级别后，对降维编码求集合并集 |

ProbMinHash 投影保留 `source_register_count` 和原始 Top-L 参数，并将采样族写入
兼容性指纹。从 1024 个寄存器投影得到的 256 个寄存器，不能当作独立构建的
256 寄存器 Sketch 混用；同一个原始采样族投影后的对象可查询、合并和索引。
投影对象封存，避免后续更新错误地使用较小规模的随机竞争参数。

ProbMinHash 合并沿用原生的最小竞争值语义，不能代替同一 k-mer 在两批输入中的
丰度相加。需要合计丰度时，应在带权构建前聚合频数。

SetSketch 的完整分区最大值包含寄存器等级为零的分区，也处理同等级内的更大
哈希值，因此投影不会仅凭量化后的寄存器等级猜测目标状态。

## 3. 构建状态与内存

`SketchConfig.retain_operation_state` 默认为 `true`：

- HLL 保留候选见证哈希，供倒排索引使用；folding 后保留支持折叠寄存器的见证。
- SetSketch 保留每个分区的完整最大哈希值，供投影与候选索引使用。
- BinDash 在 `finalize` 后保留未填充的完整桶最小值，供合并和投影使用。

这些辅助数组分别需要约每寄存器或每槽位 8 字节额外空间。仅需紧凑比较状态时，
可以在高层构建配置中关闭 `retain_operation_state`。原生 BinDash 构造器也提供
`retain_operations=False`；原生 HLL、SetSketch 的 C++ 构造器保留见证开关。

关闭或丢弃相应辅助状态后，HLL/SetSketch 的见证索引、SetSketch 的降精度投影、
BinDash 的合并和降槽位投影会报出明确的状态缺失错误。通过 `fromPacked` 或
`fromRegisters` 恢复的历史紧凑状态也适用此限制，不能从有损摘要恢复已丢弃信息。

## 4. store_hashes 的数据含义

| 算法 | 返回内容 |
|---|---|
| MinHash、FracMinHash | 有序保留哈希值 |
| KSSD | 有序降维编码；Python 的统一入口兼容 32 位与 64 位编码 |
| HLL、SetSketch | 寄存器等级，按原始位置排列 |
| ProbMinHash | 每个寄存器的 winner 标识，保留位置及空值 |
| BinDash | 位打包后的签名数组 |

此接口沿用原 KSSD 的“返回摘要数组”用法，不写文件，也不等于完整状态序列化。
例如 ProbMinHash 的竞争值、SetSketch 的分区完整最大值和 BinDash 的未压缩桶值
不在这个数组中；不能仅凭导出数组承诺恢复全部状态操作。

## 5. index_dict 的使用

```python
import rabbitsketch as rs

a = rs.FastKMV(size=1024, kmer_size=21)
b = rs.FastKMV(size=1024, kmer_size=21)
a.update(sequence_a)
b.update(sequence_b)
a.finalize()
b.finalize()

combined = a.merge(b)
smaller = combined.project(512)
values = smaller.store_hashes()

index = rs.index_dict([a, b], num_threads=4)
# 等价入口：rs.FastKMV.index_dict([a, b], num_threads=4)
for match in index.pairs(min_similarity=0.2, num_threads=4):
    print(match.left, match.right, match.estimate.value)
```

索引提供 `postings(value, position=0)`、`candidates(sample)`、
`search(sample, min_similarity)`、`pairs(min_similarity, num_threads)`，以及
`key_count`、`posting_count`、`approximate_candidates` 属性。样本编号对应输入
集合中的位置。查询限于该集合；构建时删除单例键，不作为外部新样本查询数据库。

构建复用分桶排序、单例删除和 CSR 连续存储。匹配使用完整的值与位置组合键，
不会通过 XOR 合并位置而引入无法区分的键碰撞。批量查询采用线程局部时间戳计数，
仅处理共享索引键的候选，并调用原生估计器验证；`pairs` 只输出一次每个无序对。

HLL、SetSketch 使用实际输入见证哈希生成候选，`approximate_candidates=True`，
不承诺找出全部超过阈值的样本对。其他五种算法的共享键覆盖正相似度的原生比较，
但 `min_similarity=0` 也不会枚举没有共享键的零相似度样本对。
带权 ProbMinHash 自动使用 winner 匹配和 `WeightedJaccard`；其余模式使用普通
Jaccard。阈值为相似度阈值，不是 Mash 距离阈值。

旧的 KSSD 文件工作流 `index_dict` 重载仍然保留。

## 6. 实现与验证入口

- [统一构建与操作](../src/api/SketchBuilder.h)
- [索引 API](../src/api/SketchIndex.h)、[CSR 索引执行](../src/api/SketchIndex.cpp)
- [原生类型适配](../src/api/SketchOperations.h)
- [Python 绑定](../src/pybind.h)
- [C++ 操作回归](../tests/sketch_operations_test.cpp)
- [Python 操作回归](../tests/test_sketch_operations.py)
