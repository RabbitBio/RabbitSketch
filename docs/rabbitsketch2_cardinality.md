# RabbitSketch 基数估计 API

更新日期：2026-11-17。

下表六种算法的原生 C++ 类、Python 类和高层 `BuiltSketch` 均提供 `cardinality()`，返回
不同 canonical k-mer 数量的浮点估计值。输入频次不会被当作不同元素数量。
本轮实验中的 MinHash 是 `FastKMV`；Python 的历史 `MinHash` 保留其原有估计器。

| 算法 / Python 类 | 方法 | 适用条件与返回值 |
|---|---|---|
| HLL / `HyperLogLog` | 原生 Ertl MLE | `finalize()` 计算并缓存，后续读取 O(1) |
| `SetSketch` | 修正后的分区寄存器估计器 | `finalize()` 计算并缓存，后续读取 O(1) |
| MinHash / `FastKMV` | 满 K 时 `(K-1) × KEY_MAX / tau_K` | 未满时返回保留的不同哈希数；哈希碰撞除外 |
| `FracMinHash` | 保留哈希数 / 实际采样概率 | 首次查询会 finalize；空样本返回 0 |
| KSSD / `Kssd` | 保留降维编码数 / `admissionProbability()` | 使用当前 dr 的实际采样概率；支持 32/64 位编码，空样本返回 0 |
| `BinDash` | `-B × log(E/B)`，E 为填充前的空桶数 | 空输入返回 0；E=0 返回正无穷，表示饱和、估计不可用 |

`capabilities_for(algorithm).supports(Capability.Cardinality)` 表示该算法有此
入口，不保证每个对象都有可用估计。ProbMinHash、OrderMinHash 和历史 MinHash
的统一 Cardinality capability 为 false，调用其 `BuiltSketch.cardinality()` 会抛出
`RuntimeError`（C++ `logic_error`）。历史 MinHash 保留原有的原生基数方法。
ProbMinHash 不提供基数接口；`total_weight` 和 winner 数不能代替不同元素数量。

## 空样本、饱和及空间

FracMinHash/KSSD 的 0 是空**样本**的估计值，不足以证明原始集合为空；样本数和
实际空间应同时报告。BinDash 的 infinity 应统计为不可用；不能过滤这些记录后将剩余样本误差当作全部输入的精度。

查询数据、完整对象状态、构建峰值内存是不同概念。HLL/SetSketch 的见证数据和
BinDash 的未压缩桶最小值会增加完整状态空间；FracMinHash/KSSD 参数不保证固定
字节预算。该接口不改变原生状态表示。

## 调用示例

```python
import math
import rabbitsketch as rs

sketch = rs.FastKMV(size=1024, kmer_size=21)
sketch.update("ACGTACGTACGTACGTACGTACGTACGT")
sketch.finalize()
estimate = sketch.cardinality()
assert math.isfinite(estimate)

# 高层 builder 返回的 BuiltSketch 使用同一入口。
config = rs.SketchConfig()  # 默认 FastKMV，k=21
for result in rs.build_fastx("sample.fa", [config]):
    print(result.sketch.cardinality())
```

```cpp
#include "fastkmv.h"
#include "api/SketchBuilder.h"

Sketch::FastKMV sketch(1024, 21, 42);
const char sequence[] = "ACGTACGTACGTACGTACGTACGTACGT";
sketch.update(sequence, sizeof(sequence) - 1);
sketch.finalize();
double estimate = sketch.cardinality();
auto built = Sketch::API::BuiltSketch::fromFastKMV(sketch);
double same_estimate = built.cardinality();
```

总耗时应测量同一个对象从构造到首次 `cardinality()` 返回的完整过程。
HLL/SetSketch 已在 `finalize()` 中计算基数，不应再重复加一次冷查询时间。
