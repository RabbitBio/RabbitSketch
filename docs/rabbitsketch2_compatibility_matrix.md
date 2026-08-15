# RabbitSketch 2 生产兼容矩阵

更新时间：2026-08-15

本表是调用方在不读取原始序列时可以依赖的比较规则。`supported` 表示已有原生
估计器；`projection` 表示可以从较高信息量状态精确得到共同状态；`same only`
表示参数不同就拒绝。不同 backend 一律 `unsupported`。

| backend | 同 backend | resolution 变化 | 主要限制 |
|---|---|---|---|
| Legacy MinHash | supported | `min(K)` projection | k、seed、canonicalization 相同 |
| FastKMV | supported | `min(K)` projection | k、seed、hash profile 相同 |
| FracMinHash | supported | `max(scaled)` projection | k、seed、hash profile 相同 |
| HyperLogLog | supported | `min(p)` folding | k、seed、hash profile 相同 |
| SetSketch | supported | same precision only | `a/base`、rank orientation 相同 |
| KSSD | supported | `max(drlevel)` projection | shuffle 与 reduction 参数相同 |
| ProbMinHash | supported | same register count only | seed、Top-L、weight semantics 相同 |
| BinDash | supported | same bin count only | seed、b-bit 参数相同 |
| OrderMinHash | supported | shorter sample prefix | k、l、seed、orientation 相同 |

## 状态含义

- `Supported`：可以执行请求的 metric；
- `Unsupported`：算法或 payload 没有该估计器，或 pair 是不同 backend；
- `Incompatible`：本应属于同一算法族，但 metadata/参数不匹配。

任何非 `Supported` 状态都会令 `estimate_available=false`，数值字段保持 NaN，
具体原因写入 `plan.reason`。

## Metric 语义

- Legacy MinHash、FastKMV、FracMinHash、HLL、SetSketch 和 BinDash 面向普通集合；
- KSSD 当前提供 Jaccard、distance 和 ANI；
- unweighted ProbMinHash 面向普通集合，weighted ProbMinHash 面向明确的权重分布；
- OrderMinHash 面向顺序相似性，不与集合 Jaccard 混用。

## 接入原则

新增 backend 或 query path 时，必须同时补齐 metadata 检查、planner 规则、原生查询、
高层 `BuiltSketch` 查询和至少一个拒绝测试。没有可验证 estimator 的组合应保持
unsupported，不能用零值或临时转换掩盖。
