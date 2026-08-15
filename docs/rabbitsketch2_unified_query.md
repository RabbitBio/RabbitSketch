# RabbitSketch 2 统一查询契约

更新时间：2026-08-15

## 1. 统一的边界

统一查询层不要求九种算法产生相同数组，也不声称不同算法的 sketch 可以直接比较。
它统一三件事：

1. `RankMetadata`：输入、哈希、seed、backend、权重和 resolution 契约；
2. `Query::Request/Plan`：请求 metric，并在执行前给出支持状态和方法；
3. `Query::Result`：用同一 schema 返回结果、派生量、有效样本数和拒绝原因。

## 2. 生产规划矩阵

| 输入 | 状态 | resolution 规则 | metric 边界 |
|---|---|---|---|
| FastKMV ↔ FastKMV | supported | `min(K1,K2)` 精确 prefix | set metrics、ANI |
| Legacy MinHash ↔ Legacy MinHash | supported | `min(K1,K2)` 精确 prefix | set metrics、ANI |
| FracMinHash ↔ FracMinHash | supported | 投影到 `max(scaled1,scaled2)` | set metrics、ANI |
| HyperLogLog ↔ HyperLogLog | supported | folding 到 `min(p1,p2)` | cardinality/set metrics、ANI |
| SetSketch ↔ SetSketch | 相同 p supported | 不同 p 拒绝 | cardinality/set metrics、ANI |
| KSSD ↔ KSSD | supported | 投影到 `max(drlevel1,drlevel2)` | Jaccard/distance/ANI |
| ProbMinHash ↔ ProbMinHash | 相同 register 数 supported | 不降精度 | unweighted set 或 weighted metrics |
| BinDash ↔ BinDash | 相同 bin 数 supported | 不降精度 | set metrics、ANI |
| OrderMinHash ↔ OrderMinHash | supported | 共同较短 sample prefix | order similarity/distance |
| 不同 backend | unsupported | 不执行 | 无 |

同 backend 仍必须匹配 k-mer size、seed、canonicalization、hash profile、weight
semantics 和 backend 参数指纹。metadata 不匹配返回 `incompatible`；缺少已支持的
估计器或 metric 返回 `unsupported`。

## 3. 请求和返回值

```cpp
Sketch::Query::Request request;
request.metric = Sketch::Query::Metric::Jaccard;
request.confidence_level = 0.95;

auto result = Sketch::Query::query(left, right, request);
if (!result.estimate_available)
    throw std::runtime_error(result.plan.reason);
```

`Result.value` 对应所请求 metric。Jaccard、distance、ANI、cardinality、intersection、
containment、weighted 和 order 字段只在算法能够推导时填充，否则保持 NaN。
调用者必须先检查 `estimate_available`，不能把拒绝结果当作零相似度。

planner 的执行方法只有九种 native path：`FastKmvNative`、`HllNative`、
`SetSketchNative`、`KssdNative`、`ProbMinHashNative`、`BinDashNative`、
`LegacyMinHashNative`、`FracMinHashNative` 和 `OrderMinHashNative`。

## 4. 特殊 metric

- weighted ProbMinHash 与普通集合估计不是同一 estimand。frequency 或用户权重对象
  只能请求 weighted Jaccard/containment；unweighted 对象可请求普通 set metrics。
- OrderMinHash 只接受 order similarity/distance。
- KSSD 当前不提供依赖 cardinality 的 intersection/containment。
- ANI 和 Mash distance 是由 k-mer similarity 模型换算的派生量，不是 alignment。

## 5. 高层内存对象

`buildFastx/build_fastx` 返回的 `BuiltSketch` 保存算法原生状态并复用同一个 planner：

```python
request = rs.QueryRequest()
request.metric = rs.QueryMetric.Jaccard
result = left.sketch.query(right.sketch, request)
if not result.estimate_available:
    raise ValueError(result.plan.reason)
```
