# RabbitSketch 2 实用 API 指南

更新时间：2026-08-15

本文描述当前高层生产接口：统一配置、一次 FASTX 读取、多算法内存构建和统一查询。

## 1. 数据流

```text
FASTA/FASTQ
    -> FastxReader（plain/gzip/stdin）
    -> SketchConfig[] + MultiSketchBuilder
    -> BuildResult[] / BuiltSketch
    -> Query::Result
```

高层代码优先使用 `SketchConfig`、`buildFastxFiles`、`BuiltSketch` 和
`Query::Request/Result`。原有算法类仍可用于直接 update、merge、projection 或原生
查询，但调用方不必为不同算法重复实现 FASTA 读取和输入清洗。

## 2. 算法契约

| Algorithm | sampling/resolution | 合法共同精度 | 主要 metric |
|---|---|---|---|
| `LegacyMinHash` | `BottomK`, K>=2 | 较小 K | set metrics、ANI |
| `FastKMV` | `BottomK`, K>=2 | 较小 K | set metrics、ANI |
| `FracMinHash` | `Scaled`, scaled>0 | 较大 scaled | set metrics、ANI |
| `HyperLogLog` | `RegisterPrecision`, p=4..20 | 较小 p | cardinality/set metrics、ANI |
| `SetSketch` | `RegisterPrecision`, p=4..16 | 相同 p | cardinality/set metrics、ANI |
| `KSSD` | `KssdReduction` | 较大 drlevel | Jaccard/distance/ANI |
| `ProbMinHash` | `AlgorithmNative`, registers>=2 | 相同 registers | set 或 weighted metrics |
| `BinDash` | `AlgorithmNative`, bins 为 64 的正倍数 | 相同 bins | set metrics、ANI |
| `OrderMinHash` | `AlgorithmNative`, m>=1 | 较短 sample prefix | order similarity/distance |

所有算法的 `kmer_size` 为 1..32。KSSD 还要求偶数 k、合法 half-sub-k 和
reduction level；Legacy MinHash 的 seed 必须能放入 32 bit。构建前可用
`capabilitiesFor(algorithm)` 查询能力，最终 pair 是否可执行仍由 query planner 决定。

## 3. FASTX 输入语义

- `FastxReader(path)` 流式读取多行 FASTA/FASTQ，接受未压缩、gzip 和 `"-"` 标准输入。
- 截断 FASTQ、序列/质量长度不等和非法 Phred+33 字节会明确报错。
- 当前高层 builder 支持 DNA/RNA；RNA 的 `U` 转为 `T`。其他 molecule 枚举暂不接收。
- 输入先大写，再处理 ambiguous policy、homopolymer compression 和质量过滤。
- `RejectRecord` 遇到非 A/C/G/T 立即拒绝；`SkipKmer` 跳过跨越歧义碱基的窗口。
- `minimum_base_quality>=0` 跳过包含低质量碱基的整个 k-mer；FASTA 不能请求质量阈值。
- abundance 模式仅用于 ProbMinHash：先聚合 canonical k-mer 频次，再应用 abundance
  上下限并以 frequency weight 更新。
- OrderMinHash 必须使用 `OneSketchPerRecord + RejectRecord`，且不能使用会切断序列的
  质量过滤。

`BuildStats` 报告 records、input bases、candidate/accepted k-mers、歧义/质量/丰度
跳过量和 distinct k-mers。

## 4. 聚合和结果

- `OneSketchPerFile`：每个文件、每个对应 config 一个结果；
- `OneSketchPerRecord`：每条 record、每个对应 config 一个结果；
- `OneSketchPerCollection`：全部输入形成一个逻辑 collection。

`buildFastxFiles(paths, configs)` 保持输入与 config 的稳定顺序。调用方应读取
`BuildResult.source`、`record_name`、`label` 和 `config`，不要从 label 反解析关系。

```cpp
#include "api/SketchBuilder.h"

Sketch::API::SketchConfig config;
config.algorithm = Sketch::API::Algorithm::FastKMV;
config.sampling = Sketch::API::SamplingMode::BottomK;
config.kmer_size = 21;
config.resolution = 1024;

auto results = Sketch::API::buildFastxFiles({"a.fq.gz", "b.fq.gz"}, {config});
auto identity = results.at(0).sketch.query(results.at(0).sketch);
```

`MultiSketchBuilder` 是 move-only；一次 `finish()` 后 sealed，再次 update/finish 会抛出
`logic_error`。它适用于调用方已经持有 `FastxRecord` 的流式管线。

## 5. 查询和错误处理

```cpp
Sketch::Query::Request request;
request.metric = Sketch::Query::Metric::Jaccard;
auto result = left.query(right, request);
if (!result.estimate_available)
    throw std::runtime_error(result.plan.reason);
```

`Result.value` 是所请求 metric；可推导的 Jaccard、distance、ANI、cardinality、
containment、weighted 或 order 字段会同时填充，其余保持 NaN。拒绝结果不是零相似度。

planner 检查 backend、k、seed、canonicalization、hash profile、weight semantics、
backend 参数和 resolution。同 backend 不代表必然兼容，不同 backend 明确拒绝。
frequency/user-supplied ProbMinHash 应请求 weighted metric；OrderMinHash 只接受 order
metric。

## 6. Python、线程和部署

Python 对应接口为 `SketchConfig`、`FastxReader`、`MultiSketchBuilder`、
`build_fastx/build_fastx_files`、`BuiltSketch` 和 `QueryRequest/QueryResult`。耗时的
native update、build 和 query 会释放 GIL，但同一对象不承诺可被多个线程同时修改。

默认 Release 不使用 `-march=native`。x86 热路径运行时选择 scalar、SSE2、AVX2 或
AVX-512BW，其他架构使用 scalar。Python 的 `runtime_info()` 可查看选择，环境变量
`RABBITSKETCH_SIMD=scalar` 可强制安全路径。OpenMP 可选；缺少时保持正确的串行行为。

CMake 安装后使用 `find_package(RabbitSketch 2 CONFIG REQUIRED)` 和
`RabbitSketch::rabbitsketch`，也可使用 `rabbitsketch.pc`。Python 包采用 PEP 517 构建，
没有额外 Python runtime dependency。
