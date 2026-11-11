#!/usr/bin/env bash
# bench_233_index_vs_noindex.sh
#
# 对 233.list（全量或截断），逐算法对比：
#   allpairs (-l)  vs  inverted-index (-l --index)
# 采集 wall time（秒）和 max RSS（KiB），结果汇成一张 CSV。
#
# 用法（在 examples/ 目录下执行）：
#   bash bench_233_index_vs_noindex.sh [选项]
#
# 常用选项：
#   --exe     <path>   rabbitsketch 可执行（默认 ./rabbitsketch）
#   --list    <path>   文件列表（默认 ./233.list）
#   --out     <path>   CSV 输出（默认 ./bench_233.csv）
#   --limit   <int>    取列表前 N 行（默认 0=全量）
#   --threads <int>    OMP 线程数（默认 nproc）
#   --dist    <float>  Mash 距离阈值（默认 0.05）
#   --kmer    <int>    k-mer 长度（默认 21）
#   --seed    <int>    哈希种子（默认 42）
#   --repeat  <int>    每项重复次数（默认 1）
#
# 示例：
#   cd /home/user_home/tongz/Finish/RabbitSketch/examples
#   bash bench_233_index_vs_noindex.sh --threads 64
#
# bindash 的 index 模式使用 bin 位置与压缩值构建倒排索引。
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

# ── 默认值 ───────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="${SCRIPT_DIR}/rabbitsketch"
LIST="${SCRIPT_DIR}/233.list"
OUT="/home/bigssd/zt/bench_233.csv"
LIMIT=0
THREADS="$(nproc 2>/dev/null || echo 4)"
DIST="0.05"
KMER="21"
SEED="42"
REPEAT=1

# ── 参数解析 ─────────────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe)     EXE="$2";     shift 2 ;;
    --list)    LIST="$2";    shift 2 ;;
    --out)     OUT="$2";     shift 2 ;;
    --limit)   LIMIT="$2";   shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --dist)    DIST="$2";    shift 2 ;;
    --kmer)    KMER="$2";    shift 2 ;;
    --seed)    SEED="$2";    shift 2 ;;
    --repeat)  REPEAT="$2";  shift 2 ;;
    -h|--help) sed -n '3,35p' "$0"; exit 0 ;;
    *) die "未知参数: $1" ;;
  esac
done

# ── 前置检查 ─────────────────────────────────────────────────────────────────
[[ -x "$EXE" ]]  || die "可执行文件不存在或无执行权限: $EXE"
[[ -f "$LIST" ]] || die "文件列表不存在: $LIST"
[[ -x /usr/bin/time ]] || die "缺少 /usr/bin/time（GNU time，用于测 RSS）"

# list 内路径是相对路径（233/xxx），必须在 examples/ 目录下运行
cd "$SCRIPT_DIR"

# ── 准备输出目录 ──────────────────────────────────────────────────────────────
WORK_DIR="/home/bigssd/zt"
mkdir -p "$WORK_DIR"
LOG_DIR="${WORK_DIR}/bench_233_logs"
mkdir -p "$LOG_DIR"

SUB_LIST="$(readlink -f "$LIST")"
if [[ "$LIMIT" -gt 0 ]]; then
  SUB_LIST="${LOG_DIR}/sub.list"
  head -n "$LIMIT" "$(readlink -f "$LIST")" > "$SUB_LIST"
fi
ACTUAL_N="$(wc -l < "$SUB_LIST" | tr -d ' ')"

# ── 输出参数概要 ──────────────────────────────────────────────────────────────
echo "========================================"
echo "  RabbitSketch bench: --index vs allpairs"
echo "========================================"
echo "  exe     : $EXE"
echo "  list    : $LIST"
  echo "  N       : $ACTUAL_N 条序列$( [[ "$LIMIT" -gt 0 ]] && echo "（截断至 $LIMIT）" || echo "（全量）" )"
echo "  threads : $THREADS"
echo "  dist    : $DIST   kmer: $KMER   seed: $SEED"
echo "  repeat  : $REPEAT"
echo "  out csv : $OUT"
echo "  logs    : $LOG_DIR/"
echo "========================================"

# ── CSV 表头（仅在文件不存在时写；已有则续跑） ────────────────────────────────
if [[ ! -f "$OUT" ]]; then
  echo "algo,mode,n_seqs,threads,dist,kmer,seed,repeat,exit_code,elapsed_s,max_rss_mib,out_lines" > "$OUT"
  echo "  (新建 CSV)"
else
  echo "  (已有 CSV，跳过已成功的行，继续未完成项)"
fi

# 加载已成功完成的 (algo,mode,repeat) 集合
declare -A DONE
while IFS=, read -r _algo _mode _ _ _ _ _ _rep ec _rest; do
  [[ "$_algo" == "algo" ]] && continue          # 跳过表头
  [[ "$ec" == "0" ]] && DONE["${_algo}:${_mode}:${_rep}"]=1
done < "$OUT"

# ── 单次运行函数 ──────────────────────────────────────────────────────────────
run_one() {
  local algo="$1"   # e.g. hll
  local mode="$2"   # allpairs | index
  local rep="$3"

  # 已成功跑过则跳过（exit_code=0 在 CSV 中已有记录）
  if [[ "${DONE[${algo}:${mode}:${rep}]+set}" == "set" ]]; then
    printf "  %-16s %-10s rep=%d  already done, skip\n" "$algo" "$mode" "$rep"
    return 0
  fi

  local tag="${algo}_${mode}_r${rep}"
  local out_path="${WORK_DIR}/res.${tag}.tsv"   # 距离文件写到 bigssd，跑完删
  local log_path="${LOG_DIR}/run.${tag}.log"
  local time_log="${LOG_DIR}/time.${tag}.txt"

  # 构造命令
  local cmd=( "$EXE" "--${algo}" -i "$SUB_LIST" -l
              -o "$out_path" -t "$THREADS"
              -d "$DIST" -k "$KMER" --seed "$SEED" )
  [[ "$mode" == "index" ]] && cmd=( "${cmd[@]:0:2}" --index "${cmd[@]:2}" )

  printf "  %-16s %-10s rep=%d  ... " "$algo" "$mode" "$rep"

  set +e
  /usr/bin/time -f "ELAPSED_S=%e\nMAX_RSS_KIB=%M" -o "$time_log" \
    "${cmd[@]}" > "$log_path" 2>&1
  local ec=$?
  set -e

  # 从 time_log 里解析（最后一次出现的数字，避免程序自己打印干扰）
  local elapsed rss_kib rss_mib out_lines
  elapsed="$(  grep -oP 'ELAPSED_S=\K[\d.]+' "$time_log"  | tail -1 )"
  rss_kib="$(  grep -oP 'MAX_RSS_KIB=\K\d+'  "$time_log"  | tail -1 )"
  rss_mib="$(  awk -v k="${rss_kib:-0}" 'BEGIN{printf "%.1f", k/1024}' )"
  out_lines="$( [[ -f "$out_path" ]] && wc -l < "$out_path" | tr -d ' ' || echo 0 )"

  printf "done  elapsed=%-7s s   RSS=%-7s MiB   exit=%s\n" \
    "${elapsed:--}" "${rss_mib:--}" "$ec"

  # 写 CSV
  printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
    "$algo" "$mode" "$ACTUAL_N" "$THREADS" "$DIST" "$KMER" \
    "$SEED" "$rep" "$ec" "${elapsed:-}" "${rss_mib:-}" "$out_lines" \
    >> "$OUT"

  # 距离文件空间占用大，统计完立即删除
  [[ -f "$out_path" ]] && rm -f "$out_path"
}

# ── 算法列表（与 main.cpp 一致） ──────────────────────────────────────────────
# 第一列=算法名，第二列=是否支持 --index（1=支持）
declare -A SUPPORTS_INDEX=(
  [minhash]=1
  [kssd]=1
  [hll]=1
  [probminhash]=1
  [bindash]=1
  [setsketch]=1
  [fastkmv]=1
)
ALGO_ORDER=( minhash kssd hll probminhash bindash setsketch fastkmv )

# ── 主循环 ────────────────────────────────────────────────────────────────────
for algo in "${ALGO_ORDER[@]}"; do
  for (( r=1; r<=REPEAT; r++ )); do

    echo ""
    echo "── $algo (repeat $r) ──"

    run_one "$algo" "allpairs" "$r"

    if [[ "${SUPPORTS_INDEX[$algo]}" -eq 1 ]]; then
      run_one "$algo" "index" "$r"
    else
      printf "  %-16s %-10s rep=%d  skipped (not supported)\n" "$algo" "index" "$r"
      printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,,,%s\n" \
        "$algo" "index" "$ACTUAL_N" "$THREADS" "$DIST" "$KMER" "$SEED" "$r" "NA" "" \
        >> "$OUT"
    fi

  done
done

# ── 结尾汇总表 ────────────────────────────────────────────────────────────────
echo ""
echo "========================================"
echo "  汇总（写入 $OUT）"
echo "========================================"
printf "  %-14s %-10s %-10s %-10s\n" "algo" "mode" "elapsed_s" "rss_mib"
echo "  ─────────────────────────────────────"
tail -n +2 "$OUT" | awk -F, '{printf "  %-14s %-10s %-10s %-10s\n",$1,$2,$10,$11}'
echo "========================================"
echo "  详细日志: $LOG_DIR/"
