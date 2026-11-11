#!/usr/bin/env bash
# Thread-scaling benchmark for RabbitSketch --index on 233.list.
#
# Defaults are intentionally minimal: all algorithm-specific parameters and the
# distance threshold use rabbitsketch's built-in defaults. BinDash --index
# uses the native (bin position, b-bit value) inverted index.
#
# Usage:
#   bash bench_233_index_scaling.sh
#   bash bench_233_index_scaling.sh --threads "1 2 4 8 16 32 64"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="${SCRIPT_DIR}/rabbitsketch"
LIST="${SCRIPT_DIR}/233.list"
WORK_DIR="/home/bigssd/zt/rabbitsketch_index_scaling"
OUT="${WORK_DIR}/bench_233_index_scaling.csv"
LOG_DIR="${WORK_DIR}/logs"
THREADS_LIST="1 2 4 8 16 32 64"

die() { echo "ERROR: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe)     EXE="$2"; shift 2 ;;
    --list)    LIST="$2"; shift 2 ;;
    --out)     OUT="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; LOG_DIR="${2}/logs"; shift 2 ;;
    --threads) THREADS_LIST="$2"; shift 2 ;;
    -h|--help)
      sed -n '3,14p' "$0"
      exit 0
      ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -x "$EXE" ]] || die "missing executable: $EXE"
[[ -f "$LIST" ]] || die "missing list: $LIST"
[[ -x /usr/bin/time ]] || die "missing /usr/bin/time"

mkdir -p "$WORK_DIR" "$LOG_DIR"
cd "$SCRIPT_DIR"

N_SEQS="$(wc -l < "$LIST" | tr -d ' ')"

if [[ ! -f "$OUT" ]]; then
  echo "algo,index_supported,n_seqs,threads,exit_code,elapsed_s,max_rss_mib,out_lines,cmd" > "$OUT"
fi

declare -A DONE
while IFS=, read -r algo _support _n threads ec _rest; do
  [[ "$algo" == "algo" ]] && continue
  [[ "$ec" == "0" ]] && DONE["${algo}:${threads}"]=1
done < "$OUT"

ALGO_ORDER=(fastkmv probminhash kssd setsketch hll bindash)
declare -A SUPPORTS_INDEX=(
  [minhash]=1
  [fastkmv]=1
  [probminhash]=1
  [kssd]=1
  [setsketch]=1
  [hll]=1
  [bindash]=1
)

echo "RabbitSketch --index thread scaling"
echo "  exe      : $EXE"
echo "  list     : $LIST"
echo "  n_seqs   : $N_SEQS"
echo "  threads  : $THREADS_LIST"
echo "  output   : $OUT"
echo "  logs     : $LOG_DIR"
echo

run_one() {
  local algo="$1"
  local threads="$2"
  local support="${SUPPORTS_INDEX[$algo]}"

  if [[ "${DONE[${algo}:${threads}]+set}" == "set" ]]; then
    printf "skip %-12s T=%-3s (already recorded)\n" "$algo" "$threads"
    return 0
  fi

  if [[ "$support" -eq 0 ]]; then
    printf "skip %-12s T=%-3s (--index unsupported)\n" "$algo" "$threads"
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%q\n" \
      "$algo" "$support" "$N_SEQS" "$threads" "NA" "" "" "" \
      "$EXE --$algo --index -i $LIST -l -t $threads" >> "$OUT"
    return 0
  fi

  local tag="${algo}_T${threads}"
  local dist_out="${WORK_DIR}/${tag}.dist.tsv"
  local log_path="${LOG_DIR}/${tag}.log"
  local time_path="${LOG_DIR}/${tag}.time"
  local cmd=( "$EXE" "--${algo}" --index -i "$LIST" -l -t "$threads" -o "$dist_out" )

  printf "run  %-12s T=%-3s ... " "$algo" "$threads"
  set +e
  /usr/bin/time -f "ELAPSED_S=%e\nMAX_RSS_KIB=%M" -o "$time_path" \
    "${cmd[@]}" > "$log_path" 2>&1
  local ec=$?
  set -e

  local elapsed rss_kib rss_mib out_lines
  elapsed="$(awk -F= '$1=="ELAPSED_S"{v=$2} END{print v}' "$time_path" 2>/dev/null || true)"
  rss_kib="$(awk -F= '$1=="MAX_RSS_KIB"{v=$2} END{print v}' "$time_path" 2>/dev/null || true)"
  rss_mib="$(awk -v k="${rss_kib:-0}" 'BEGIN{printf "%.1f", k/1024}')"
  out_lines="$( [[ -f "$dist_out" ]] && wc -l < "$dist_out" | tr -d ' ' || echo 0 )"

  printf "exit=%s elapsed=%ss RSS=%sMiB lines=%s\n" \
    "$ec" "${elapsed:-}" "$rss_mib" "$out_lines"

  printf "%s,%s,%s,%s,%s,%s,%s,%s,%q\n" \
    "$algo" "$support" "$N_SEQS" "$threads" "$ec" "${elapsed:-}" \
    "$rss_mib" "$out_lines" "${cmd[*]}" >> "$OUT"

  rm -f "$dist_out"
}

for algo in "${ALGO_ORDER[@]}"; do
  echo "== $algo =="
  for threads in $THREADS_LIST; do
    run_one "$algo" "$threads"
  done
done

echo
echo "Done: $OUT"
