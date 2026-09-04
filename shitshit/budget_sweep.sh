#!/usr/bin/env bash
# Window-cap sweep under a wall-clock budget and a memory watchdog.
#
# Each point sets `window_cap` on MS1 and/or MS2 (`-` = library default) and
# runs a full timsseek search. A watchdog polls `footprint -p` (physical
# footprint: resident + compressed + swapped) every 2 s and kills the run past
# MEM_LIMIT_MB. Points run in priority order until the budget cannot fit
# another one.
set -u
cd "$(dirname "$0")/.."
RAW=~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d
FASTA=~/fasta/hela_gt20peps.fasta
BIN=./target/release/timsseek
OUT=shitshit/budget_sweep_results.tsv
BUDGET_S=3600
MEM_LIMIT_MB=8192
POINT_ESTIMATE_S=200   # stop when less than this remains

printf "tag\tms1_cap\tms2_cap\tstatus\tpeak_fp_mb\twall_s\tidx_s\tq01\tq05\tno_signal\n" > "$OUT"
T0=$(date +%s)

# tag  ms1_N ms1_W  ms2_N ms2_W
POINTS=(
  "ms2_1000_100    -    -    1000 100"
  "ms2_2000_100    -    -    2000 100"
  "ms2_3000_100    -    -    3000 100"
  "ms2_1500_100    -    -    1500 100"
  "ms1_1000_100    1000 100  -    -"
  "ms1_2000_100    2000 100  -    -"
  "ms1_3000_100    3000 100  -    -"
  "ms2_1000_50     -    -    1000 50"
  "ms2_4000_200    -    -    4000 200"
  "ms2_500_50      -    -    500  50"
  "ms2_250_25      -    -    250  25"
  "both_3000_2000  3000 100  2000 100"
  "both_2000_1500  2000 100  1500 100"
  "ms2_750_100     -    -    750  100"
  "ms1_1500_100    1500 100  -    -"
  "ms2_4000_100    -    -    4000 100"
  "both_3000_3000  3000 100  3000 100"
  "ms2_2000_50     -    -    2000 50"
  "ms2_8000_200    -    -    8000 200"
  "ms1_4000_100    4000 100  -    -"
)

run_point() {
  local tag=$1 m1n=$2 m1w=$3 m2n=$4 m2w=$5
  local cfg=shitshit/cfg_${tag}.toml log=/tmp/sweep_${tag}.log res=shitshit/res_${tag}
  $BIN --print-default-config > "$cfg"
  [ "$m1n" != "-" ] && printf '\n[indexing.ms1]\nwindow_cap = { max_peaks = %s, window_da = %s }\n' "$m1n" "$m1w" >> "$cfg"
  [ "$m2n" != "-" ] && printf '\n[indexing.ms2]\nwindow_cap = { max_peaks = %s, window_da = %s }\n' "$m2n" "$m2w" >> "$cfg"

  local start; start=$(date +%s)
  $BIN --config "$cfg" --fasta $FASTA --output-dir "$res" --raw-inputs $RAW -O > "$log" 2>&1 &
  local pid=$!
  local peak=0 status=ok fp
  while kill -0 $pid 2>/dev/null; do
    fp=$(footprint -p $pid 2>/dev/null | grep -o "Footprint: [0-9]* MB" | grep -o "[0-9]*")
    if [ -n "$fp" ]; then
      [ "$fp" -gt "$peak" ] && peak=$fp
      if [ "$fp" -gt "$MEM_LIMIT_MB" ]; then
        kill -9 $pid 2>/dev/null
        status=OOM_KILLED
        echo "WATCHDOG: footprint ${fp} MB > ${MEM_LIMIT_MB} MB, killed" >> "$log"
      fi
    fi
    sleep 2
  done
  wait $pid 2>/dev/null; local rc=$?
  [ "$status" = ok ] && [ $rc -ne 0 ] && status="exit_$rc"
  local wall=$(( $(date +%s) - start ))

  local idx q01 q05 nosig
  idx=$(grep -o "Loading index [. ]* [0-9.]*s" "$log" | grep -o "[0-9.]*s$" | tr -d s)
  nosig=$(grep -o "no_observed_signal=[0-9]*" "$log" | head -1 | cut -d= -f2)
  q01=""; q05=""
  if [ "$status" = ok ]; then
    local parquet; parquet=$(ls "$res"/*/results.parquet 2>/dev/null | head -1)
    [ -n "$parquet" ] && read -r q01 q05 < <(uv run python -c "
import polars as pl
t = pl.read_parquet('$parquet').filter(pl.col('is_target'))
print(t.filter(pl.col('qvalue')<=0.01).height, t.filter(pl.col('qvalue')<=0.05).height)
" 2>/dev/null)
  fi
  local c1="-" c2="-"
  [ "$m1n" != "-" ] && c1="${m1n}/${m1w}"
  [ "$m2n" != "-" ] && c2="${m2n}/${m2w}"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$c1" "$c2" "$status" "$peak" "$wall" "${idx:-}" "${q01:-}" "${q05:-}" "${nosig:-}" >> "$OUT"
  # Results parquet is small; drop the rest of the run dir to save disk.
  find "$res" -name "*.log*" -delete 2>/dev/null
}

for p in "${POINTS[@]}"; do
  elapsed=$(( $(date +%s) - T0 ))
  if [ $(( BUDGET_S - elapsed )) -lt $POINT_ESTIMATE_S ]; then
    echo "budget exhausted after ${elapsed}s" >> "$OUT"
    break
  fi
  # shellcheck disable=SC2086
  run_point $p
done
echo "SWEEP_DONE total_s=$(( $(date +%s) - T0 ))" >> "$OUT"
