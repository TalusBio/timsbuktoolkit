#!/usr/bin/env bash
# Same configs as the full-run sweeps, on a 2-minute RT slice. If the slice
# ranks them the way the full run did, iteration can happen on the slice.
set -u
cd "$(dirname "$0")/.."
RAW=~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d
FASTA=~/fasta/hela_gt20peps.fasta
LIB=shitshit/hela_gt20peps.mzspeclib.txt.gz
BIN=./target/release/timsseek
OUT=shitshit/slice_results.tsv
SLICE="${SLICE:-[600.0, 720.0]}"

printf "tag\tfull_run_q01\tstatus\tpeak_fp_mb\twall_s\tidx_s\tq01\tq05\tno_signal\n" > "$OUT"

run_point() {
  local tag=$1 full=$2 m1=$3 m2=$4
  local cfg=shitshit/cfg_slice_${tag}.toml log=/tmp/slice_${tag}.log res=shitshit/res_slice_${tag}
  $BIN --print-default-config > "$cfg"
  printf '\n[indexing]\nrt_range_seconds = %s\n' "$SLICE" >> "$cfg"
  [ -n "$m1" ] && printf '\n[indexing.ms1]\n%s\n' "$(echo -e "$m1")" >> "$cfg"
  [ -n "$m2" ] && printf '\n[indexing.ms2]\n%s\n' "$(echo -e "$m2")" >> "$cfg"
  local start; start=$(date +%s)
  local libflag=(--speclib-uri "$LIB")
  [ -f "$LIB" ] || libflag=(--speclib-uri "$LIB" --build-if-missing)
  $BIN --config "$cfg" --fasta $FASTA "${libflag[@]}" --output-dir "$res" --raw-inputs $RAW -O > "$log" 2>&1 &
  local pid=$! peak=0 status=ok fp
  while kill -0 $pid 2>/dev/null; do
    fp=$(footprint -p $pid 2>/dev/null | grep -o "Footprint: [0-9]* MB" | grep -o "[0-9]*")
    [ -n "$fp" ] && [ "$fp" -gt "$peak" ] && peak=$fp
    sleep 2
  done
  wait $pid 2>/dev/null; local rc=$?
  [ $rc -ne 0 ] && status="exit_$rc"
  local wall=$(( $(date +%s) - start )) idx nosig q01="" q05=""
  idx=$(grep -o "Loading index [. ]* [0-9.]*s" "$log" | grep -o "[0-9.]*s$" | tr -d s)
  nosig=$(grep -o "no_observed_signal=[0-9]*" "$log" | head -1 | cut -d= -f2)
  if [ "$status" = ok ]; then
    local parquet; parquet=$(ls "$res"/*/results.parquet 2>/dev/null | head -1)
    [ -n "$parquet" ] && read -r q01 q05 < <(uv run python -c "
import polars as pl
t = pl.read_parquet('$parquet').filter(pl.col('is_target'))
print(t.filter(pl.col('qvalue')<=0.01).height, t.filter(pl.col('qvalue')<=0.05).height)
" 2>/dev/null)
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$full" "$status" "$peak" "$wall" "${idx:-}" "$q01" "$q05" "${nosig:-}" >> "$OUT"
}

CAP='window_cap = { max_peaks = 500, window_da = 50.0 }'
run_point baseline      4036 "" ""
run_point cap500_50     4481 "" "$CAP"
run_point cap_nt        4544 "transitive = false" "$CAP\ntransitive = false"
run_point ms2_20k       4128 "" "max_peaks = 20000"
run_point cap750_100    4268 "" 'window_cap = { max_peaks = 750, window_da = 100.0 }'
run_point cap_bins1     4101 "" "$CAP\nmz_tol = { Bins = 1 }"
run_point ms2_10k       3880 "" "max_peaks = 10000"
echo SLICE_DONE >> "$OUT"
