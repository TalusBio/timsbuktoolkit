#!/usr/bin/env bash
# Mobility ablation: same centroiding configs as before, with the mobility
# axis treated as absent in every search phase. If the ID ranking survives,
# the centroiding effect is on m/z position (or count), not on the search's
# mobility window fitting itself to the centroids.
set -u
cd "$(dirname "$0")/.."
RAW=~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d
LIB=shitshit/hela_gt20peps.mzspeclib.txt.gz
BIN=./target/release/timsseek
OUT=shitshit/sweep6_results.tsv
MEM_LIMIT_MB=8192
printf "tag\twith_mobility_q01\tstatus\tpeak_fp_mb\twall_s\tidx_s\tq01\tq05\tno_signal\n" > "$OUT"

run_point() {
  local tag=$1 ref=$2 m1=$3 m2=$4
  local cfg=shitshit/cfg_${tag}.toml log=/tmp/sweep_${tag}.log res=shitshit/res_${tag}
  $BIN --print-default-config > "$cfg"
  printf '\n[indexing]\nignore_mobility = true\n' >> "$cfg"
  printf '\n[indexing.ms1]\n%s\n' "$(echo -e "$m1")" >> "$cfg"
  printf '\n[indexing.ms2]\nwindow_cap = { max_peaks = 500, window_da = 50.0 }\n%s\n' "$(echo -e "$m2")" >> "$cfg"
  local start; start=$(date +%s)
  $BIN --config "$cfg" --speclib-uri "$LIB" --output-dir "$res" --raw-inputs $RAW -O > "$log" 2>&1 &
  local pid=$! peak=0 status=ok fp
  while kill -0 $pid 2>/dev/null; do
    fp=$(footprint -p $pid 2>/dev/null | grep -o "Footprint: [0-9]* MB" | grep -o "[0-9]*")
    if [ -n "$fp" ]; then
      [ "$fp" -gt "$peak" ] && peak=$fp
      [ "$fp" -gt "$MEM_LIMIT_MB" ] && { kill -9 $pid; status=OOM_KILLED; }
    fi
    sleep 5
  done
  wait $pid 2>/dev/null; local rc=$?
  [ "$status" = ok ] && [ $rc -ne 0 ] && status="exit_$rc"
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
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$ref" "$status" "$peak" "$wall" "${idx:-}" "$q01" "$q05" "${nosig:-}" >> "$OUT"
}

run_point nomob_ref          4481 "" ""
run_point nomob_nt_both_im1  4624 "transitive = false\nim_tol = { Pct = 1.0 }" "transitive = false\nim_tol = { Pct = 1.0 }"
run_point nomob_bins1        4101 "mz_tol = { Bins = 1 }" "mz_tol = { Bins = 1 }"
run_point nomob_nt_bins1     4269 "transitive = false\nmz_tol = { Bins = 1 }" "transitive = false\nmz_tol = { Bins = 1 }"
run_point nomob_im4_5        ""   "im_tol = { Pct = 4.0 }" "im_tol = { Pct = 5.0 }"
echo SWEEP_DONE >> "$OUT"
