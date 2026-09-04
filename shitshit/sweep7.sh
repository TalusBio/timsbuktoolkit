#!/usr/bin/env bash
# Phase 3 window test. The derived windows are ~60 ppm / 14 % because the
# calibrant residual MAD is inflated (13.6 ppm, 3.2 %). Shrinking the sigma
# multipliers is a blunt way to run Phase 3 at ~10 ppm / ~2.4 % and see what
# the loose windows cost. MS2 cap 500/50, default centroiding.
set -u
cd "$(dirname "$0")/.."
RAW=~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d
LIB=shitshit/hela_gt20peps.mzspeclib.txt.gz
BIN=./target/release/timsseek
OUT=shitshit/sweep7_results.tsv
MEM_LIMIT_MB=8192
printf "tag\tstatus\tpeak_fp_mb\twall_s\tphase3_s\tmz_win\tmob_win\tq01\tq05\tno_signal\n" > "$OUT"

run_point() {
  local tag=$1 calib=$2
  local cfg=shitshit/cfg_${tag}.toml log=/tmp/sweep_${tag}.log res=shitshit/res_${tag}
  $BIN --print-default-config > "$cfg"
  # [calibration] already exists in the template with n_calibrants etc.; TOML
  # forbids redefining a table, so patch the sigma lines in place.
  [ -n "$calib" ] && while IFS='=' read -r k v; do
    [ -z "$k" ] && continue
    sed -i '' "s/^${k} *=.*/${k} =${v}/" "$cfg"
  done <<< "$(echo -e "$calib")"
  printf '\n[indexing.ms2]\nwindow_cap = { max_peaks = 500, window_da = 50.0 }\n' >> "$cfg"
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
  local wall=$(( $(date +%s) - start )) p3 mzw mobw nosig q01="" q05=""
  local clean; clean=$(sed 's/\x1b\[[0-9;]*m//g' "$log")
  mzw=$(echo "$clean" | grep -o "m/z: ([0-9.]*, [0-9.]*) ppm" | head -1 | grep -o "([^)]*)")
  mobw=$(echo "$clean" | grep -o "mobility: ([0-9.]*, [0-9.]*) %" | head -1 | grep -o "([^)]*)")
  p3=$(echo "$clean" | grep -A1 "^Phase 3" | grep -o "^ *[0-9.]*s " | tr -d ' s')
  nosig=$(echo "$clean" | grep -o "no_observed_signal=[0-9]*" | head -1 | cut -d= -f2)
  if [ "$status" = ok ]; then
    local parquet; parquet=$(ls "$res"/*/results.parquet 2>/dev/null | head -1)
    [ -n "$parquet" ] && read -r q01 q05 < <(uv run python -c "
import polars as pl
t = pl.read_parquet('$parquet').filter(pl.col('is_target'))
print(t.filter(pl.col('qvalue')<=0.01).height, t.filter(pl.col('qvalue')<=0.05).height)
" 2>/dev/null)
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$status" "$peak" "$wall" "${p3:-}" "${mzw:-}" "${mobw:-}" "$q01" "$q05" "${nosig:-}" >> "$OUT"
}

run_point sig3_3      ""
run_point sig0p5_3    "mz_sigma= 0.5"
run_point sig3_0p5    "mobility_sigma= 0.5"
run_point sig0p5_0p5  "mz_sigma= 0.5\nmobility_sigma= 0.5"
run_point sig1_1      "mz_sigma= 1.0\nmobility_sigma= 1.0"
echo SWEEP_DONE >> "$OUT"
