#!/usr/bin/env bash
# Does non-transitive clustering change IDs? MS2 window cap fixed at 500/50.
# Points are extra TOML lines appended under [indexing.ms1] / [indexing.ms2].
set -u
cd "$(dirname "$0")/.."
RAW=~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d
FASTA=~/fasta/hela_gt20peps.fasta
BIN=./target/release/timsseek
OUT=shitshit/sweep4_results.tsv
MEM_LIMIT_MB=8192

printf "tag\tms1_extra\tms2_extra\tstatus\tpeak_fp_mb\twall_s\tidx_s\tq01\tq05\tno_signal\n" > "$OUT"

run_point() {
  local tag=$1 m1=$2 m2=$3
  local cfg=shitshit/cfg_${tag}.toml log=/tmp/sweep_${tag}.log res=shitshit/res_${tag}
  $BIN --print-default-config > "$cfg"
  printf '\n[indexing.ms1]\n%s\n' "$(echo -e "$m1")" >> "$cfg"
  printf '\n[indexing.ms2]\nwindow_cap = { max_peaks = 500, window_da = 50.0 }\n%s\n' "$(echo -e "$m2")" >> "$cfg"
  local start; start=$(date +%s)
  $BIN --config "$cfg" --fasta $FASTA --output-dir "$res" --raw-inputs $RAW -O > "$log" 2>&1 &
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
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$(echo "$m1" | tr '\n' ';')" "$(echo "$m2" | tr '\n' ';')" "$status" "$peak" "$wall" "${idx:-}" "$q01" "$q05" "${nosig:-}" >> "$OUT"
}

run_point nt_ppm           "transitive = false"                               "transitive = false"
run_point nt_bins1         "transitive = false\nmz_tol = { Bins = 1 }"        "transitive = false\nmz_tol = { Bins = 1 }"
run_point nt_bins1_im4     "transitive = false\nmz_tol = { Bins = 1 }\nim_tol = { Pct = 4.0 }" "transitive = false\nmz_tol = { Bins = 1 }\nim_tol = { Pct = 5.0 }"
run_point nt_ppm_im4       "transitive = false\nim_tol = { Pct = 4.0 }"       "transitive = false\nim_tol = { Pct = 5.0 }"
echo SWEEP_DONE >> "$OUT"
