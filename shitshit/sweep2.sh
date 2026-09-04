#!/usr/bin/env bash
# Window-cap sweep, MS2 only. MS1 stays at the library default so the effect of
# the MS2 cap is read on its own.
set -u
cd "$(dirname "$0")/.."
RAW=~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d
FASTA=~/fasta/hela_gt20peps.fasta
BIN=./target/release/timsseek
OUT=shitshit/sweep2_results.tsv

printf "tag\tN\tW\tpeak_gb\twall_s\tidx_s\tq01\tq05\tno_signal\n" > "$OUT"

run_point() {
  local n=$1 w=$2 tag="ms2only_w${2}_n${1}"
  local cfg=shitshit/cfg_${tag}.toml log=/tmp/sweep_${tag}.log res=shitshit/res_${tag}
  $BIN --print-default-config > "$cfg"
  printf '\n[indexing.ms2]\nwindow_cap = { max_peaks = %s, window_da = %s }\n' "$n" "$w" >> "$cfg"
  /usr/bin/time -l $BIN --config "$cfg" --fasta $FASTA --output-dir "$res" --raw-inputs $RAW -O > "$log" 2>&1
  local peak wall idx q01 q05 nosig
  peak=$(awk '/peak memory footprint/ {printf "%.2f", $1/1e9}' "$log")
  wall=$(awk '/ real / {printf "%.0f", $1}' "$log")
  idx=$(grep -o "Loading index [. ]* [0-9.]*s" "$log" | grep -o "[0-9.]*s$" | tr -d s)
  nosig=$(grep -o "no_observed_signal=[0-9]*" "$log" | head -1 | cut -d= -f2)
  local parquet; parquet=$(ls "$res"/*/results.parquet | head -1)
  read -r q01 q05 < <(uv run python -c "
import polars as pl
t = pl.read_parquet('$parquet').filter(pl.col('is_target'))
print(t.filter(pl.col('qvalue')<=0.01).height, t.filter(pl.col('qvalue')<=0.05).height)
" 2>/dev/null)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$n" "$w" "$peak" "$wall" "$idx" "$q01" "$q05" "$nosig" >> "$OUT"
}

run_point 1000 100
run_point 2000 100
run_point 1500 100
echo SWEEP_DONE >> "$OUT"
