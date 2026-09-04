#!/usr/bin/env bash
# Centroiding config sweep on the single-window isotope fixture.
# One line per (window, config): MS1 precursor envelope metrics, MS2 fragment
# envelope metrics, null rates, centroiding time. Seconds per config.
set -u
cd "$(dirname "$0")/.."
RAW=~/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_off_DIA.d
SEEDS=shitshit/seeds_all.csv
P=./target/release/examples/isotope_frames
OUT=shitshit/probe_sweep.tsv
CENTERS="${CENTERS:-450 660}"
WINDOW="${WINDOW:-15}"

printf "center\tconfig\tseeds\tframes\tms\tms1_vis\tms1_cpa\tms1_split\tms1_null\tms1_cos50\tms1_cos95\tms2_frags\tms2_cpa\tms2_split\tms2_null\tms2_cos50\tms2_cos95\n" > "$OUT"

run() {
  local tag=$1; shift
  for c in $CENTERS; do
    local log=/tmp/probe_${tag}_${c}.log
    $P $RAW $SEEDS --rt-center "$c" --rt-window "$WINDOW" --halo 1 "$@" > "$log" 2>&1 || { echo "FAILED $tag $c"; continue; }
    awk -v center="$c" -v tag="$tag" '
      /^seeds:/        { seeds=$2 }
      /^centroided/    { frames=$2; ms=$5 }
      /^##### MS1/     { sec="ms1" }
      /^##### MS2/     { sec="ms2"; frags=$5 }
      /^-- isotope positions/ { blk="on" }
      /^-- null positions/    { blk="null" }
      /^ *M\+0 / && blk=="on"   { if (sec=="ms1") { v1=$2; c1=$4; s1=$5 } else { c2=$4; s2=$5 } }
      /^ *M\+0 / && blk=="null" { if (sec=="ms1") n1=$2; else n2=$2 }
      /^averagine cosine/ {
        match($0, /p50 [0-9.]+/); p50=substr($0, RSTART+4, RLENGTH-4)
        match($0, />=0.95 [0-9.]+/); p95=substr($0, RSTART+7, RLENGTH-7)
        if (sec=="ms1") { cos1=p50; q1=p95 } else { cos2=p50; q2=p95 }
      }
      END { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            center, tag, seeds, frames, ms, v1, c1, s1, n1, cos1, q1, frags, c2, s2, n2, cos2, q2 }
    ' "$log" >> "$OUT"
  done
}

CAP=(--ms2-cap 500/50)
run ref            "${CAP[@]}"
run nocap
run cap1000_100    --ms2-cap 1000/100
run cap250_25      --ms2-cap 250/25
run nt             "${CAP[@]}" --transitive false
run bins1          "${CAP[@]}" --ms1-mz-bins 1 --ms2-mz-bins 1
run bins2          "${CAP[@]}" --ms1-mz-bins 2 --ms2-mz-bins 2
run nt_bins1       "${CAP[@]}" --transitive false --ms1-mz-bins 1 --ms2-mz-bins 1
run nt_bins2       "${CAP[@]}" --transitive false --ms1-mz-bins 2 --ms2-mz-bins 2
run nt_bins3       "${CAP[@]}" --transitive false --ms1-mz-bins 3 --ms2-mz-bins 3
run nt_bins2_im3   "${CAP[@]}" --transitive false --ms1-mz-bins 2 --ms2-mz-bins 2 --ms1-im-pct 3
run ms1_im1        "${CAP[@]}" --ms1-im-pct 1
run ms1_im3        "${CAP[@]}" --ms1-im-pct 3
run ms1_im4        "${CAP[@]}" --ms1-im-pct 4
run ms2_im2        "${CAP[@]}" --ms2-im-pct 2
run ms2_im5        "${CAP[@]}" --ms2-im-pct 5
run nt_ms1_im3     "${CAP[@]}" --transitive false --ms1-im-pct 3
run nt_ms1_im4     "${CAP[@]}" --transitive false --ms1-im-pct 4
run ms1_im0p5      "${CAP[@]}" --ms1-im-pct 0.5
run nt_ms1_im1     "${CAP[@]}" --transitive false --ms1-im-pct 1
run ms2_im1        "${CAP[@]}" --ms2-im-pct 1
run both_im1       "${CAP[@]}" --ms1-im-pct 1 --ms2-im-pct 1
run nt_both_im1    "${CAP[@]}" --transitive false --ms1-im-pct 1 --ms2-im-pct 1
run ms1_es0        "${CAP[@]}" --ms1-early-stop 0
run ms1_cap60k     "${CAP[@]}" --ms1-max-peaks 60000
run ms1_cap10k     "${CAP[@]}" --ms1-max-peaks 10000
echo DONE >> "$OUT"
