#!/bin/bash

set -x
set -e
set -o pipefail

# This generates the `ubb_elution_groups.json` file.
# RAW_FILE=../230510_PRTC_13_S1-B1_1_12817.d/
# ./sage --annotate-matches --fasta ~/fasta/20231030_LINEARIZED_UP000005640_9606.fasta check_config.json $RAW_FILE
# uv run build.py

../../target/release/timsquery query-index --raw-file-path ../230510_PRTC_13_S1-B1_1_12817.d --tolerance-settings-path "../../templates/tolerance_settings_narrow_rt.json" --elution-groups-path "ubb_elution_groups.json" --output-path "query_results" --index expanded-raw-frame-index --aggregator multi-cmg-stats --format pretty-json
uv run plot.py
