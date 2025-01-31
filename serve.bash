#!/bin/bash

cargo run --bin timsseek_rts --release -- \
    --config ./tolconfig.json \
    --dotd-file /Users/sebastianpaez/git/2025_dev_engine_deatchmatch/raw_data/MSR2963_SET5REP2D7_DMSO_DIA_S4-D7_1_7173.d &
pid1=$!
trap "kill -2 $pid1" SIGINT

uv run --project timsseek_rts/python/ --verbose streamlit run timsseek_rts/python/receiver.py
wait
