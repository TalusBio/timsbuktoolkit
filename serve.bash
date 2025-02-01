#!/bin/bash

RAW_FILE=$1

cargo run --bin timsseek_rts --release -- \
    --config ./tolconfig.json \
    --dotd-file $RAW_FILE &
SERVER_PID=$!

uv run --project timsseek_rts/python/ --verbose streamlit run timsseek_rts/python/receiver.py
kill $SERVER_PID
wait
