#!/bin/bash

RAW_FILE=$1
if [ -z "$RAW_FILE" ]; then
    echo "Please provide a raw file"
    exit 1
fi

cargo run --bin timsseek_rts --release -- \
    --config ./tolconfig.json \
    --dotd-file $RAW_FILE &
SERVER_PID=$!

uv run --verbose python -m streamlit run python/timsseek_rts_receiver/timsseek_rts_receiver/receiver.py
kill $SERVER_PID
wait
