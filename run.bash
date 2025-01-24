
FASTA_FILE="asdasd.fasta"
SPECLIB_NAME="asdad.ndjson"
DOTD_FILE="/Users/sebastianpaez/git/ionmesh/benchmark/240402_PRTC_01_S1-A1_1_11342.d"
RESULTS_DIR="full_human_predicted"

uv run speclib_build_fasta \
    --fasta_file $FASTA_FILE \
    --decoy_strategy REVERSE \
    --max_ions 10 \
    --outfile $SPECLIB_NAME \
    --model onnx

cargo run --release --bin timsseek -- \
    --config $SPECLIB_NAME \
    --output-dir $RESULTS_DIR \
    --dotd-file $DOTD_FILE
    
uv run -s showscores.py --results_dir $RESULTS_DIR