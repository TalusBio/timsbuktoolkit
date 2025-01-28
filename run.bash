
set -x
set -e

FASTA_FILE="/Users/sebastianpaez/fasta/hela_gt20peps.fasta"
SPECLIB_NAME="asdad.ndjson"
DOTD_FILE="/Users/sebastianpaez/git/ionmesh/benchmark/240402_PRTC_01_S1-A1_1_11342.d"
RESULTS_DIR="top_proteins_hela"
SUMMARY_DIR="top_proteins_hela_summary"

cat << EOF > config_use.json
{
    "analysis": {
        "chunk_size": 20000,
        "tolerance": {
            "ms": {"ppm":  [15.0, 15.0]},
            "mobility": {"percent": [3.0, 3.0]},
            "quad": {"absolute": [0.1, 0.1]}
        }
    }
}
EOF

# uv run speclib_build_fasta \
#     --fasta_file $FASTA_FILE \
#     --decoy_strategy REVERSE \
#     --max_ions 10 \
#     --outfile $SPECLIB_NAME \
#     --model onnx

cargo run --release --bin timsseek -- \
    --config config_use.json \
    --speclib-file $SPECLIB_NAME \
    --output-dir $RESULTS_DIR \
    --dotd-file $DOTD_FILE

# Technically does T/D competition
uv run -s showscores.py --results_dir $RESULTS_DIR --output_dir $SUMMARY_DIR