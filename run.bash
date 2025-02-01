
set -x
set -e

if [ -n "${FULL_RUN}" ]; then
    # Full run
    echo "Full run"
    sleep 2
    FASTA_FILE="$HOME/fasta/20231030_LINEARIZED_UP000005640_9606.fasta"
    SPECLIB_NAME="20231030_LINEARIZED_UP000005640_9606.ndjson"
    DOTD_FILE="/Users/sebastianpaez/git/ionmesh/benchmark/240402_PRTC_01_S1-A1_1_11342.d"
    RESULTS_DIR="hela_search_results"
    SUMMARY_DIR="hela_search_summary"
    EXTRAS=""
elif [ -n "${FULL_MCCOSS}" ]; then
    echo "Bo data run"
    sleep 2
    DOTD_FILE="$HOME/data/bo_maccoss/N20211212chenc_WOSP00101_DIA_60min_K562_rep1_1_Slot2-37_1_9898.d"
    FASTA_FILE="$HOME/fasta/20231030_LINEARIZED_UP000005640_9606.fasta"
    SPECLIB_NAME="20231030_LINEARIZED_UP000005640_9606.ndjson"
    RESULTS_DIR="mccoss_search_results"
    SUMMARY_DIR="mccoss_search_summary"
    EXTRAS=""
elif [ -n "${VIMENTIN_ONLY}" ]; then
    echo "VIM only"
    sleep 2
    DOTD_FILE="$HOME/git/ionmesh/benchmark/240402_PRTC_01_S1-A1_1_11342.d"
    FASTA_FILE="$HOME/fasta/VIMENTIN.fasta"
    SPECLIB_NAME="vimentin.ndjson"
    RESULTS_DIR="vimentin_search_results"
    SUMMARY_DIR="vimentin_search_summary"
    EXTRAS="--full-output"
else
    # Quick run
    echo "Quick run"
    sleep 2
    FASTA_FILE="$HOME/fasta/hela_gt20peps.fasta"
    SPECLIB_NAME="asdad.ndjson"
    DOTD_FILE="$HOME/git/ionmesh/benchmark/240402_PRTC_01_S1-A1_1_11342.d"
    RESULTS_DIR="top_proteins_hela"
    SUMMARY_DIR="top_proteins_hela_summary"
    EXTRAS=""
fi

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

if [ -f "$SPECLIB_NAME" ]; then
    sleep 2
    echo "Speclib already built"
else
    echo "Building speclib"
    sleep 2
    uv run speclib_build_fasta \
        --fasta_file $FASTA_FILE \
        --decoy_strategy REVERSE \
        --max_ions 10 \
        --outfile $SPECLIB_NAME \
        --model onnx
fi

cargo run --release --bin timsseek -- \
    --config config_use.json \
    --speclib-file $SPECLIB_NAME \
    --output-dir $RESULTS_DIR \
    --dotd-file $DOTD_FILE $EXTRAS

# Technically does T/D competition
uv run -s showscores.py --results_dir $RESULTS_DIR --output_dir $SUMMARY_DIR