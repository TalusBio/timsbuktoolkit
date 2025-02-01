# Timsbuktoolkit

A high-performance toolkit for processing and analyzing timsTOF mass spectrometry data.

> ⚠️ **Development Status - Alpha**: This project is currently under active development. APIs and features may change frequently. While functional, it should be considered experimental software.

The main intent of this project is to provide a platform performant
and transparent way to query and analyze timsTOF mass spectrometry data.

## Overview

timsseek is a collection of Rust-based tools designed for efficient processing and analysis of timsTOF
mass spectrometry data. The project consists of several components:

- `timsquery`
    - Implements a series of modular indices+aggregators+queries that can be used to query timsTOF data.
    - It also compoiles to a cli that can be used to query the data.
- `timsseek`: Implement spectral library reading+build and core logic to score peptide-data matches.
- `timsseek_cli`: Command-line interface for a peptide-centric search engine.
- `timsseek_rts`
    - Command-line program that starts a server where on-demand search of peptides can be performed.
    - It also incluides an example receiver server in srteamlit (python) to show how to interface with it.

## Installation

### Prerequisites

- Rust (latest stable version)
- UV (for all python-related tasks)

### Building from Source

1. Clone the repository:
```bash
git clone https://github.com/TalusBio/timsbuktoolkit.git
cd timsbuktoolkit
```

2. Build the Rust components:
```bash
cargo build --release
```

## Usage

Each component has a different usage pattern.

### Command Line Interface

#### Timsseek

To run timsquery we need a spectral library and a configuration file and a raw
data file.

The current implementation of the speclib is an ndjson file
(we also have a builder for the library ... I am happy to
integrate other sources of predictions for it.)

```bash
DOTD_FILE="$HOME/data/my_data.d"
FASTA_FILE="$HOME/fasta/VIMENTIN.fasta"
SPECLIB_NAME="vimentin.ndjson"
RESULTS_DIR="vimentin_search_results"
SUMMARY_DIR="vimentin_search_summary"

# Write the config file
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

# Build the spectral lib
# Rn the models for RT+mobility are pretty rudimentary and
# hard-coded for a 22 min gradient, we can improve them in the future.
uv run speclib_build_fasta \
    --fasta_file $FASTA_FILE \
    --decoy_strategy REVERSE \
    --max_ions 10 \
    --outfile $SPECLIB_NAME \
    --model onnx

# Run timsseek using the generated speclib + config
cargo run --release --bin timsseek -- \
    --config config_use.json \
    --speclib-file $SPECLIB_NAME \
    --output-dir $RESULTS_DIR \
    --dotd-file $DOTD_FILE $EXTRAS

# Rn this is kind of an ugly script that runs some summary plotting
# and target-decoy competitions.
uv run -s showscores.py --results_dir $RESULTS_DIR --output_dir $SUMMARY_DIR
```

#### On-Demand Search

```bash
RAW_FILE=$HOME/data/mysupercoolfile.d

# Write the config file
cat << EOF > tolconfig.json
{
    "ms": {"ppm":  [15.0, 15.0]},
    "mobility": {"percent": [10.0, 10.0]},
    "quad": {"absolute": [0.1, 0.1]}
}
EOF

# This initializes the server from the file.
# Depending on the system/file it might take ~7-30 seconds.
# To index the data
cargo run --bin timsseek_rts --release -- \
    --config ./tolconfig.json \
    --dotd-file $RAW_FILE &
SERVER_PID=$!

# To start the receiver, this sample app allows typing a peptide
# and visualizing the scores
uv run --project timsseek_rts/python/ --verbose streamlit run timsseek_rts/python/receiver.py
kill $SERVER_PID
wait

```

## Development

### Setting up the Development Environment

TODO

### Common Tasks

Most common tasks are defined in the `Taskfile.yml` file and can be run using the `task` command:

```bash
# Run tests
task test

# Format code
task fmt

# Run linter
task clippy

# Check dependencies
task license_check

# Run benchmarks
task bench
```

## License

This project is licensed under the Apache License, Version 2.0.

## Authors

- Sebastian Paez


## Contributing

Contrubutions are welcome and not all of them have to be code!
Some of the forms of contributing to the current state of the project could be:

- Requesting documentation
    - Since we wrote the project, it is very hard to see it from an user perspective
      so having people reminding us to document something is incredibly helpful.
- Docs
    - We are still working on the docs, but we welcome any help to improve them.
      Even suggestions on how to host/serve them would be very welcome!
- Reporting bugs
    - Since we are still in early development, there it little expectation of
    correctness or completeness but there might be several use cases-edge cases
    that have not been considered yet, we appreciate you reporting them.
- Ideas
    - If you have any idea how to improve the project, please let us know!
      We are more than happy to discuss whether it fits the scope of the project
      and evaluate how viable it would be to implement it!
- Code
    - We welcome pull requests! We would really appreciate if an issue is open 
      to discuss potential changes before they are merged.


