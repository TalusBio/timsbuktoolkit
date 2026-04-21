# Timsbuktoolkit

A high-performance toolkit for processing and analyzing timsTOF mass spectrometry data.

> ⚠️ **Development Status - Alpha**: This project is currently under active development. APIs and features may change frequently. While functional, it should be considered experimental software.

The main intent of this project is to provide a platform performant
and transparent way to query and analyze timsTOF mass spectrometry data.

## Overview

Timsbuktoolkit is a collection of Rust-based tools designed for efficient processing and analysis of timsTOF
mass spectrometry data. The project consists of several components:

- `timsquery`: Library implementing a series of modular aggregators+queries that can be used to query timsTOF data.
- `timsquery_cli`: Command-line interface for querying timsTOF data using the timsquery library.
- `timsseek`: Implement spectral library reading+build and core logic to score peptide-data matches.
- `timsseek_cli`: Command-line interface for a peptide-centric search engine.
- `speclib_build_cli`: CLI that builds a spectral library from FASTA using Koina.
- `timscentroid`: Internal library for indexing and centroiding timsTOF data.
- `calibrt`: Internal library for retention time calibration.
- `alloc_track`: Dev-only tracking global allocator. Opt-in via `--features track-alloc` on `timsseek_cli`; emits per-phase allocation deltas on stderr.

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

To run timsseek we need a spectral library and a configuration file and a raw
data file.

The current implementation of the speclib is an ndjson file
(we also have a builder for the library ... I am happy to
integrate other sources of predictions for it.)

```bash
DOTD_FILE="$HOME/data/my_data.d"
FASTA_FILE="$HOME/fasta/VIMENTIN.fasta"
SPECLIB_NAME="vimentin.ndjson"
RESULTS_DIR="vimentin_search_results"

# Build the spectral lib using Koina (Prosit) for fragment/RT prediction.
# Requires network access to https://koina.wilhelmlab.org or a local Koina server.
cargo run --release -p speclib_build_cli -- \
    --fasta $FASTA_FILE \
    --fixed-mod "C[U:4]" \
    --max-ions 10 \
    -o $SPECLIB_NAME

# Run timsseek. Config is optional; defaults work for most runs.
# To tweak tolerances: `timsseek --write-default-config config.toml`, edit, pass with `-c`.
# TOML and JSON both accepted (sniffed by extension).
cargo run --release --bin timsseek -- \
    --speclib-file $SPECLIB_NAME \
    --output-dir $RESULTS_DIR \
    --dotd-files $DOTD_FILE
```

## S3 inputs

Both CLIs accept `s3://` URIs anywhere a path is accepted (AWS / MinIO / R2). `.d` can be a directory, `.tar`, or S3 prefix; `.idx` sidecars short-circuit staging.

```bash
timsseek --raw-inputs s3://bkt/sample.d.tar \
         --speclib-uri s3://bkt/lib.msgpack.zst \
         --output-uri s3://bkt/runs/out

speclib_build_cli --fasta s3://bkt/proteome.fasta \
                  --output s3://bkt/lib.msgpack.zst
```

Auth via AWS default chain. MinIO/R2: set `AWS_ENDPOINT_URL`. See `docs/development.md` for `[staging]` config + env var list.

## Development

See [docs/development.md](docs/development.md) for dev utilities, compile flags, env vars, Taskfile targets, and scripts.

## License

This project is licensed under the Apache License, Version 2.0.

## Authors

- Sebastian Paez


## Contributing

Contributions are welcome and not all of them have to be code!
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


