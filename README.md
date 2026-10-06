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

Libraries are written as mzSpecLib. DIA-NN (`.speclib`, TSV, parquet),
Spectronaut, Skyline and JSON libraries are all read as well.

```bash
DOTD_FILE="$HOME/data/my_data.d"
FASTA_FILE="$HOME/fasta/VIMENTIN.fasta"
SPECLIB_NAME="vimentin.mzspeclib.txt.gz"
RESULTS_DIR="vimentin_search_results"

# Predict the library locally. No network and no server: the model is compiled
# into the binary. The `--out` suffix picks the format, and `.gz` compresses in
# the writer thread, so the uncompressed library never lands on disk.
cargo run --release --bin timsseek -- build-library \
    --fasta $FASTA_FILE \
    --fixed-mod "C[UNIMOD:4]" \
    --max-fragments 10 \
    -o $SPECLIB_NAME

# Or predict exact modified peptides and charge states from a TSV.
# protein_ids separates multiple IDs with semicolons.
cat > peptides.tsv <<'EOF'
proforma	protein_ids	decoy	decoy_group
PEC[UNIMOD:4]TIDEK/2	P1	false	cam
PECTIDEK/3	P1	false	plain
PEK[UNIMOD:259]TIDEK/2	PRTC	false	heavy
PECTDIEK/3	P1	true	plain
EOF
timsseek build-library --peptides peptides.tsv \
    --decoy-method shuffle --decoy-seed 42 \
    -o peptides.mzspeclib.txt.gz

# Run timsseek. Config is optional; defaults work for most runs.
# To tweak tolerances: `timsseek --write-default-config config.toml`, edit, pass with `-c`.
# TOML and JSON both accepted (sniffed by extension).
cargo run --release --bin timsseek -- \
    --speclib-uri $SPECLIB_NAME \
    --output-uri $RESULTS_DIR \
    --raw-inputs $DOTD_FILE
```

For a one-off search, skip the library file entirely. `--fasta` with no
`--speclib-uri` predicts the library into memory and searches it, so the only
thing written is the results:

```bash
cargo run --release --bin timsseek -- \
    --fasta $FASTA_FILE \
    --output-uri $RESULTS_DIR \
    --raw-inputs $DOTD_FILE
```

The `[library]` section of the configuration file supplies the prediction
settings, exactly as it does for `build-library`, and the model and FASTA digests
land in `config_used.json` so a run with no library file is still traceable.
Naming both a library and a FASTA loads the library and compares its recorded
build settings with the requested prediction settings. Differences are warned
about but do not prevent the search. Remote libraries stage their conventional
`.config.json` sidecar for this comparison when the format has no provenance
header.

To predict the library on the first run and read it on every run after, add
`--build-if-missing`. It writes the library to the path `--speclib-uri` names
and then searches it; a library that is already there is simply opened, so the
same command is safe to leave in a script:

```bash
cargo run --release --bin timsseek -- \
    --speclib-uri $SPECLIB_NAME \
    --fasta $FASTA_FILE \
    --output-uri $RESULTS_DIR \
    --raw-inputs $DOTD_FILE \
    --build-if-missing
```

With no `--speclib-uri`, the library is written to
`$RESULTS_DIR/<FASTA stem>.mzspeclib.txt.gz`, with the same `.config.json`
sidecar `build-library` writes, and an earlier run's library there is reported
as a collision unless `--overwrite` is passed. Without the flag, a library that
is not there is an error rather than a prediction: a mistyped path should not
cost minutes of prediction. Libraries are written through the filesystem, so a
remote `--speclib-uri` or output directory that would have to be built into is
rejected by name.

An overwrite compares the existing library's provenance with the settings about
to replace it and reports whether they match before publishing the new artifact.

## S3 inputs

A search accepts `s3://` URIs anywhere a path is accepted (AWS / MinIO / R2). `.d` can be a directory, `.tar`, or S3 prefix; `.idx` sidecars short-circuit staging.

```bash
timsseek --raw-inputs s3://bkt/sample.d.tar \
         --speclib-uri s3://bkt/lib.mzspeclib.txt.gz \
         --output-uri s3://bkt/runs/out
```

A proteome build takes minutes and reports as it goes: a progress bar per phase on a terminal, a log line every 30s when stderr is redirected. Per-phase wall times land in the `.config.json` sidecar.

Not supported by `build-library` right now, having been dropped along with the separate `speclib_build_cli`:

- **Remote paths.** `--fasta`, `--peptides` and `--out` are filesystem paths; a remote URI is rejected by name rather than staged. Build locally and copy.
- **Acquisition and chromatography context.** A build uses the model artifact's own defaults; picking a different one is a decision about the model.
- **Fragment and precursor filters** beyond `--min-intensity` and `--max-fragments`: no minimum transition count, and no precursor or fragment *m/z* bounds.
- **Decoy generation.** FASTA decoys are pseudo-reversed. Peptide TSV builds support `--decoy-method pseudo-reverse|shuffle` and `--decoy-seed`. Generated decoys are predicted by default; `--no-decoys` disables generation but retains supplied decoys. Pseudo-reverse retries with successively shorter interior spans; shuffle retries with its seed. When no collision-free decoy exists, the builder warns and keeps the target unpaired. Supplied target and decoy rows must share a `decoy_group` and charge; each group has one target and at most one supplied decoy. Peptide TSV builds require an mzSpecLib output suffix so the groups survive serialization.

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
