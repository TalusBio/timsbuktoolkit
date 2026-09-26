# Development

Reference for working on timsbuktoolkit. All binaries ship `--help` for CLI flags. Environment variables are NOT listed by `--help`; see the Env vars table below.

## Binaries

| Binary | Crate | Purpose |
|--------|-------|---------|
| `timsseek` | `timsseek_cli` | Peptide-centric search |
| `timsquery_cli` | `timsquery_cli` | Low-level timsTOF query |
| `timsquery_viewer` | `timsquery_viewer` | GUI viewer for query results |
| `calib_dash` | `calib_dash` | Replay a saved `calibration.json` in the RT-calibration dashboard |

Run any with `--help` for the full flag list.

## Library inputs

`timsquery_cli` uses the shared library reader registry: DIA-NN `.speclib`,
TSV/TXT and Parquet; Spectronaut TSV; Skyline CSV; mzSpecLib (including gzip);
and target JSON. Format detection also inspects contents; the
[reader registry](../rust/timsquery/src/serde/library_file.rs) defines dispatch.

`timsseek` and `timsquery_viewer` use the same registry through `ReferenceLibrary`.
Every loaded scoring library contains geometry usable for query extraction.
Scoring requires retained fragments and aligned reference intensities. mzSpecLib
can supply these without sequence or fragment annotations. Opaque peaks disable
fragment-isotope scores and automatic mass-shift decoy generation library-wide;
supplied decoys remain usable. Without decoys, search scores the full acquisition
RT range and writes raw scores, omitting competition, discriminant-score and
q-value columns (`result_mode=raw` Parquet metadata). It bypasses calibration,
rescoring and q-value filtering. Supplied decoys do not by themselves validate a
decoy strategy for a new analyte class.

Library RT is declared once per library: seconds, normalized index, unspecified
coordinates, or absent. Measured minutes convert to seconds; normalized indices
retain their values, including zero and negative values. Declared normalized-scale
metadata is preserved from mzSpecLib headers and prediction provenance; records
select the actual axis. Mixed RT availability or
axes are rejected. An RT-free library searches unrestricted on RT, omits RT
prediction/residual features, and retains decoy competition, rescoring and q-values
when decoys are available. RT, m/z and mobility calibrate independently. Without
an RT fit, score-selected apexes still supply m/z and mobility measurements,
without ridge filtering. RT-free prescoring retains bounded candidates across
observed RT bands. Mobility calibration requires a searchable run axis; FAIMS and
absent mobility do not contribute sentinel measurements. An axis without measurements keeps its configured tolerance
in primary and secondary extraction. A failed RT fit never reinterprets an index
as seconds; initial extraction searches unrestricted RT. Secondary queries still
center on detected apex RT and observed mobility.

Results format 5 includes `library_rt_axis` Parquet metadata. Unavailable library
RT, calibrated RT and RT residuals are NaN; observed apex RT remains seconds.
Calibration format v4 records the input library axis and effective tolerance
enums; an empty RT snapshot means no RT fit, while residuals can still contain
m/z and mobility calibration. Older calibration files
must be regenerated. Explicit `rt_seconds` JSON/Python inputs remain seconds. JSON `rt_axis` is
preserved whether optional precursor/fragment labels are supplied or filled in.
Rust source geometry exposes `LibraryRT<f32>` through `Target` and `OwnedTarget`.
`ExtractionQuery` borrows that geometry with a separately resolved acquisition
window. `RtSelection::Centered(ObservedRTSeconds(...))` or `FullRun` resolves
against the acquisition extent once. Collectors copy that window and reuse their
geometry/intensity buffers; they retain no source borrow. Index queries accept
`PeakTolerance` for m/z, mobility and quadrupole only. Isotope extraction shifts
collector buffers directly, without an owned target or library-axis clone.
The viewer searches full-run without an RT fit, including libraries measured in
seconds. Direct-query CLI restricted RT inputs explicitly mean acquisition
seconds; normalized library indices require calibration first. Python target
`rt_seconds` likewise specifies an acquisition coordinate.

`run_report.json` records the resolved `scoring_plan` once for the search library,
plus `calibration_scoring_plan` when a separate calibration library is supplied.
These are the same plans serialized in Parquet metadata: sequence-operation
coverage, precursor-isotope method/reasons, fragment-isotope availability, and
`decoys` (requested policy, resolved strategy, reason, stored-decoy count).
Per-file `pipeline.raw_scores` records whether raw scoring was used.


Programmatic `TargetTable::Str` accepts an optional intensity sidecar. Scoring
assigns unique packed unknown keys (currently at most 255 peaks per entry),
preserving the original opaque labels separately; even a string such as `y3`
is not interpreted as chemistry. The query reader's target-list JSON schemas
(`Target` and `ElutionGroupInput` arrays) still supply geometry without reference
intensities, so those reader routes remain extraction-only.

Standalone `calib_dash` reads saved `calibration.json`, not a spectral library.

[Analyte module documentation](../rust/timsquery/src/chemistry/analyte.rs) documents chemistry storage, reader mappings,
library-wide sequence eligibility, and results format version 5.

Precursor isotope envelopes retain the three-bin C/S approximation. Scoring
finalization includes known modification C/S deltas or an explicitly based
molecular formula. If any stored target or decoy lacks usable counts, every
entry uses mass-estimated C/S. Generated mass-shift decoys reuse the parent's
envelope. The selected method and unavailable-count reasons appear in the
shared CLI/viewer plan report and Parquet scoring-plan metadata. Other elements
are ignored by this approximation; isotope-labelled C/S and unspecified formula
bases cannot supply its composition counts.

## Cargo features

| Feature | Crate | Effect | Use case | Enable |
|---------|-------|--------|----------|--------|
| `parallel` / `rayon` | `timsseek_cli` / `timsseek` | Rayon parallel scoring | Default; fastest wall-time | On by default |
| `instrumentation` | `timsseek_cli` / `timsseek` | `tracing-profile` perfetto spans | Perf tracing. **Requires `--no-default-features`** -- the perfetto backend captures only the main thread, so rayon worker spans are dropped entirely. Run serial or traces for the hot path are empty. | `--features instrumentation --no-default-features` |
| `track-alloc` | `timsseek_cli` | Global allocator tracking via `alloc_track` | Binary prints per-phase allocation deltas to stderr: `[alloc] <phase> d_bytes=... d_live=... churn=... peak=... hist=...`. Detect churn + memory regressions. Dev-only; do not ship. | `--features track-alloc` |
| `dashboard` | `timsseek_cli` / `rescore_dash` | Ratatui TUI of a rescoring run: score separation, per-feature histograms, FDR and calibration curves | Interactive dev inspection. Dev-only; also needs `TIMSSEEK_RESCORE_DASHBOARD` at runtime (below). | `--features dashboard` |
| `calib-dashboard` | `timsseek_cli` | Pulls in `calib_dash`, wiring an interactive terminal dashboard into Phase 1/2 of RT calibration | Step through Phase 1 prescore batches, watch the calibration curve/DP path converge, inspect the Phase 2 fit and derived tolerances. Does nothing on its own -- also requires `TIMSSEEK_CALIB_DASHBOARD=1` at runtime (see Env vars). Dev-only; do not ship. | `--features calib-dashboard` |
| `query-instr` | `timscentroid` | Per-peak atomic counters in `IndexedPeakGroup::for_each_peak` | Filter-funnel shape + pass rates. ~10× wall-time inflation -- funnel counts only, not timing. | `-p timscentroid --features query-instr` |
| `aws` / `gcp` / `azure` | `timscentroid` | `object_store` cloud backends | Read `.d` / speclib from cloud | `--features aws` (etc.) |

## Env vars

Not shown by `--help`. Read directly via `std::env::var` / `var_os`.

| Env var | Binaries | Default | Purpose |
|---------|----------|---------|---------|
| `RUST_LOG` | all | `info` (where defaulted) | `tracing-subscriber` EnvFilter. Examples: `RUST_LOG=debug`, `RUST_LOG=timsseek=trace,timscentroid=debug`. |
| `BUCKET_SIZE` | `timsseek` | `256` | Overrides peak-index rebucket size after load. Raw `.d` files ship with `bucket_size=4096`, too large for tight mz tolerances. Lower → faster Phase 1/3 (~−24% wall at 256). For perf experiments. |
| `TIMSCENTROID_WORKER_THREADS` | any using `timscentroid` (`timsseek`, `timsquery_cli`, `timsquery_viewer`) | `8` | Tokio runtime worker threads for cloud (`object_store` / S3) reads. Bump for higher remote concurrency. |
| `TIMSSEEK_RESCORE_DASHBOARD` | `timsseek`, built with `--features dashboard` | unset (off) | Any value but `0`/`false` opens the rescore dashboard TUI after the Phase 6 write; without `--features dashboard` the check is compiled out entirely. Warns and skips when stdout is not a terminal. Invocation: `cargo run -r --bin timsseek --features dashboard -- ...`. |
| `TIMSSEEK_CALIB_DASHBOARD` | `timsseek_cli` | unset | Set to `1`, `true` or `yes` to open the interactive RT-calibration dashboard. Inert unless the binary was built `--features calib-dashboard`. |
| `CALIB_DASH_FRAME_BUDGET_MB` | `timsseek_cli` | `64` | Caps the dashboard's retained-frame slab, in megabytes. |

## Taskfile

`task --list-all` enumerates everything. Non-obvious ones:

- `task test`, `task fmt`, `task clippy` -- `task fmt` runs nightly rustfmt + ruff. Do not use `cargo fmt` (stable silently drops nightly-only opts).
- `task docker` -- cross-builds linux/amd64 images.
- `task license_check`, `task todos`, `task bumpver`, `task build_python`.

Per-crate: `rust/timsseek/Taskfile.yml` adds a watch loop (`task timsseek`) -- rebuild + test + fmt + clippy on source change.

## S3 staging

`timsseek_cli` config:

```toml
[staging]
tempdir_root = "/scratch/timsseek"   # default: system temp
max_prefix_keys = 256
save_sidecar = false                  # write .idx next to raw input
stale_sweep_age_hours = 24            # 0 disables startup sweep
```

Startup sweeps `timsseek-staging-*` subdirs older than threshold that lack a `.lock` sentinel. Reclaims tempdirs from SIGKILL'd/crashed runs.

Env vars (via `object_store` default chain): `AWS_{ACCESS_KEY_ID,SECRET_ACCESS_KEY,SESSION_TOKEN,REGION}`, `AWS_ENDPOINT_URL` (MinIO/R2), `AWS_S3_FORCE_PATH_STYLE` (auto-on with endpoint).

Enable with `--features aws` on `timscentroid` / `tims_stage`. Default build omits AWS SDK.

MinIO smoke test (needs pre-seeded bucket):

```bash
MINIO_TEST_ENDPOINT=http://localhost:9000 \
AWS_ACCESS_KEY_ID=minioadmin AWS_SECRET_ACCESS_KEY=minioadmin \
MINIO_TEST_BUCKET=tims-stage-ci \
cargo test -p tims_stage --features aws --test minio_smoke
```

## Tracked scripts

| Path | Purpose |
|------|---------|
| `bench/wandb_bench.py` | wandb-logged benchmark runner |
| `bench/entrapment.py` | Generate peptide entrapments, run timsseek, estimate FDP; see [entrapment.md](entrapment.md) |
| `scripts/release.sh` | Release cut helper |
| `Dockerfile` | Multi-stage image (used by `task docker`) |
