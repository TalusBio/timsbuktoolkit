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
The current scoring bridge requires ion-annotated fragments and reference
fragment intensities; extraction also accepts opaque fragment labels. The
query reader's target-list JSON schemas (`Target` and `ElutionGroupInput` arrays)
currently load geometry without a reference-intensity sidecar. This describes
those reader paths, not a restriction of JSON as an encoding.

Standalone `calib_dash` reads saved `calibration.json`, not a spectral library.

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
| `scripts/release.sh` | Release cut helper |
| `Dockerfile` | Multi-stage image (used by `task docker`) |
