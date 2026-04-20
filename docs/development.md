# Development

Reference for working on timsbuktoolkit. All binaries ship `--help` for CLI flags. Environment variables are NOT listed by `--help`; see the Env vars table below.

## Binaries

| Binary | Crate | Purpose |
|--------|-------|---------|
| `timsseek` | `timsseek_cli` | Peptide-centric search |
| `timsseek_sample_speclib` | `timsseek_cli` | Speclib sampling utility |
| `timsquery_cli` | `timsquery_cli` | Low-level timsTOF query |
| `speclib_build` | `speclib_build_cli` | Build speclib from FASTA via Koina |
| `timsquery_viewer` | `timsquery_viewer` | GUI viewer for query results |

Run any with `--help` for the full flag list.

## Cargo features

| Feature | Crate | Effect | Use case | Enable |
|---------|-------|--------|----------|--------|
| `parallel` / `rayon` | `timsseek_cli` / `timsseek` | Rayon parallel scoring | Default; fastest wall-time | On by default |
| `instrumentation` | `timsseek_cli` / `timsseek` | `tracing-profile` perfetto spans | Perf tracing. **Requires `--no-default-features`** — the perfetto backend captures only the main thread, so rayon worker spans are dropped entirely. Run serial or traces for the hot path are empty. | `--features instrumentation --no-default-features` |
| `track-alloc` | `timsseek_cli` | Global allocator tracking via `alloc_track` | Binary prints per-phase allocation deltas to stderr: `[alloc] <phase> d_bytes=... d_live=... churn=... peak=... hist=...`. Detect churn + memory regressions. Dev-only; do not ship. | `--features track-alloc` |
| `query-instr` | `timscentroid` | Per-peak atomic counters in `IndexedPeakGroup::for_each_peak` | Filter-funnel shape + pass rates. ~10× wall-time inflation — funnel counts only, not timing. | `-p timscentroid --features query-instr` |
| `aws` / `gcp` / `azure` | `timscentroid` | `object_store` cloud backends | Read `.d` / speclib from cloud | `--features aws` (etc.) |

## Env vars

Not shown by `--help`. Read directly via `std::env::var`.

| Env var | Binaries | Default | Purpose |
|---------|----------|---------|---------|
| `RUST_LOG` | all | `info` (where defaulted) | `tracing-subscriber` EnvFilter. Examples: `RUST_LOG=debug`, `RUST_LOG=timsseek=trace,timscentroid=debug`. |
| `BUCKET_SIZE` | `timsseek` | `256` | Overrides peak-index rebucket size after load. Raw `.d` files ship with `bucket_size=4096`, too large for tight mz tolerances. Lower → faster Phase 1/3 (~−24% wall at 256). For perf experiments. |
| `TIMSCENTROID_WORKER_THREADS` | any using `timscentroid` (`timsseek`, `timsquery_cli`, `timsquery_viewer`) | `8` | Tokio runtime worker threads for cloud (`object_store` / S3) reads. Bump for higher remote concurrency. |

## Taskfile

`task --list-all` enumerates everything. Non-obvious ones:

- `task test`, `task fmt`, `task clippy` — `task fmt` runs nightly rustfmt + ruff. Do not use `cargo fmt` (stable silently drops nightly-only opts).
- `task speclib:build -- <args>` — wrapper around `speclib_build`.
- `task speclib:local-koina` / `task speclib:stop-koina` — local Koina docker. First run downloads all models (~10-30 min).
- `task docker` — cross-builds linux/amd64 images.
- `task license_check`, `task todos`, `task bumpver`, `task build_python`.

Per-crate: `rust/timsseek/Taskfile.yml` adds a watch loop (`task timsseek`) — rebuild + test + fmt + clippy on source change.

## Tracked scripts

| Path | Purpose |
|------|---------|
| `bench/wandb_bench.py` | wandb-logged benchmark runner |
| `scripts/release.sh` | Release cut helper |
| `Dockerfile` | Multi-stage image (used by `task docker`) |
