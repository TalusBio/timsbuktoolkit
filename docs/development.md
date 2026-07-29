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
| `dashboard` | `timsseek_cli` / `rescore_dash` | Ratatui TUI of a rescoring run: score separation, per-feature histograms, FDR and calibration curves | Interactive dev inspection. Dev-only; also needs `TIMSSEEK_RESCORE_DASHBOARD` at runtime (below). | `--features dashboard` |
| `query-instr` | `timscentroid` | Per-peak atomic counters in `IndexedPeakGroup::for_each_peak` | Filter-funnel shape + pass rates. ~10× wall-time inflation — funnel counts only, not timing. | `-p timscentroid --features query-instr` |
| `aws` / `gcp` / `azure` | `timscentroid` | `object_store` cloud backends | Read `.d` / speclib from cloud | `--features aws` (etc.) |

## Env vars

Not shown by `--help`. Read directly via `std::env::var`.

| Env var | Binaries | Default | Purpose |
|---------|----------|---------|---------|
| `RUST_LOG` | all | `info` (where defaulted) | `tracing-subscriber` EnvFilter. Examples: `RUST_LOG=debug`, `RUST_LOG=timsseek=trace,timscentroid=debug`. |
| `BUCKET_SIZE` | `timsseek` | `256` | Overrides peak-index rebucket size after load. Raw `.d` files ship with `bucket_size=4096`, too large for tight mz tolerances. Lower → faster Phase 1/3 (~−24% wall at 256). For perf experiments. |
| `TIMSCENTROID_WORKER_THREADS` | any using `timscentroid` (`timsseek`, `timsquery_cli`, `timsquery_viewer`) | `8` | Tokio runtime worker threads for cloud (`object_store` / S3) reads. Bump for higher remote concurrency. |
| `TIMSSEEK_RESCORE_DASHBOARD` | `timsseek`, built with `--features dashboard` | unset (off) | Any value but `0`/`false` opens the rescore dashboard TUI after the Phase 6 write; without `--features dashboard` the check is compiled out entirely. Warns and skips when stdout is not a terminal. Invocation: `cargo run -r --bin timsseek --features dashboard -- ...`. |
| `TIMSSEEK_LDA_DUMP` | `timsseek` | unset (off) | `=/some/prefix` writes the raw LINEAR-lane matrix for offline feature engineering: `<prefix>.f64` (`nrows`/`ncols` header then row-major LE `f64`), `<prefix>.labels` (`u8`, 1 = target), `<prefix>.names.txt`. Best-effort — I/O errors log and are skipped. |
| `TIMSSEEK_MLP_*` | `timsseek` | unset (compiled defaults) | Dev-only hyperparameter sweep overrides for `MlpConfig`, which is deliberately NOT TOML-exposed. See the table below. |

### `TIMSSEEK_MLP_*` (dev-only MLP hyperparameter overrides)

**Not a supported interface.** An escape hatch for sweeps, so a hyperparameter
experiment is not a recompile. Read in exactly one place,
`MlpConfig::from_env` (`rust/timsseek/src/ml/mlp.rs`), which `MlpConfig::default()`
calls in non-test builds. Only the three MLP rescore models (`mlp`, `mlp_all`,
`hybrid_mlp`) read them; the config file stays the authority for everything else.

| Env var | Field | Format |
|---------|-------|--------|
| `TIMSSEEK_MLP_EPOCHS` | `epochs` (upper bound, not a target) | positive integer — `60` |
| `TIMSSEEK_MLP_LR` | `lr` | finite positive float — `1e-3` |
| `TIMSSEEK_MLP_BATCH_SIZE` | `batch_size` | positive integer — `512` |
| `TIMSSEEK_MLP_HIDDEN` | `hidden` | comma-separated positive integers — `128,64`; or `none` for no hidden layer (bare linear model) |
| `TIMSSEEK_MLP_WEIGHT_DECAY` | `weight_decay` | finite non-negative float — `0`, `1e-3` |
| `TIMSSEEK_MLP_PATIENCE` | `early_stopping_patience` | positive integer — `8`; or `none` / `off` to disable early stopping entirely |
| `TIMSSEEK_MLP_SEED` | `seed` | `u64`, decimal or `0x` hex — `0x2545F4914F6CDD1D` |

- **Unset changes nothing.** With none set the config is bit-identical to the
  compiled default and nothing is logged.
- **`TIMSSEEK_MLP_PATIENCE` also changes the TRAINING-SET SIZE on the hybrid
  path.** `hybrid_mlp` cross-fits its `mlp_score` column through the `crossfit`
  driver, which has no fold to spare for validation, so a set patience carves 20%
  of each fold's train rows off for the stopping decision and the optimizer sees
  80%; `off` trains on 100%. A sweep arm comparing `PATIENCE=off` to
  `PATIENCE=8` on `hybrid_mlp` therefore varies two things at once. `mlp` /
  `mlp_all` are unaffected — they early-stop on the next fold, so every train row
  still trains.
- **A malformed value ABORTS the run** with the variable, the value and the
  expected format — it never falls back to the default, because a warned-past
  variable produces a sweep row labelled with a value that never trained. The
  CLI checks at startup (when an MLP model is selected) so the abort lands in
  the first second, not at Phase 5.
- **Any override active** logs the full effective config once per run at `info`,
  prefixed `MLP config: DEV OVERRIDE ACTIVE`.
- `loss` is deliberately absent: `MlpLoss::Focal` is two coupled floats that also
  interact with the per-row target/decoy sample weights, and wants a real config
  surface.

Per-fold results are one `info` line each, so a sweep is greppable:

```bash
grep -o 'MLP fold summary:.*' run.log
# MLP fold summary: fold=0 epochs_run=15 epoch_budget=60 kept_epoch=7 \
#   best_held_out_loss=0.467983 final_train_loss=0.437800 restored=true \
#   train_rows=38046 fit_rows=38046 inputs=92 lane_columns=101
```

`kept_epoch` is the 1-based epoch whose weights were restored; `kept_epoch=none` /
`best_held_out_loss=none` mean there was no held-out set. `inputs` < `lane_columns`
means columns were culled (also a `warn`).

`RUST_LOG=timsseek::ml::mlp=debug` additionally traces train/held-out loss every
epoch (`epochs` lines per fold) — for looking at a curve, not for building a table.

## Taskfile

`task --list-all` enumerates everything. Non-obvious ones:

- `task test`, `task fmt`, `task clippy` — `task fmt` runs nightly rustfmt + ruff. Do not use `cargo fmt` (stable silently drops nightly-only opts).
- `task speclib:build -- <args>` — wrapper around `speclib_build`.
- `task speclib:local-koina` / `task speclib:stop-koina` — local Koina docker. First run downloads all models (~10-30 min).
- `task docker` — cross-builds linux/amd64 images.
- `task license_check`, `task todos`, `task bumpver`, `task build_python`.

Per-crate: `rust/timsseek/Taskfile.yml` adds a watch loop (`task timsseek`) — rebuild + test + fmt + clippy on source change.

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
