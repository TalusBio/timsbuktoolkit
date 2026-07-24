//! Benchmark index loading through each `tims_stage::resolve()` branch
//! without the full search pipeline.
//!
//! Reads one URI from the `LOAD_BENCH_URI` env var. Resolves, stages (if
//! needed), loads/builds the index, prints explicit timings for each
//! phase, and exits. No scoring, no output writes.
//!
//! Use this to measure:
//! - S3 tar staging time vs prefix staging time vs local-idx load vs
//!   remote-idx load.
//! - Per-phase breakdown (resolve probe, stage, build/load).
//!
//! Build + run:
//!   cargo build --release --features aws -p timsseek_cli --example load_bench
//!   LOAD_BENCH_URI='s3://bkt/sample.d.tar' \
//!     ./target/release/examples/load_bench
//!
//! Env:
//!   LOAD_BENCH_URI (required) — URI of the raw input to resolve+load.
//!   RUST_LOG (optional)       — default "info,tims_stage=info,timscentroid=info".

use std::time::Instant;

use tims_stage::{
    PerRunTempdir,
    Resolved,
    StagingBackend,
    StagingConfig,
};
use timscentroid::{
    IndexedTimstofPeaks,
    IndexingCentroidingConfig,
    StorageLocation,
};
use tracing::info;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let filter = std::env::var("RUST_LOG")
        .unwrap_or_else(|_| "info,tims_stage=info,timscentroid=info".to_string());
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(true)
        .with_writer(std::io::stderr)
        .init();

    let uri = std::env::var("LOAD_BENCH_URI").map_err(|_| "LOAD_BENCH_URI env var required")?;

    info!(%uri, "load_bench start");
    let t_total = Instant::now();

    // Phase 1: resolve (pure + one HEAD for sidecar on non-.idx inputs).
    let t = Instant::now();
    let resolved = tims_stage::resolve(&uri)?;
    let resolve_ms = t.elapsed().as_millis();
    info!(
        elapsed_ms = resolve_ms,
        kind = kind_of(&resolved),
        "resolve done"
    );

    // Phase 2 + 3: load directly (idx) or build via the unified raw core.
    let cfg = IndexingCentroidingConfig::default();
    match resolved {
        Resolved::Idx { loc } => {
            let t = Instant::now();
            let _idx = IndexedTimstofPeaks::load_from_storage(loc)?;
            info!(
                elapsed_ms = t.elapsed().as_millis(),
                "load_from_storage done"
            );
        }
        Resolved::Raw { uri } => {
            let backend = PerRunTempdir::new(StagingConfig::default())?;
            let t = Instant::now();
            let raw = tims_stage::load_raw(&uri, &backend, &cfg)?;
            info!(
                elapsed_ms = t.elapsed().as_millis(),
                reader = raw.reader_name,
                "load_raw done"
            );
        }
        Resolved::Tar { spec } => {
            let backend = PerRunTempdir::new(StagingConfig::default())?;
            let t = Instant::now();
            let staged = backend.stage(&spec)?;
            info!(elapsed_ms = t.elapsed().as_millis(), dotd = ?staged.as_ref(), "stage (tar) done");

            let t = Instant::now();
            let raw = tims_stage::load_raw(staged.as_ref().to_str().unwrap(), &backend, &cfg)?;
            info!(
                elapsed_ms = t.elapsed().as_millis(),
                reader = raw.reader_name,
                "load_raw (staged) done"
            );
        }
    };

    info!(total_ms = t_total.elapsed().as_millis(), "load_bench end");
    Ok(())
}

fn kind_of(r: &Resolved) -> &'static str {
    match r {
        Resolved::Idx { .. } => "Idx",
        Resolved::Raw { .. } => "Raw",
        Resolved::Tar { spec } => match spec {
            tims_stage::SourceSpec::S3Tar { .. } => "Tar(S3)",
            tims_stage::SourceSpec::LocalTar { .. } => "Tar(Local)",
        },
    }
}

// Silence unused import if StorageLocation is re-exported only when feature is on.
#[allow(dead_code)]
fn _unused(_: Option<StorageLocation>) {}
