mod cli;
mod config;
mod errors;
mod processing;

use clap::Parser;
use timsquery::utils::TupleRange;
use timsquery::{
    IndexedTimstofPeaks,
    IndexingCentroidingConfig,
    load_index,
};
use timsseek::scoring::Scorer;
use timsseek::scoring::timings::TimedStep;
use tracing::{
    error,
    info,
    info_span,
};
use tracing_subscriber::filter::EnvFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use tracing_subscriber::fmt::{
    self,
};
use tracing_subscriber::prelude::*;
use tracing_subscriber::{
    self,
};

#[cfg(feature = "instrumentation")]
use tracing_profile::{
    PrintTreeConfig,
    PrintTreeLayer,
};

use cli::Cli;
use config::{
    Config,
    InputConfig,
    OutputConfig,
};
// use tracing_profile::PerfettoLayer;

#[cfg(all(target_os = "windows", not(feature = "track-alloc")))]
use mimalloc::MiMalloc;

#[cfg(all(target_os = "windows", not(feature = "track-alloc")))]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

#[cfg(feature = "track-alloc")]
#[global_allocator]
static GLOBAL: alloc_track::TrackingAllocator = alloc_track::TrackingAllocator::new();

struct ValidatedInputs {
    raw_inputs: Vec<String>,
    speclib_uri: String,
    calib_lib_uri: Option<String>,
    output_uri: String,
    overwrite: bool,
}

use tims_stage::{
    expand_local_uri,
    is_remote_uri,
};

fn validate_inputs(
    config: &Config,
    args: &Cli,
) -> std::result::Result<ValidatedInputs, errors::CliError> {
    info!("Validating inputs and outputs before processing...");

    // `~` is not expanded by `Path::exists` / `File::open`; canonicalise once
    // so validation, staging and probes all see the same expanded form.
    let raw_inputs: Vec<String> = match config.analysis.raw_inputs.clone() {
        Some(files) => files.iter().map(|f| expand_local_uri(f)).collect(),
        None => {
            return Err(errors::CliError::Config {
                source: "No raw files provided, please provide raw_inputs in either the config file or with the --raw-inputs flag".to_string(),
            });
        }
    };

    let speclib_uri: String = match &config.input {
        Some(InputConfig::Speclib { uri }) => expand_local_uri(uri),
        None => {
            return Err(errors::CliError::Config {
                source: "No input specified".to_string(),
            });
        }
    };

    let output_uri: String = match &config.output {
        Some(output_config) => expand_local_uri(&output_config.uri),
        None => {
            return Err(errors::CliError::Config {
                source: "No output directory specified".to_string(),
            });
        }
    };

    // Local only: remote resolution happens at open.
    if !is_remote_uri(&speclib_uri) && !std::path::Path::new(&speclib_uri).exists() {
        return Err(errors::CliError::Io {
            source: "Speclib file does not exist".to_string(),
            path: Some(speclib_uri.clone()),
        });
    }
    info!("✓ Speclib URI: {}", speclib_uri);

    let calib_lib_uri: Option<String> = args.calib_lib.as_deref().map(expand_local_uri);
    if let Some(ref uri) = calib_lib_uri {
        if !is_remote_uri(uri) && !std::path::Path::new(uri).exists() {
            return Err(errors::CliError::Io {
                source: "Calibration library file does not exist".to_string(),
                path: Some(uri.clone()),
            });
        }
        info!("✓ Calibration library URI: {}", uri);
    }

    // Local-path existence only; remote URIs are resolved by tims_stage at
    // staging time.
    for raw_uri in &raw_inputs {
        if !is_remote_uri(raw_uri) && !std::path::Path::new(raw_uri.as_str()).exists() {
            return Err(errors::CliError::Io {
                source: "Raw file does not exist".to_string(),
                path: Some(raw_uri.clone()),
            });
        }
    }
    info!("✓ All {} raw input(s) validated", raw_inputs.len());

    // Remote outputs skip the writability probe — the upload is the write test.
    if !is_remote_uri(&output_uri) {
        let output_dir_path = std::path::Path::new(&output_uri);

        match std::fs::create_dir_all(output_dir_path) {
            Ok(_) => {
                info!("✓ Output directory ready");
            }
            Err(e) => {
                return Err(errors::CliError::Io {
                    source: format!("Failed to create output directory: {}", e),
                    path: Some(output_uri.clone()),
                });
            }
        };

        let test_file = output_dir_path.join(".write_test");
        match std::fs::write(&test_file, "test") {
            Ok(_) => {
                let _ = std::fs::remove_file(&test_file);
                info!("✓ Output directory is writable");
            }
            Err(e) => {
                return Err(errors::CliError::Io {
                    source: format!("Output directory is not writable: {}", e),
                    path: Some(output_uri.clone()),
                });
            }
        }
    }

    // Probe every artifact up-front so a collision fails before the heavy
    // analysis rather than after it.
    //
    // IMPORTANT — the two artifact lists below must stay in sync with the
    // writer sites (search for `ARTIFACT-LIST`):
    //   per-sample: processing.rs, main.rs overwrite-cleanup block
    //   run-level:  main.rs run report, OutputSink::finalize_run call site
    if !args.overwrite {
        let mut collisions: Vec<String> = Vec::new();
        for raw_uri in &raw_inputs {
            let sample = sample_name_from_uri(raw_uri).ok_or_else(|| errors::CliError::Io {
                source: "Unable to extract file stem".to_string(),
                path: Some(raw_uri.clone()),
            })?;
            // ARTIFACT-LIST (per-sample)
            for artifact in ["results.parquet", "performance_report.json"] {
                let uri = join_output_uri(&output_uri, &format!("{sample}/{artifact}"));
                if probe_uri_exists(&uri)? {
                    collisions.push(uri);
                }
            }
        }
        // ARTIFACT-LIST (run-level)
        for artifact in ["run_report.json", "config_used.json"] {
            let uri = join_output_uri(&output_uri, artifact);
            if probe_uri_exists(&uri)? {
                collisions.push(uri);
            }
        }
        if !collisions.is_empty() {
            let list = collisions
                .iter()
                .take(8)
                .map(|s| format!("  - {s}"))
                .collect::<Vec<_>>()
                .join("\n");
            let more = if collisions.len() > 8 {
                format!("\n  ... and {} more", collisions.len() - 8)
            } else {
                String::new()
            };
            return Err(errors::CliError::Config {
                source: format!(
                    "{} output artifact(s) already exist; pass --overwrite/-O to replace them:\n{list}{more}",
                    collisions.len()
                ),
            });
        }
        info!("✓ No output collisions (checked local + remote artifacts)");
    } else {
        info!("✓ --overwrite set; existing output artifacts will be replaced");
    }

    info!("All validations passed! Starting processing...");

    Ok(ValidatedInputs {
        raw_inputs,
        speclib_uri,
        calib_lib_uri,
        output_uri,
        overwrite: args.overwrite,
    })
}

/// Join an artifact path onto a base output URI. Remote URIs get a plain
/// `base/rel` concat; local paths go through `PathBuf::join` so OS-specific
/// separators are respected.
fn join_output_uri(base: &str, rel: &str) -> String {
    if is_remote_uri(base) {
        format!("{}/{}", base.trim_end_matches('/'), rel)
    } else {
        std::path::Path::new(base)
            .join(rel)
            .to_string_lossy()
            .to_string()
    }
}

/// Cheap existence probe: local `Path::exists` or remote HEAD.
fn probe_uri_exists(uri: &str) -> Result<bool, errors::CliError> {
    tims_stage::uri_exists(uri).map_err(|e| errors::CliError::Io {
        source: format!("existence probe {uri}: {e}"),
        path: Some(uri.to_string()),
    })
}

/// `sample.d.tar`, `sample.d/`, `sample.d.idx/` all collapse to `sample`.
fn sample_name_from_uri(uri: &str) -> Option<String> {
    let trimmed = uri.trim_end_matches('/');
    let mut stem = trimmed.rsplit('/').next()?;
    // Loop so chained suffixes collapse fully. Order matters: `.idx`/`.tar`
    // come off before `.d` so they can't leave a bare `.d` behind.
    loop {
        let before = stem;
        for ext in [".idx", ".tar", ".d"] {
            if let Some(s) = stem.strip_suffix(ext) {
                stem = s;
            }
        }
        if stem == before {
            break;
        }
    }
    if stem.is_empty() {
        None
    } else {
        Some(stem.to_string())
    }
}

fn get_frag_range_from_index(
    index: &IndexedTimstofPeaks,
) -> Result<TupleRange<f64>, errors::CliError> {
    // `fragmented_range` panics only when there are no ms2 window groups —
    // catch it and surface a readable "not a DIA run" error instead.
    let result =
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| index.fragmented_range()));
    match result {
        Ok(r) => Ok(r),
        Err(_) => Err(errors::CliError::DataReading {
            source: "Index has no MS2 window groups — is this a DIA run?".to_string(),
        }),
    }
}

fn process_single_file(
    raw_uri: &str,
    backend: &tims_stage::PerRunTempdir,
    save_sidecar: bool,
    speclib: &timsseek::data_sources::speclib::Speclib,
    calib_lib: Option<&timsseek::data_sources::speclib::Speclib>,
    config: &Config,
    base_output_dir: &std::path::Path,
    overwrite: bool,
    max_qvalue: f32,
    no_feature_stats: bool,
) -> std::result::Result<timsseek::scoring::PipelineReport, errors::CliError> {
    let file_name = std::path::Path::new(raw_uri)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(raw_uri);
    let _file_span = info_span!("file", name = file_name).entered();
    info!("Processing raw input: {}", raw_uri);

    let step = TimedStep::begin("Loading index");
    let (index, index_source) = load_index(
        raw_uri,
        backend,
        save_sidecar,
        IndexingCentroidingConfig::default(),
    )
    .map_err(|e| errors::CliError::Io {
        source: format!("load_index({raw_uri}): {e}"),
        path: Some(raw_uri.to_string()),
    })?;
    // Cache reuse vs raw build is otherwise invisible to the user.
    let load_index_ms = step.finish_with(format_args!("{index_source}")).as_millis() as u64;
    match index.mobility_kind() {
        timscentroid::MobilityKind::Ook0 => info!("Mobility axis: TIMS 1/K0 (searchable)"),
        timscentroid::MobilityKind::Absent => {
            info!("Mobility axis: none — mobility filter + score features disabled")
        }
        timscentroid::MobilityKind::Unsupported(label) => info!(
            "Mobility axis: unsupported [{label}] — mobility filter + score features disabled"
        ),
    }
    alloc_track::snap!("Loading index");

    // On-disk caches use `bucket_size=4096`, too large for the tight mz
    // tolerances of Phase 1 / Phase 3 (measured: ~−24% wall at 256).
    const REBUCKET_LEN: usize = 256;
    let step = TimedStep::begin(format_args!("Rebucket at {}", REBUCKET_LEN));
    let index = index.rebucket(REBUCKET_LEN);
    step.finish();

    let fragmented_range = get_frag_range_from_index(&index)?;

    let pipeline = Scorer {
        index,
        broad_tolerance: config.analysis.tolerance.clone(),
        fragmented_range,
    };

    let file_stem = sample_name_from_uri(raw_uri).ok_or_else(|| errors::CliError::Io {
        source: "Unable to derive sample name from URI".to_string(),
        path: Some(raw_uri.to_string()),
    })?;
    let file_output_dir = base_output_dir.join(&file_stem);

    std::fs::create_dir_all(&file_output_dir).map_err(|e| errors::CliError::Io {
        source: format!("Failed to create output subdirectory: {}", e),
        path: Some(file_output_dir.to_string_lossy().to_string()),
    })?;

    // ARTIFACT-LIST (per-sample): keep in sync with validate_inputs.
    if overwrite {
        let results_file = file_output_dir.join("results.parquet");
        if results_file.exists() {
            std::fs::remove_file(&results_file).map_err(|e| errors::CliError::Io {
                source: format!("Failed to remove existing results file: {}", e),
                path: Some(results_file.to_string_lossy().to_string()),
            })?;
        }

        let perf_report_file = file_output_dir.join("performance_report.json");
        if perf_report_file.exists() {
            std::fs::remove_file(&perf_report_file).map_err(|e| errors::CliError::Io {
                source: format!("Failed to remove existing performance report: {}", e),
                path: Some(perf_report_file.to_string_lossy().to_string()),
            })?;
        }
    }

    let file_output_config = OutputConfig {
        uri: file_output_dir.to_string_lossy().to_string(),
    };

    let report = processing::run_pipeline(
        speclib,
        calib_lib,
        &pipeline,
        config.analysis.chunk_size,
        &file_output_config,
        max_qvalue,
        load_index_ms,
        &config.calibration,
        no_feature_stats,
        config.analysis.rescore_model,
    )
    .map_err(|e| errors::CliError::DataReading {
        source: format!("{}", e),
    })?;

    info!("Successfully processed {}", raw_uri);
    Ok(report)
}

/// Local-vs-remote output routing: a remote `dest_uri` (s3://, gs://, az://)
/// writes into a tempdir and uploads per-sample; a local one writes directly.
struct OutputSink {
    dest_uri: String,
    working_dir: std::path::PathBuf,
    remote: bool,
    _tempdir: Option<tempfile::TempDir>,
}

impl OutputSink {
    fn new(dest_uri: &str) -> Result<Self, errors::CliError> {
        if is_remote_uri(dest_uri) {
            let td = tempfile::Builder::new()
                .prefix("timsseek-output-")
                .tempdir()
                .map_err(|e| errors::CliError::Io {
                    source: format!("output tempdir: {e}"),
                    path: None,
                })?;
            let working_dir = td.path().to_path_buf();
            Ok(Self {
                dest_uri: dest_uri.to_string(),
                working_dir,
                remote: true,
                _tempdir: Some(td),
            })
        } else {
            std::fs::create_dir_all(dest_uri).map_err(|e| errors::CliError::Io {
                source: format!("create output dir: {e}"),
                path: Some(dest_uri.to_string()),
            })?;
            Ok(Self {
                dest_uri: dest_uri.to_string(),
                working_dir: std::path::PathBuf::from(dest_uri),
                remote: false,
                _tempdir: None,
            })
        }
    }

    fn sample_dir(&self, sample: &str) -> std::path::PathBuf {
        self.working_dir.join(sample)
    }

    fn root(&self) -> &std::path::Path {
        &self.working_dir
    }

    /// Where a sample's files *end up*, not the working tempdir. Use for
    /// user-facing output so users don't see a tempdir that will be wiped.
    fn dest_uri_for_sample(&self, sample: &str) -> String {
        if self.remote {
            format!("{}/{}", self.dest_uri.trim_end_matches('/'), sample)
        } else {
            self.sample_dir(sample).to_string_lossy().into_owned()
        }
    }

    /// Upload and remove a per-sample subdir after the sample has finished
    /// writing; no-op for local destinations.
    fn finalize_sample(&self, sample: &str) -> Result<(), errors::CliError> {
        if !self.remote {
            return Ok(());
        }
        let local = self.sample_dir(sample);
        let sample_dest = self.dest_uri_for_sample(sample);
        for entry in std::fs::read_dir(&local).map_err(|e| errors::CliError::Io {
            source: format!("read sample dir: {e}"),
            path: Some(local.to_string_lossy().to_string()),
        })? {
            let entry = entry.map_err(|e| errors::CliError::Io {
                source: format!("read dir entry: {e}"),
                path: None,
            })?;
            let bn = entry.file_name().to_string_lossy().to_string();
            let dest = format!("{sample_dest}/{bn}");
            tims_stage::upload_file(&entry.path(), &dest).map_err(|e| errors::CliError::Io {
                source: format!("upload {dest}: {e}"),
                path: None,
            })?;
        }
        std::fs::remove_dir_all(&local).map_err(|e| errors::CliError::Io {
            source: format!("cleanup sample dir: {e}"),
            path: Some(local.to_string_lossy().to_string()),
        })?;
        Ok(())
    }

    /// Upload named top-level files (run_report.json, config_used.json)
    /// that exist in the working dir; no-op for local destinations.
    fn finalize_run(&self, files: &[&str]) -> Result<(), errors::CliError> {
        if !self.remote {
            return Ok(());
        }
        for bn in files {
            let local = self.working_dir.join(bn);
            if !local.exists() {
                continue;
            }
            let dest = format!("{}/{}", self.dest_uri.trim_end_matches('/'), bn);
            tims_stage::upload_file(&local, &dest).map_err(|e| errors::CliError::Io {
                source: format!("upload {dest}: {e}"),
                path: None,
            })?;
        }
        Ok(())
    }
}

/// Load a Speclib from a local path or remote URI.
///
/// Remote URIs are streamed (not buffered — speclibs can be multi-GB) to a
/// tempfile keeping the original basename, because `Speclib::from_file`
/// sniffs the format by extension.
fn speclib_from_uri(
    uri: &str,
    decoy_strategy: timsseek::DecoyPolicy,
) -> Result<
    (
        timsseek::data_sources::speclib::Speclib,
        Option<tempfile::TempDir>,
    ),
    errors::CliError,
> {
    if !is_remote_uri(uri) {
        let path = std::path::Path::new(uri);
        let lib = timsseek::data_sources::speclib::Speclib::from_file(path, decoy_strategy)
            .map_err(|e| errors::CliError::Config {
                source: format!("Failed to load speclib {uri}: {:?}", e),
            })?;
        return Ok((lib, None));
    }

    let trimmed = uri.trim_end_matches('/');
    let fname = trimmed
        .rsplit('/')
        .next()
        .filter(|s| !s.is_empty())
        .unwrap_or("speclib.bin");
    let td = tempfile::Builder::new()
        .prefix("timsseek-speclib-")
        .tempdir()
        .map_err(|e| errors::CliError::Io {
            source: format!("speclib tempdir: {e}"),
            path: None,
        })?;
    let local = td.path().join(fname);
    tims_stage::download_to_file(uri, &local).map_err(|e| errors::CliError::Io {
        source: format!("download speclib {uri}: {e}"),
        path: Some(uri.to_string()),
    })?;
    let lib = timsseek::data_sources::speclib::Speclib::from_file(&local, decoy_strategy).map_err(
        |e| errors::CliError::Config {
            source: format!("Failed to load speclib {uri}: {:?}", e),
        },
    )?;
    Ok((lib, Some(td)))
}

/// Holds the `instrumentation` flush guard, which flushes aggregated spans on
/// drop and so must outlive every traced span. Type-erased via `Box<dyn Any>`
/// to avoid naming the `tracing_profile` type outside the feature gate.
struct TracingHandle {
    #[cfg(feature = "instrumentation")]
    _tree_guard: Box<dyn std::any::Any>,
    #[cfg(not(feature = "instrumentation"))]
    _empty: (),
}

/// Install the tracing subscriber. Log file resolution order:
///   1. `--log-path -`            → stderr-only, no file
///   2. `--log-path <p>`          → that path verbatim
///   3. default, local output     → `<output_dir>/timsseek-<ts>.log`
///   4. default, no/remote output → `./timsseek-<ts>.log` in CWD
///
/// The `YYYYMMDDTHHMMSS` local-time suffix avoids clobbering previous runs
/// sharing the directory. Logs never reach the terminal unless `--log-path -`
/// opts in.
fn init_tracing(args: &Cli, config: &Config) -> TracingHandle {
    let log_file_path: Option<std::path::PathBuf> = match args.log_path {
        Some(ref p) if p.to_str() == Some("-") => None,
        Some(ref p) => Some(p.clone()),
        None => {
            let base: std::path::PathBuf = args
                .output_uri
                .as_ref()
                .or(config.output.as_ref().map(|o| &o.uri))
                .filter(|d| !is_remote_uri(d.as_str()))
                .map(std::path::PathBuf::from)
                .unwrap_or_else(|| {
                    std::env::current_dir().unwrap_or_else(|_| std::path::PathBuf::from("."))
                });
            let ts = chrono::Local::now().format("%Y%m%dT%H%M%S");
            Some(base.join(format!("timsseek-{ts}.log")))
        }
    };

    let build_env_filter = || {
        EnvFilter::builder()
            .with_default_directive(
                args.log_level
                    .parse()
                    .unwrap_or_else(|_| "info".parse().unwrap()),
            )
            .from_env_lossy()
            .add_directive("forust_ml=warn".parse().unwrap())
            .add_directive("timscentroid::storage=warn".parse().unwrap())
    };

    let (file_layer, stderr_warn_layer, stderr_all_layer) =
        if let Some(ref log_path) = log_file_path {
            if let Some(parent) = log_path.parent() {
                let _ = std::fs::create_dir_all(parent);
            }
            let log_file = std::fs::File::create(log_path).expect("Failed to create log file");
            let fl = fmt::layer()
                .with_writer(std::sync::Mutex::new(log_file))
                .with_filter(build_env_filter());
            let sl = fmt::layer()
                .with_writer(std::io::stderr)
                .without_time()
                .with_filter(tracing_subscriber::filter::LevelFilter::WARN);
            (Some(fl), Some(sl), None)
        } else {
            let sl = fmt::layer()
                .with_writer(std::io::stderr)
                .with_filter(build_env_filter());
            (None, None, Some(sl))
        };

    let spans_layer = log_file_path.as_ref().map(|log_path| {
        let mut spans_path = log_path.clone();
        let fname = spans_path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("timsseek.log")
            .to_string();
        spans_path.set_file_name(format!("{fname}.spans.jsonl"));
        let spans_file = std::fs::File::create(&spans_path).expect("Failed to create spans file");
        fmt::layer()
            .json()
            .with_span_events(FmtSpan::NEW | FmtSpan::CLOSE)
            .with_writer(std::sync::Mutex::new(spans_file))
            .with_filter(build_env_filter())
    });

    #[cfg(feature = "instrumentation")]
    let perf_filter = EnvFilter::builder()
        .with_default_directive("trace".parse().unwrap())
        .with_env_var("RUST_PERF_LOG")
        .from_env_lossy()
        .add_directive("forust_ml::gradientbooster=warn".parse().unwrap());

    #[cfg(feature = "instrumentation")]
    let events_filter = tracing_subscriber::filter::filter_fn(|metadata| !metadata.is_event());

    #[cfg(feature = "instrumentation")]
    let (tree_layer, _tree_guard) = PrintTreeLayer::new(PrintTreeConfig {
        attention_above_percent: 25.0,
        relevant_above_percent: 2.5,
        hide_below_percent: 0.0,
        display_unaccounted: true,
        no_color: false,
        accumulate_spans_count: false,
        accumulate_events: false,
        aggregate_similar_siblings: true,
    });
    #[cfg(feature = "instrumentation")]
    let tree_layer = tree_layer
        .with_filter(perf_filter)
        .with_filter(events_filter);

    let reg = tracing_subscriber::registry()
        .with(file_layer)
        .with(spans_layer)
        .with(stderr_warn_layer)
        .with(stderr_all_layer);

    #[cfg(feature = "instrumentation")]
    let reg = reg.with(tree_layer);

    reg.init();

    println!("timsseek v{}", env!("CARGO_PKG_VERSION"));
    match log_file_path {
        Some(ref log_path) => println!("Log: {}", log_path.display()),
        None => println!("Log: <stderr> (--log-path -)"),
    }
    println!();

    TracingHandle {
        #[cfg(feature = "instrumentation")]
        _tree_guard: Box::new(_tree_guard),
        #[cfg(not(feature = "instrumentation"))]
        _empty: (),
    }
}

fn main() {
    if let Err(e) = run() {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run() -> std::result::Result<(), errors::CliError> {
    let args = Cli::parse();

    if args.print_default_config {
        print!("{}", config::DEFAULT_CONFIG_TOML);
        return Ok(());
    }
    if let Some(ref path) = args.write_default_config {
        std::fs::write(path, config::DEFAULT_CONFIG_TOML).map_err(|e| errors::CliError::Io {
            source: e.to_string(),
            path: Some(path.to_string_lossy().to_string()),
        })?;
        eprintln!("Wrote default config to {}", path.display());
        return Ok(());
    }

    let mut config = match args.config {
        Some(ref config_path) => {
            let text = std::fs::read_to_string(config_path).map_err(|e| errors::CliError::Io {
                source: e.to_string(),
                path: Some(config_path.to_string_lossy().to_string()),
            })?;
            let is_toml = config_path
                .extension()
                .and_then(|e| e.to_str())
                .map(|e| e.eq_ignore_ascii_case("toml"))
                .unwrap_or(false);
            let parsed: Result<Config, String> = if is_toml {
                toml::from_str(&text).map_err(|e| e.to_string())
            } else {
                serde_json::from_str(&text).map_err(|e| e.to_string())
            };
            parsed.map_err(|e| errors::CliError::ParseError {
                msg: format!(
                    "Failed to parse config file {}: {e}\n\n\
                     Run `timsseek --print-default-config` for a reference template, \
                     or `--write-default-config <path>` to drop one to disk.\n\n\
                     Reference default:\n```toml\n{}```\n",
                    config_path.display(),
                    config::DEFAULT_CONFIG_TOML,
                ),
            })?
        }
        None => {
            info!("No config file provided, using default configuration");
            Config::default_config()
        }
    };

    // CLI args win over the config file.
    if !args.raw_inputs.is_empty() {
        config.analysis.raw_inputs = Some(args.raw_inputs.clone());
    }
    if let Some(ref speclib_uri) = args.speclib_uri {
        config.input = Some(InputConfig::Speclib {
            uri: speclib_uri.clone(),
        });
    }
    if config.input.is_none() {
        return Err(errors::CliError::Config {
            source: "No input provided, please provide one in either the config file or with the --speclib-uri flag".to_string(),
        });
    }
    if let Some(ref output_uri) = args.output_uri {
        config.output = Some(OutputConfig {
            uri: output_uri.clone(),
        });
    }

    if let Some(strategy) = args.decoy_strategy {
        config.analysis.decoy_strategy = strategy;
    }
    if let Some(model) = args.rescore_model {
        config.analysis.rescore_model = model.into();
    }

    // Held in `run()`'s scope so the instrumentation flush guard drops after
    // all work completes.
    let _tracing = init_tracing(&args, &config);

    info!("Parsed configuration: {:#?}", config.clone());
    alloc_track::snap!("start");

    let validated = validate_inputs(&config, &args)?;

    // The stale-tempdir sweep runs inside `PerRunTempdir::new`.
    let staging_cfg = config.staging.clone().unwrap_or_default();
    let save_sidecar_flag = staging_cfg.save_sidecar;
    let backend = tims_stage::PerRunTempdir::new(tims_stage::StagingConfig {
        tempdir_root: staging_cfg.tempdir_root.clone(),
        max_prefix_keys: staging_cfg.max_prefix_keys,
        stale_sweep_age_hours: staging_cfg.stale_sweep_age_hours,
    })
    .map_err(|e| errors::CliError::Io {
        source: format!("staging backend: {e}"),
        path: None,
    })?;

    let sink = OutputSink::new(&validated.output_uri)?;

    // ARTIFACT-LIST (run-level): keep in sync with validate_inputs.
    let config_output_path = sink.root().join("config_used.json");

    // Local-only; a remote destination just overwrites on upload.
    if validated.overwrite && config_output_path.exists() {
        std::fs::remove_file(&config_output_path).map_err(|e| errors::CliError::Io {
            source: format!("Failed to remove existing config file: {}", e),
            path: Some(config_output_path.to_string_lossy().to_string()),
        })?;
    }

    let config_json =
        serde_json::to_string_pretty(&config).map_err(|e| errors::CliError::ParseError {
            msg: format!("Failed to serialize config: {}", e),
        })?;
    std::fs::write(&config_output_path, config_json).map_err(|e| errors::CliError::Io {
        source: e.to_string(),
        path: Some(config_output_path.to_string_lossy().to_string()),
    })?;
    info!("Wrote final configuration to {:?}", config_output_path);

    let mut run_report = timsseek::scoring::RunReport::default();
    let mut failed_files: Vec<(String, errors::CliError)> = Vec::new();
    let mut successful_files: Vec<String> = Vec::new();

    // Loaded once, shared across all files.
    let step = TimedStep::begin("Loading speclib");
    info!(
        "Building database from speclib URI {}",
        validated.speclib_uri
    );
    info!(
        "Decoy generation strategy: {}",
        config.analysis.decoy_strategy
    );
    let (speclib, _speclib_td) =
        speclib_from_uri(&validated.speclib_uri, config.analysis.decoy_strategy)?;
    let load_speclib_ms = step
        .finish_with(format_args!(
            "{} entries, {:.1} frags/entry",
            speclib.len(),
            speclib.avg_fragments()
        ))
        .as_millis() as u64;
    alloc_track::snap!("Loading speclib");

    let (calib_lib, _calib_td, load_calib_lib_ms) = match &validated.calib_lib_uri {
        Some(uri) => {
            let step = TimedStep::begin("Loading calib lib");
            info!("Loading calibration library from {}", uri);
            let (lib, td) = speclib_from_uri(uri, config.analysis.decoy_strategy)?;
            let ms = step
                .finish_with(format_args!(
                    "{} entries, {:.1} frags/entry",
                    lib.len(),
                    lib.avg_fragments()
                ))
                .as_millis() as u64;
            alloc_track::snap!("Loading calib lib");
            (Some(lib), td, ms)
        }
        None => (None, None, 0),
    };

    run_report.load_speclib_ms = load_speclib_ms;
    run_report.load_calib_lib_ms = load_calib_lib_ms;
    run_report.speclib_entries = speclib.len();
    run_report.calib_lib_entries = calib_lib.as_ref().map_or(0, |l| l.len());

    let total_files = validated.raw_inputs.len();
    info!("Processing {} raw input(s)", total_files);

    for (idx, raw_uri) in validated.raw_inputs.iter().enumerate() {
        info!(
            "Processing input {} of {}: {}",
            idx + 1,
            total_files,
            raw_uri
        );

        let sample_name = match sample_name_from_uri(raw_uri) {
            Some(s) => s,
            None => {
                let e = errors::CliError::Io {
                    source: "Unable to derive sample name from URI".to_string(),
                    path: Some(raw_uri.clone()),
                };
                error!("Failed to process {}: {}", raw_uri, e);
                failed_files.push((raw_uri.clone(), e));
                continue;
            }
        };

        // Header/footer pair; `processing::run_pipeline` phase output lands
        // between them, so batched runs still show per-input wall time.
        println!("=== [{}/{}] {} ===", idx + 1, total_files, sample_name);
        let file_start = std::time::Instant::now();
        let sample_dest = sink.dest_uri_for_sample(&sample_name);

        match process_single_file(
            raw_uri,
            &backend,
            save_sidecar_flag,
            &speclib,
            calib_lib.as_ref(),
            &config,
            sink.root(),
            validated.overwrite,
            args.max_qvalue,
            args.no_feature_stats,
        ) {
            Ok(report) => {
                if let Err(e) = sink.finalize_sample(&sample_name) {
                    error!("Failed to finalize sample {}: {}", sample_name, e);
                    println!(
                        "=== [{}/{}] {} failed upload after {:?} ===",
                        idx + 1,
                        total_files,
                        sample_name,
                        file_start.elapsed()
                    );
                    run_report.status = timsseek::scoring::timings::RunStatus::Aborted;
                    run_report.abort_reason =
                        Some(format!("upload failure on sample {sample_name}: {e}"));
                    failed_files.push((raw_uri.clone(), e));
                    error!("Aborting batch due to upload failure");
                    break;
                }
                println!("Output: {sample_dest}");
                println!(
                    "=== [{}/{}] {} done in {:?} ===",
                    idx + 1,
                    total_files,
                    sample_name,
                    file_start.elapsed()
                );
                successful_files.push(raw_uri.clone());
                let mut outputs = vec![format!("{sample_dest}/results.parquet")];
                if !args.no_feature_stats {
                    outputs.push(format!(
                        "{sample_dest}/{}",
                        processing::FEATURE_STATS_FILENAME
                    ));
                    outputs.push(format!(
                        "{sample_dest}/{}",
                        processing::FEATURE_IMPORTANCE_FILENAME
                    ));
                }
                run_report.files.push(timsseek::scoring::FileReport {
                    file_name: raw_uri.clone(),
                    pipeline: report,
                    outputs,
                });
            }
            Err(e) => {
                error!("Failed to process {}: {}", raw_uri, e);
                println!(
                    "=== [{}/{}] {} FAILED after {:?}: {} ===",
                    idx + 1,
                    total_files,
                    sample_name,
                    file_start.elapsed(),
                    e
                );
                // I/O errors are likely systemic (disk full, permissions) —
                // abort the batch instead of failing every remaining file.
                if matches!(e, errors::CliError::Io { .. }) {
                    run_report.status = timsseek::scoring::timings::RunStatus::Aborted;
                    run_report.abort_reason = Some(format!("I/O error on {raw_uri}: {e}"));
                    failed_files.push((raw_uri.clone(), e));
                    error!("Aborting batch due to I/O error");
                    break;
                }
                failed_files.push((raw_uri.clone(), e));
            }
        }
    }

    // ARTIFACT-LIST (run-level): keep in sync with validate_inputs.
    let finalize_step = TimedStep::begin("Finalize run");
    // Must happen before serialization so the report self-describes where to
    // fetch everything.
    let dest_root = sink
        .dest_uri_for_sample("")
        .trim_end_matches('/')
        .to_string();
    run_report.artifacts = vec![
        format!("{dest_root}/run_report.json"),
        format!("{dest_root}/config_used.json"),
    ];
    let run_report_path = sink.root().join("run_report.json");
    if let Ok(json) = serde_json::to_string_pretty(&run_report) {
        let _ = std::fs::write(&run_report_path, json);
        info!("Wrote run report to {:?}", run_report_path);
    }

    // ARTIFACT-LIST (run-level): keep in sync with validate_inputs.
    sink.finalize_run(&["run_report.json", "config_used.json"])?;
    finalize_step.finish();

    info!("Successfully processed {} file(s)", successful_files.len());
    if !failed_files.is_empty() {
        error!("Failed to process {} file(s):", failed_files.len());
        for (file, err) in &failed_files {
            error!("  {}: {}", file, err);
        }
        return Err(errors::CliError::Config {
            source: format!("Failed to process {} file(s)", failed_files.len()),
        });
    }

    Ok(())
}

#[cfg(test)]
mod sample_name_tests {
    use super::sample_name_from_uri;
    #[test]
    fn local_dotd_plain() {
        assert_eq!(sample_name_from_uri("/data/run.d").as_deref(), Some("run"));
    }
    #[test]
    fn local_dotd_trailing_slash() {
        assert_eq!(sample_name_from_uri("/data/run.d/").as_deref(), Some("run"));
    }
    #[test]
    fn s3_tar_collapses_both_suffixes() {
        assert_eq!(
            sample_name_from_uri("s3://bkt/run.d.tar").as_deref(),
            Some("run")
        );
    }
    #[test]
    fn s3_idx_directory() {
        assert_eq!(
            sample_name_from_uri("s3://bkt/run.d.idx/").as_deref(),
            Some("run")
        );
    }
}
