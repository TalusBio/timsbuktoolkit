mod cli;
mod config;
mod errors;
mod processing;

use clap::Parser;
use timsquery::utils::TupleRange;
use timsquery::{
    CentroidingConfig,
    IndexedTimstofPeaks,
    load_index,
};
use timsseek::scoring::Scorer;
use timsseek::scoring::timings::TimedStep;
use tracing::{
    error,
    info,
};
use tracing_subscriber::filter::EnvFilter;
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

/// Validated inputs ready for processing
struct ValidatedInputs {
    raw_inputs: Vec<String>,
    speclib_uri: String,
    calib_lib_uri: Option<String>,
    output_uri: String,
    overwrite: bool,
}

fn is_remote_uri(uri: &str) -> bool {
    uri.starts_with("s3://") || uri.starts_with("gs://") || uri.starts_with("az://")
}

/// Validates all inputs and outputs before processing begins.
/// Returns ValidatedInputs on success, or an error if any validation fails.
fn validate_inputs(
    config: &Config,
    args: &Cli,
) -> std::result::Result<ValidatedInputs, errors::CliError> {
    info!("Validating inputs and outputs before processing...");

    // Get list of raw files to process
    let raw_inputs = match config.analysis.raw_inputs.clone() {
        Some(files) => files,
        None => {
            return Err(errors::CliError::Config {
                source: "No raw files provided, please provide raw_inputs in either the config file or with the --raw-inputs flag".to_string(),
            });
        }
    };

    // Get speclib URI
    let speclib_uri = match &config.input {
        Some(InputConfig::Speclib { uri }) => uri.clone(),
        None => {
            return Err(errors::CliError::Config {
                source: "No input specified".to_string(),
            });
        }
    };

    // Get output URI (local path or s3:// etc.)
    let output_uri: String = match &config.output {
        Some(output_config) => output_config.uri.clone(),
        None => {
            return Err(errors::CliError::Config {
                source: "No output directory specified".to_string(),
            });
        }
    };

    // Validate speclib exists (local only — remote resolution happens at open)
    if !is_remote_uri(&speclib_uri) && !std::path::Path::new(&speclib_uri).exists() {
        return Err(errors::CliError::Io {
            source: "Speclib file does not exist".to_string(),
            path: Some(speclib_uri.clone()),
        });
    }
    info!("✓ Speclib URI: {}", speclib_uri);

    // Validate calib lib if provided
    let calib_lib_uri: Option<String> = args
        .calib_lib
        .as_ref()
        .map(|p| p.to_string_lossy().to_string());
    if let Some(ref uri) = calib_lib_uri {
        if !is_remote_uri(uri) && !std::path::Path::new(uri).exists() {
            return Err(errors::CliError::Io {
                source: "Calibration library file does not exist".to_string(),
                path: Some(uri.clone()),
            });
        }
        info!("✓ Calibration library URI: {}", uri);
    }

    // Validate all raw inputs: local-path existence only. Remote URIs are
    // resolved by tims_stage at staging time.
    for raw_uri in &raw_inputs {
        if !is_remote_uri(raw_uri) && !std::path::Path::new(raw_uri.as_str()).exists() {
            return Err(errors::CliError::Io {
                source: "Raw file does not exist".to_string(),
                path: Some(raw_uri.clone()),
            });
        }
    }
    info!("✓ All {} raw input(s) validated", raw_inputs.len());

    // Local-output path checks: existence and write probe. Remote outputs
    // (s3://, gs://, az://) skip these — the upload itself is the write test.
    if !is_remote_uri(&output_uri) {
        let output_dir_path = std::path::Path::new(&output_uri);

        if output_dir_path.exists() && !args.overwrite {
            return Err(errors::CliError::Config {
                source: format!(
                    "Output directory {:?} already exists. Use --overwrite/-O to overwrite.",
                    output_uri
                ),
            });
        }

        match std::fs::create_dir_all(output_dir_path) {
            Ok(_) => {
                if args.overwrite {
                    info!("✓ Using existing output directory (overwrite mode)");
                } else {
                    info!("✓ Created output directory");
                }
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

        // Validate per-file output subdirectories (local only)
        for raw_uri in &raw_inputs {
            let file_stem = sample_name_from_uri(raw_uri).ok_or_else(|| errors::CliError::Io {
                source: "Unable to extract file stem".to_string(),
                path: Some(raw_uri.clone()),
            })?;
            let file_output_dir = output_dir_path.join(&file_stem);

            if file_output_dir.exists() && !args.overwrite {
                return Err(errors::CliError::Config {
                    source: format!(
                        "Output subdirectory {:?} already exists. Use --overwrite/-O to overwrite.",
                        file_output_dir
                    ),
                });
            }
        }
        info!("✓ All output subdirectories validated");
    } else {
        info!(
            "✓ Remote output URI, skipping local path checks: {}",
            output_uri
        );
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

/// Extract a sample name from a raw-input URI. Handles trailing slashes and
/// `.d`/`.tar`/`.idx` suffixes so the result is stable across input kinds.
fn sample_name_from_uri(uri: &str) -> Option<String> {
    let trimmed = uri.trim_end_matches('/');
    let last = trimmed.rsplit('/').next()?;
    let stem = last
        .strip_suffix(".d")
        .or_else(|| last.strip_suffix(".tar"))
        .or_else(|| last.strip_suffix(".idx"))
        .unwrap_or(last);
    if stem.is_empty() {
        None
    } else {
        Some(stem.to_string())
    }
}

fn get_frag_range_from_index(
    index: &IndexedTimstofPeaks,
) -> Result<TupleRange<f64>, errors::CliError> {
    // `IndexedTimstofPeaks::fragmented_range` already folds over every
    // ms2 window-group via `QuadrupoleIsolationScheme::fragmented_range`,
    // which is defined on the ring-shape geometry (AABB/Trapezoid/Polygon).
    // It panics only if there are no window groups — treat that as a
    // non-DIA run and surface a readable error instead.
    //
    // NOTE: the task spec proposed reaching into raw `isolation_mz` /
    // `isolation_width` Vec<f32> fields, but `QuadrupoleIsolationScheme`
    // stores classified ring shapes (not the raw arrays). Using the
    // existing public accessor is both correct and narrower.
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
) -> std::result::Result<timsseek::scoring::PipelineReport, errors::CliError> {
    info!("Processing raw input: {}", raw_uri);

    let step = TimedStep::begin("Loading index");
    let index =
        load_index(raw_uri, backend, save_sidecar, CentroidingConfig::default()).map_err(|e| {
            errors::CliError::Io {
                source: format!("load_index({raw_uri}): {e}"),
                path: Some(raw_uri.to_string()),
            }
        })?;
    let load_index_ms = step.finish().as_millis() as u64;
    alloc_track::snap!("Loading index");

    // Rebucket to the scoring-optimal size. Existing on-disk caches
    // are written with `bucket_size=4096`, which is too large for the
    // tight mz tolerances used by Phase 1 / Phase 3 (measured: ~−24%
    // wall at bucket_size=256). `BUCKET_SIZE` env var overrides for
    // experiments.
    let new_bucket_size: usize = std::env::var("BUCKET_SIZE")
        .ok()
        .and_then(|s| s.parse().ok())
        .filter(|&bs| bs > 0)
        .unwrap_or(256);
    let step = TimedStep::begin(format_args!("Rebucket at {}", new_bucket_size));
    let index = index.rebucket(new_bucket_size);
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

    // If overwrite mode, delete the specific files we're about to write
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
    )
    .map_err(|e| errors::CliError::DataReading {
        source: format!("{}", e),
    })?;

    info!("Successfully processed {}", raw_uri);
    Ok(report)
}

/// Handles local-vs-remote output routing. When `dest_uri` is an s3:// /
/// gs:// / az:// URL, writes into a tempdir and uploads per-sample; when
/// local, writes directly.
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

    /// Upload and remove a per-sample subdir after the sample has finished
    /// writing; no-op for local destinations.
    fn finalize_sample(&self, sample: &str) -> Result<(), errors::CliError> {
        if !self.remote {
            return Ok(());
        }
        let local = self.sample_dir(sample);
        for entry in std::fs::read_dir(&local).map_err(|e| errors::CliError::Io {
            source: format!("read sample dir: {e}"),
            path: Some(local.to_string_lossy().to_string()),
        })? {
            let entry = entry.map_err(|e| errors::CliError::Io {
                source: format!("read dir entry: {e}"),
                path: None,
            })?;
            let bn = entry.file_name().to_string_lossy().to_string();
            let dest = format!("{}/{}/{}", self.dest_uri.trim_end_matches('/'), sample, bn);
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

/// Load a Speclib from a local path or remote URI. Remote URIs are
/// downloaded to a tempfile (preserving the original filename so
/// `SpeclibFormat::detect_from_path` still works on extension), then
/// handed to the existing path-based `Speclib::from_file`.
///
/// NOTE: the task spec suggested wiring `tims_stage::open_reader` directly
/// into the speclib loader, but `Speclib::from_file` takes a `Path` (it
/// sniffs format by extension and re-opens the file) rather than a
/// `Read`. Rerouting through a tempfile keeps the change scoped to
/// `timsseek_cli` and preserves format auto-detection.
fn speclib_from_uri(
    uri: &str,
    decoy_strategy: timsseek::DecoyStrategy,
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

    // Remote: stream into a tempfile preserving the filename so
    // format-by-extension detection still works.
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
    {
        use std::io::Write;
        let mut reader = tims_stage::open_reader(uri).map_err(|e| errors::CliError::Io {
            source: format!("open speclib {uri}: {e}"),
            path: Some(uri.to_string()),
        })?;
        let mut f = std::fs::File::create(&local).map_err(|e| errors::CliError::Io {
            source: format!("create tempfile: {e}"),
            path: Some(local.to_string_lossy().to_string()),
        })?;
        std::io::copy(&mut reader, &mut f).map_err(|e| errors::CliError::Io {
            source: format!("write tempfile: {e}"),
            path: Some(local.to_string_lossy().to_string()),
        })?;
        f.flush().map_err(|e| errors::CliError::Io {
            source: format!("flush tempfile: {e}"),
            path: Some(local.to_string_lossy().to_string()),
        })?;
    }
    let lib = timsseek::data_sources::speclib::Speclib::from_file(&local, decoy_strategy).map_err(
        |e| errors::CliError::Config {
            source: format!("Failed to load speclib {uri}: {:?}", e),
        },
    )?;
    Ok((lib, Some(td)))
}

fn main() {
    if let Err(e) = run() {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run() -> std::result::Result<(), errors::CliError> {
    // Parse command line arguments first
    let args = Cli::parse();

    // Short-circuit flags that dump the embedded default config.
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

    // Load and parse configuration, or use defaults
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

    // Override config with command line arguments if provided
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

    // Override decoy strategy if provided
    if let Some(strategy) = args.decoy_strategy {
        config.analysis.decoy_strategy = strategy;
    }

    // === Set up tracing subscriber ===
    // We defer this until after config/validation so we know the output directory for the log file.

    // Determine log file path. Remote output URIs skip the default file-log
    // path — we'd be writing to a non-local directory; use `--log-path` to
    // point at a local file explicitly in that case.
    let log_file_path = match args.log_path {
        Some(ref p) if p.to_str() == Some("-") => None, // stderr-only mode
        Some(ref p) => Some(p.clone()),
        None => args
            .output_uri
            .as_ref()
            .or(config.output.as_ref().map(|o| &o.uri))
            .filter(|d| !is_remote_uri(d.as_str()))
            .map(|d| std::path::Path::new(d).join("timsseek.log")),
    };

    // Build the env filter for the main logging layer
    let env_filter = EnvFilter::builder()
        .with_default_directive(
            args.log_level
                .parse()
                .unwrap_or_else(|_| "info".parse().unwrap()),
        )
        .from_env_lossy()
        .add_directive("forust_ml=warn".parse().unwrap())
        .add_directive("timscentroid::storage=warn".parse().unwrap());

    // Use Option layers so we can build a single subscriber type regardless
    // of whether we're writing to a log file or to stderr.
    let (file_layer, stderr_warn_layer, stderr_all_layer) =
        if let Some(ref log_path) = log_file_path {
            // File mode: log file gets env_filter, stderr gets WARN+ only
            if let Some(parent) = log_path.parent() {
                let _ = std::fs::create_dir_all(parent);
            }
            let log_file = std::fs::File::create(log_path).expect("Failed to create log file");
            let fl = fmt::layer()
                .with_writer(std::sync::Mutex::new(log_file))
                .with_filter(env_filter);
            let sl = fmt::layer()
                .with_writer(std::io::stderr)
                .without_time()
                .with_filter(tracing_subscriber::filter::LevelFilter::WARN);
            (Some(fl), Some(sl), None)
        } else {
            // stderr-only mode (--log-path -)
            let sl = fmt::layer()
                .with_writer(std::io::stderr)
                .with_filter(env_filter);
            (None, None, Some(sl))
        };

    #[cfg(feature = "instrumentation")]
    let perf_filter = EnvFilter::builder()
        .with_default_directive("trace".parse().unwrap())
        .with_env_var("RUST_PERF_LOG")
        .from_env_lossy()
        .add_directive("forust_ml::gradientbooster=warn".parse().unwrap());

    #[cfg(feature = "instrumentation")]
    let events_filter = tracing_subscriber::filter::filter_fn(|metadata| !metadata.is_event());

    #[cfg(feature = "instrumentation")]
    let (tree_layer, _guard) = PrintTreeLayer::new(PrintTreeConfig {
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
        .with(stderr_warn_layer)
        .with(stderr_all_layer);

    #[cfg(feature = "instrumentation")]
    let reg = reg.with(tree_layer);

    reg.init();

    // Print version and log path to stdout
    if let Some(ref log_path) = log_file_path {
        println!("timsseek v{}", env!("CARGO_PKG_VERSION"));
        println!("Log: {}", log_path.display());
        println!();
    }

    info!("Parsed configuration: {:#?}", config.clone());
    alloc_track::snap!("start");

    let validated = validate_inputs(&config, &args)?;

    // Build staging backend once from [staging] config (falls back to
    // defaults when absent). The sweep runs in `PerRunTempdir::new`.
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

    // Route outputs through a local tempdir when the destination is remote
    // (s3://, gs://, az://). Per-sample subdirs are uploaded + removed
    // after each sample finishes; run-level files are uploaded after the
    // batch.
    let sink = OutputSink::new(&validated.output_uri)?;

    let config_output_path = sink.root().join("config_used.json");

    // If overwrite mode, delete existing config file (local-only; remote
    // will simply overwrite on upload).
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

    // Load speclib once (shared across all files)
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
        .finish_with(format_args!("{} entries", speclib.len()))
        .as_millis() as u64;
    alloc_track::snap!("Loading speclib");

    // Load calibration library once (if provided)
    let (calib_lib, _calib_td, load_calib_lib_ms) = match &validated.calib_lib_uri {
        Some(uri) => {
            let step = TimedStep::begin("Loading calib lib");
            info!("Loading calibration library from {}", uri);
            let (lib, td) = speclib_from_uri(uri, config.analysis.decoy_strategy)?;
            let ms = step
                .finish_with(format_args!("{} entries", lib.len()))
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
        ) {
            Ok(report) => {
                if let Err(e) = sink.finalize_sample(&sample_name) {
                    error!("Failed to finalize sample {}: {}", sample_name, e);
                    failed_files.push((raw_uri.clone(), e));
                    error!("Aborting batch due to upload failure");
                    break;
                }
                successful_files.push(raw_uri.clone());
                run_report.files.push(timsseek::scoring::FileReport {
                    file_name: raw_uri.clone(),
                    pipeline: report,
                });
            }
            Err(e) => {
                error!("Failed to process {}: {}", raw_uri, e);
                // I/O errors are likely systemic (disk full, permissions) —
                // abort the batch instead of failing every remaining file.
                if matches!(e, errors::CliError::Io { .. }) {
                    failed_files.push((raw_uri.clone(), e));
                    error!("Aborting batch due to I/O error");
                    break;
                }
                failed_files.push((raw_uri.clone(), e));
            }
        }
    }

    // Write run-level report into the sink's working dir.
    let run_report_path = sink.root().join("run_report.json");
    if let Ok(json) = serde_json::to_string_pretty(&run_report) {
        let _ = std::fs::write(&run_report_path, json);
        info!("Wrote run report to {:?}", run_report_path);
    }

    // Upload run-level artifacts for remote destinations (no-op locally).
    sink.finalize_run(&["run_report.json", "config_used.json"])?;

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
