mod cli;
mod config;
mod errors;
mod processing;

use clap::Parser;
use timsquery::TimsTofPath;
use timsquery::serde::load_index_auto;
use timsquery::utils::TupleRange;
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

#[cfg(target_os = "windows")]
use mimalloc::MiMalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

/// Validated inputs ready for processing
struct ValidatedInputs {
    dotd_files: Vec<std::path::PathBuf>,
    speclib_path: std::path::PathBuf,
    calib_lib_path: Option<std::path::PathBuf>,
    output_directory: std::path::PathBuf,
    overwrite: bool,
}

/// Validates all inputs and outputs before processing begins.
/// Returns ValidatedInputs on success, or an error if any validation fails.
fn validate_inputs(
    config: &Config,
    args: &Cli,
) -> std::result::Result<ValidatedInputs, errors::CliError> {
    info!("Validating inputs and outputs before processing...");

    // Get list of raw files to process
    let dotd_files = match config.analysis.dotd_files.clone() {
        Some(files) => files,
        None => {
            return Err(errors::CliError::Config {
                source: "No raw files provided, please provide dotd_files in either the config file or with the --dotd-files flag".to_string(),
            });
        }
    };

    // Get speclib path
    let speclib_path = match &config.input {
        Some(InputConfig::Speclib { path }) => path.clone(),
        None => {
            return Err(errors::CliError::Config {
                source: "No input specified".to_string(),
            });
        }
    };

    // Get output directory
    let output_directory = match &config.output {
        Some(output_config) => output_config.directory.clone(),
        None => {
            return Err(errors::CliError::Config {
                source: "No output directory specified".to_string(),
            });
        }
    };

    // Validate speclib exists
    if !speclib_path.exists() {
        return Err(errors::CliError::Io {
            source: "Speclib file does not exist".to_string(),
            path: Some(speclib_path.to_string_lossy().to_string()),
        });
    }
    info!("✓ Speclib file exists: {:?}", speclib_path);

    // Validate calib lib if provided
    let calib_lib_path = args.calib_lib.clone();
    if let Some(ref path) = calib_lib_path {
        if !path.exists() {
            return Err(errors::CliError::Io {
                source: "Calibration library file does not exist".to_string(),
                path: Some(path.to_string_lossy().to_string()),
            });
        }
        info!("✓ Calibration library exists: {:?}", path);
    }

    // Validate all raw files exist
    for dotd_file in &dotd_files {
        if !dotd_file.exists() {
            return Err(errors::CliError::Io {
                source: "Raw file does not exist".to_string(),
                path: Some(dotd_file.to_string_lossy().to_string()),
            });
        }
    }
    info!("✓ All {} raw file(s) exist", dotd_files.len());

    // Check if output directory exists and handle based on overwrite flag
    if output_directory.exists() && !args.overwrite {
        return Err(errors::CliError::Config {
            source: format!(
                "Output directory {:?} already exists. Use --overwrite/-O to overwrite.",
                output_directory
            ),
        });
    }

    // Create output directory and test write permissions
    match std::fs::create_dir_all(&output_directory) {
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
                path: Some(output_directory.to_string_lossy().to_string()),
            });
        }
    };

    // Test write permissions by creating a test file
    let test_file = output_directory.join(".write_test");
    match std::fs::write(&test_file, "test") {
        Ok(_) => {
            let _ = std::fs::remove_file(&test_file);
            info!("✓ Output directory is writable");
        }
        Err(e) => {
            return Err(errors::CliError::Io {
                source: format!("Output directory is not writable: {}", e),
                path: Some(output_directory.to_string_lossy().to_string()),
            });
        }
    }

    // Validate per-file output subdirectories
    for dotd_file in &dotd_files {
        let file_stem = dotd_file
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| errors::CliError::Io {
                source: "Unable to extract file stem".to_string(),
                path: Some(dotd_file.to_string_lossy().to_string()),
            })?;
        let file_output_dir = output_directory.join(file_stem);

        // Check if per-file directory exists when not in overwrite mode
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

    info!("All validations passed! Starting processing...");

    Ok(ValidatedInputs {
        dotd_files,
        speclib_path,
        calib_lib_path,
        output_directory,
        overwrite: args.overwrite,
    })
}

fn get_frag_range(file: &TimsTofPath) -> Result<TupleRange<f64>, errors::CliError> {
    let reader = file
        .load_frame_reader()
        .map_err(|e| errors::CliError::DataReading {
            source: format!("Failed to load frame reader: {:?}", e),
        })?;
    let dia_windows = reader
        .dia_windows
        .as_ref()
        .ok_or_else(|| errors::CliError::DataReading {
            source: "File does not contain DIA windows — is this a DIA run?".to_string(),
        })?;

    let upper_mz = dia_windows
        .iter()
        .map(|w| {
            w.isolation_mz
                .iter()
                .zip(w.isolation_width.iter())
                .map(|(imz, iw)| imz + (iw / 2.0))
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max);

    let lower_mz = dia_windows
        .iter()
        .map(|w| {
            w.isolation_mz
                .iter()
                .zip(w.isolation_width.iter())
                .map(|(imz, iw)| imz - (iw / 2.0))
                .fold(f64::MAX, f64::min)
        })
        .fold(f64::MAX, f64::min);

    TupleRange::try_new(lower_mz, upper_mz).map_err(|e| errors::CliError::DataReading {
        source: format!("Invalid DIA m/z range: {:?}", e),
    })
}

fn process_single_file(
    dotd_file: &std::path::Path,
    speclib: &timsseek::data_sources::speclib::Speclib,
    calib_lib: Option<&timsseek::data_sources::speclib::Speclib>,
    config: &Config,
    base_output_dir: &std::path::Path,
    overwrite: bool,
    max_qvalue: f32,
) -> std::result::Result<timsseek::scoring::PipelineReport, errors::CliError> {
    info!("Processing file: {:?}", dotd_file);

    let timstofpath =
        TimsTofPath::new(dotd_file.to_str().unwrap()).map_err(|e| errors::CliError::Io {
            source: format!("Failed to open raw file: {:?}", e),
            path: Some(dotd_file.to_string_lossy().to_string()),
        })?;

    let step = TimedStep::begin("Loading index");
    let index = load_index_auto(
        dotd_file.to_str().ok_or_else(|| errors::CliError::Io {
            source: "Invalid path encoding".to_string(),
            path: None,
        })?,
        None,
    )?
    .into_eager()?;
    let load_index_ms = step.finish().as_millis() as u64;

    let fragmented_range = get_frag_range(&timstofpath)?;

    let pipeline = Scorer {
        index,
        broad_tolerance: config.analysis.tolerance.clone(),
        fragmented_range,
    };

    let file_stem = dotd_file
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| errors::CliError::Io {
            source: "Unable to extract file stem".to_string(),
            path: Some(dotd_file.to_string_lossy().to_string()),
        })?;
    let file_output_dir = base_output_dir.join(file_stem);

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
        directory: file_output_dir,
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

    info!("Successfully processed {:?}", dotd_file);
    Ok(report)
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
    if !args.dotd_files.is_empty() {
        config.analysis.dotd_files = Some(args.dotd_files.clone());
    }
    if let Some(ref speclib_file) = args.speclib_file {
        config.input = Some(InputConfig::Speclib {
            path: speclib_file.clone(),
        });
    }
    if config.input.is_none() {
        return Err(errors::CliError::Config {
            source: "No input provided, please provide one in either the config file or with the --speclib-file flag".to_string(),
        });
    }
    if let Some(ref output_dir) = args.output_dir {
        config.output = Some(OutputConfig {
            directory: output_dir.clone(),
        });
    }

    // Override decoy strategy if provided
    if let Some(strategy) = args.decoy_strategy {
        config.analysis.decoy_strategy = strategy;
    }

    // === Set up tracing subscriber ===
    // We defer this until after config/validation so we know the output directory for the log file.

    // Determine log file path
    let log_file_path = match args.log_path {
        Some(ref p) if p.to_str() == Some("-") => None, // stderr-only mode
        Some(ref p) => Some(p.clone()),
        None => args
            .output_dir
            .as_ref()
            .or(config.output.as_ref().map(|o| &o.directory))
            .map(|d| d.join("timsseek.log")),
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

    let validated = validate_inputs(&config, &args)?;

    let config_output_path = validated.output_directory.join("config_used.json");

    // If overwrite mode, delete existing config file
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
    let mut failed_files: Vec<(std::path::PathBuf, errors::CliError)> = Vec::new();
    let mut successful_files: Vec<std::path::PathBuf> = Vec::new();

    // Load speclib once (shared across all files)
    let step = TimedStep::begin("Loading speclib");
    info!(
        "Building database from speclib file {:?}",
        validated.speclib_path
    );
    info!(
        "Decoy generation strategy: {}",
        config.analysis.decoy_strategy
    );
    let speclib = timsseek::data_sources::speclib::Speclib::from_file(
        &validated.speclib_path,
        config.analysis.decoy_strategy,
    )
    .map_err(|e| errors::CliError::Config {
        source: format!("Failed to load speclib: {:?}", e),
    })?;
    let load_speclib_ms = step
        .finish_with(format_args!("{} entries", speclib.len()))
        .as_millis() as u64;

    // Load calibration library once (if provided)
    let (calib_lib, load_calib_lib_ms) = match &validated.calib_lib_path {
        Some(p) => {
            let step = TimedStep::begin("Loading calib lib");
            info!("Loading calibration library from {:?}", p);
            let lib = timsseek::data_sources::speclib::Speclib::from_file(
                p,
                config.analysis.decoy_strategy,
            )
            .map_err(|e| errors::CliError::Config {
                source: format!("Failed to load calib lib: {:?}", e),
            })?;
            let ms = step
                .finish_with(format_args!("{} entries", lib.len()))
                .as_millis() as u64;
            (Some(lib), ms)
        }
        None => (None, 0),
    };

    run_report.load_speclib_ms = load_speclib_ms;
    run_report.load_calib_lib_ms = load_calib_lib_ms;
    run_report.speclib_entries = speclib.len();
    run_report.calib_lib_entries = calib_lib.as_ref().map_or(0, |l| l.len());

    let total_files = validated.dotd_files.len();
    info!("Processing {} raw file(s)", total_files);

    for (idx, dotd_file) in validated.dotd_files.iter().enumerate() {
        info!(
            "Processing file {} of {}: {:?}",
            idx + 1,
            total_files,
            dotd_file
        );

        match process_single_file(
            dotd_file,
            &speclib,
            calib_lib.as_ref(),
            &config,
            &validated.output_directory,
            validated.overwrite,
            args.max_qvalue,
        ) {
            Ok(report) => {
                successful_files.push(dotd_file.clone());
                run_report.files.push(timsseek::scoring::FileReport {
                    file_name: dotd_file.to_string_lossy().to_string(),
                    pipeline: report,
                });
            }
            Err(e) => {
                error!("Failed to process {:?}: {}", dotd_file, e);
                // I/O errors are likely systemic (disk full, permissions) —
                // abort the batch instead of failing every remaining file.
                if matches!(e, errors::CliError::Io { .. }) {
                    failed_files.push((dotd_file.clone(), e));
                    error!("Aborting batch due to I/O error");
                    break;
                }
                failed_files.push((dotd_file.clone(), e));
            }
        }
    }

    // Write run-level report
    let run_report_path = validated.output_directory.join("run_report.json");
    if let Ok(json) = serde_json::to_string_pretty(&run_report) {
        let _ = std::fs::write(&run_report_path, json);
        info!("Wrote run report to {:?}", run_report_path);
    }

    info!("Successfully processed {} file(s)", successful_files.len());
    if !failed_files.is_empty() {
        error!("Failed to process {} file(s):", failed_files.len());
        for (file, err) in &failed_files {
            error!("  {:?}: {}", file, err);
        }
        return Err(errors::CliError::Config {
            source: format!("Failed to process {} file(s)", failed_files.len()),
        });
    }

    Ok(())
}
