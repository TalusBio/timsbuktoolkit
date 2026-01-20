mod cli;
mod config;
mod errors;
mod processing;

use clap::Parser;
use timsquery::TimsTofPath;
use timsquery::models::tolerance::RtTolerance;
use timsquery::serde::load_index_auto;
use timsquery::utils::TupleRange;
use timsseek::scoring::{
    ScoringPipeline,
    ToleranceHierarchy,
};
use tracing::{
    error,
    info,
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

#[cfg(target_os = "windows")]
use mimalloc::MiMalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

/// Validated inputs ready for processing
struct ValidatedInputs {
    dotd_files: Vec<std::path::PathBuf>,
    speclib_path: std::path::PathBuf,
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
        output_directory,
        overwrite: args.overwrite,
    })
}

fn get_frag_range(file: &TimsTofPath) -> TupleRange<f64> {
    let reader = file.load_frame_reader().unwrap();
    let upper_mz = reader
        .dia_windows
        .as_ref()
        .expect("DIA windows should be present for a dia run")
        .iter()
        .map(|w| {
            w.isolation_mz
                .iter()
                .zip(w.isolation_width.iter())
                .map(|(imz, iw)| imz + (iw / 2.0))
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max);

    let lower_mz = reader
        .dia_windows
        .expect("DIA windows should be present for a dia run")
        .iter()
        .map(|w| {
            w.isolation_mz
                .iter()
                .zip(w.isolation_width.iter())
                .map(|(imz, iw)| imz - (iw / 2.0))
                .fold(f64::MAX, f64::min)
        })
        .fold(f64::MAX, f64::min);
    TupleRange::try_new(lower_mz, upper_mz).unwrap()
}

fn process_single_file(
    dotd_file: &std::path::Path,
    speclib_path: &std::path::Path,
    config: &Config,
    base_output_dir: &std::path::Path,
    overwrite: bool,
) -> std::result::Result<(), errors::CliError> {
    info!("Processing file: {:?}", dotd_file);

    let timstofpath =
        TimsTofPath::new(dotd_file.to_str().unwrap()).map_err(|e| errors::CliError::Io {
            source: format!("Failed to open raw file: {:?}", e),
            path: Some(dotd_file.to_string_lossy().to_string()),
        })?;

    let index = load_index_auto(
        dotd_file.to_str().ok_or_else(|| errors::CliError::Io {
            source: "Invalid path encoding".to_string(),
            path: None,
        })?,
        None,
    )?
    .into_eager()?;

    let fragmented_range = get_frag_range(&timstofpath);

    let pipeline = ScoringPipeline {
        index,
        tolerances: ToleranceHierarchy {
            prescore: config.analysis.tolerance.clone(),
            secondary: config
                .analysis
                .tolerance
                .clone()
                .with_rt_tolerance(RtTolerance::Minutes((0.2, 0.2))),
        },
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

    // Process speclib
    processing::process_speclib(
        speclib_path,
        &pipeline,
        config.analysis.chunk_size,
        &file_output_config,
        config.analysis.decoy_strategy,
    )
    .unwrap();

    info!("Successfully processed {:?}", dotd_file);
    Ok(())
}

/// Converts verbosity flags to a log level string.
/// Returns the log level based on verbose/quiet counts.
/// If RUST_LOG is set, it takes precedence.
fn get_log_level(verbose: u8, quiet: u8) -> String {
    // RUST_LOG environment variable takes precedence
    if std::env::var("RUST_LOG").is_ok() {
        return std::env::var("RUST_LOG").unwrap();
    }

    // Calculate effective verbosity: positive = more verbose, negative = more quiet
    let effective = verbose as i8 - quiet as i8;

    match effective {
        2.. => "trace".to_string(),
        1 => "debug".to_string(),
        0 => "info".to_string(),
        -1 => "warn".to_string(),
        _ => "error".to_string(),
    }
}

fn main() -> std::result::Result<(), errors::CliError> {
    // Parse command line arguments first to get verbosity flags
    let args = Cli::parse();

    let log_level = get_log_level(args.verbose, args.quiet);
    let fmt_filter = EnvFilter::builder()
        .with_default_directive(log_level.parse().unwrap())
        .from_env_lossy();

    #[cfg(feature = "instrumentation")]
    let perf_filter = EnvFilter::builder()
        .with_default_directive("trace".parse().unwrap())
        .with_env_var("RUST_PERF_LOG")
        .from_env_lossy()
        .add_directive("forust_ml::gradientbooster=warn".parse().unwrap());

    // Filter out events but keep spans
    #[cfg(feature = "instrumentation")]
    let events_filter = tracing_subscriber::filter::filter_fn(|metadata| !metadata.is_event());

    // I am aware that this conditional compilation is ugly ...
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

    // let (pf_layer, pf_guard) = PerfettoLayer::new_from_env().unwrap();

    let fmt_layer = fmt::layer()
        .with_span_events(FmtSpan::CLOSE)
        .with_filter(fmt_filter);

    let reg = tracing_subscriber::registry().with(fmt_layer);

    #[cfg(feature = "instrumentation")]
    let reg = reg.with(tree_layer);

    reg.init();

    // Load and parse configuration, or use defaults
    let mut config = match args.config {
        Some(ref config_path) => {
            let conf = match std::fs::File::open(config_path) {
                Ok(x) => x,
                Err(e) => {
                    return Err(errors::CliError::Io {
                        source: e.to_string(),
                        path: Some(config_path.to_string_lossy().to_string()),
                    });
                }
            };
            let config: Result<Config, _> = serde_json::from_reader(conf);
            match config {
                Ok(x) => x,
                Err(e) => {
                    return Err(errors::CliError::ParseError { msg: e.to_string() });
                }
            }
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

    let mut failed_files: Vec<(std::path::PathBuf, errors::CliError)> = Vec::new();
    let mut successful_files: Vec<std::path::PathBuf> = Vec::new();

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
            &validated.speclib_path,
            &config,
            &validated.output_directory,
            validated.overwrite,
        ) {
            Ok(_) => {
                successful_files.push(dotd_file.clone());
            }
            Err(e) => {
                error!("Failed to process {:?}: {}", dotd_file, e);
                failed_files.push((dotd_file.clone(), e));
            }
        }
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
