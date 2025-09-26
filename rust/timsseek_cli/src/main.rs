mod cli;
mod config;
mod errors;
mod processing;

use clap::Parser;
use serde::Serialize;
use std::fs::File;
use timsquery::models::tolerance::RtTolerance;
use timsquery::utils::TupleRange;
use timsquery::{
    CentroidingConfig,
    IndexedTimstofPeaks,
    TimsTofPath,
};
use timsseek::scoring::Scorer;
use timsseek::utils::serde::load_index_caching;
use tracing::level_filters::LevelFilter;
use tracing::{
    error,
    info,
};
use tracing_subscriber::EnvFilter;

use cli::Cli;
use config::{
    Config,
    InputConfig,
    OutputConfig,
};
use std::sync::Arc;

#[cfg(target_os = "windows")]
use mimalloc::MiMalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn get_ms1_rts_as_millis(file: &TimsTofPath) -> Arc<[u32]> {
    let reader = file.load_frame_reader().unwrap();
    let mut rts: Vec<_> = reader
        .frame_metas
        .iter()
        .filter_map(|f| match f.ms_level {
            timsrust::MSLevel::MS1 => Some((f.rt_in_seconds * 1000.0).round() as u32),
            _ => None,
        })
        .collect();
    rts.sort_unstable();
    rts.dedup();
    rts.into()
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

fn main() -> std::result::Result<(), errors::CliError> {
    // Initialize logging
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::builder()
                .with_default_directive(LevelFilter::INFO.into())
                .from_env_lossy(),
        ) // This uses RUST_LOG environment variable
        .init();

    // Parse command line arguments
    let args = Cli::parse();

    // Load and parse configuration
    let conf = match std::fs::File::open(args.config.clone()) {
        Ok(x) => x,
        Err(e) => {
            return Err(errors::CliError::Io {
                source: e.to_string(),
                path: Some(args.config.to_string_lossy().to_string()),
            });
        }
    };
    let config: Result<Config, _> = serde_json::from_reader(conf);
    let mut config = match config {
        Ok(x) => x,
        Err(e) => {
            return Err(errors::CliError::ParseError { msg: e.to_string() });
        }
    };

    // Override config with command line arguments if provided
    if let Some(dotd_file) = args.dotd_file {
        config.analysis.dotd_file = Some(dotd_file);
    }
    if let Some(speclib_file) = args.speclib_file {
        config.input = Some(InputConfig::Speclib { path: speclib_file });
    }
    if config.input.is_none() {
        return Err(errors::CliError::Config {
            source: "No input provided, please provide one in either the config file or with the --speclib-file flag".to_string(),
        });
    }
    if let Some(output_dir) = args.output_dir {
        config.output = Some(OutputConfig {
            directory: output_dir,
        });
    }

    let output_config = match config.output {
        Some(ref x) => x.clone(),
        None => {
            panic!(
                "No output directory provided, please provide one in either the config file or with the --output-dir flag"
            );
        }
    };
    info!("Parsed configuration: {:#?}", config.clone());

    // Create output director
    match std::fs::create_dir_all(&output_config.directory) {
        Ok(_) => println!("Created output directory"),
        Err(e) => {
            return Err(errors::CliError::Io {
                source: e.to_string(),
                path: Some(output_config.directory.to_string_lossy().to_string()),
            });
        }
    };

    let file_loc = config.analysis.dotd_file.clone().unwrap();
    let timstofpath = match TimsTofPath::new(file_loc.to_str().unwrap()) {
        Ok(x) => x,
        Err(x) => {
            error!("Unable to find the file at path: {:?}", file_loc);
            error!("{:?}", x);
            panic!();
        }
    };

    let index = load_index_caching(file_loc).unwrap();
    let index_cycle_rt_ms = get_ms1_rts_as_millis(&timstofpath);
    let fragmented_range = get_frag_range(&timstofpath);

    // Process based on input type
    match config.input {
        // Some(InputConfig::Fasta { path, digestion }) => {
        //     // I like the idea of just converting the fasta into a speclib ...
        //     // but a part of me feels like that would be a waste of memoory ...
        //     // but it might also be a premature optmization ...
        //     // processing::process_fasta(path, &index, digestion, &config.analysis, &output_config)
        //     //     .unwrap();
        // }
        Some(InputConfig::Speclib { path }) => {
            let scorer = Scorer {
                index_cycle_rt_ms,
                index,
                tolerance: config.analysis.tolerance.clone(),
                secondary_tolerance: config
                    .analysis
                    .tolerance
                    .with_rt_tolerance(RtTolerance::Minutes((0.5, 0.5))),
                fragmented_range,
            };
            processing::process_speclib(path, &scorer, config.analysis.chunk_size, &output_config)
                .unwrap();
        }
        None => {
            return Err(errors::CliError::Config {
                source: "No input specified".to_string(),
            });
        }
    }

    Ok(())
}
