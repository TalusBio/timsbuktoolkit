mod cli;
mod config;
mod errors;
mod processing;

use clap::Parser;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsquery::models::tolerance::RtTolerance;
use timsseek::scoring::Scorer;
use tracing::level_filters::LevelFilter;
use tracing_subscriber::EnvFilter;

use cli::Cli;
use config::{
    Config,
    InputConfig,
    OutputConfig,
};

#[cfg(target_os = "windows")]
use mimalloc::MiMalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

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
    println!("{:#?}", config.clone());

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

    let dotd_file_location = &config.analysis.dotd_file;
    let index = QuadSplittedTransposedIndex::from_path_centroided(
        // let index = QuadSplittedTransposedIndex::from_path(
        dotd_file_location
            .clone()
            .unwrap() // TODO: Error handling
            .to_str()
            .expect("Path is not convertable to string"),
    )
    .unwrap();

    let fragmented_range = index.fragmented_range();

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
                index_cycle_rt_ms: index.cycle_rt_ms.clone(),
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
