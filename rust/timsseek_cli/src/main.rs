mod cli;
mod config;
mod errors;
mod processing;

use clap::Parser;
use timsquery::models::aggregators::EGCAggregator;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsseek::fragment_mass::IonAnnot;
use timsseek::utils::tdf::get_ms1_frame_times_ms;
use tracing::level_filters::LevelFilter;
use tracing_subscriber::EnvFilter;

use cli::Cli;
use config::{
    Config,
    InputConfig,
    OutputConfig,
};

// TODO: Check if this is faster in linux
// use mimalloc::MiMalloc;
//
// #[global_allocator]
// static GLOBAL: MiMalloc = MiMalloc;

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
            full_output: args.full_output,
        });
    }

    let mut output_config = match config.output {
        Some(ref x) => x.clone(),
        None => {
            panic!(
                "No output directory provided, please provide one in either the config file or with the --output-dir flag"
            );
        }
    };
    println!("{:#?}", config.clone());
    if args.full_output {
        output_config.full_output = true;
    }

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

    let tdf_path = &dotd_file_location.clone().unwrap().join("analysis.tdf");
    let ref_time_ms = get_ms1_frame_times_ms(tdf_path.to_str().unwrap()).unwrap();

    // Process based on input type
    match config.input {
        Some(InputConfig::Fasta { path, digestion }) => {
            processing::process_fasta(
                path,
                &index,
                digestion,
                &config.analysis,
                &output_config,
            )
            .unwrap();
        }
        Some(InputConfig::Speclib { path }) => {
            processing::process_speclib(
                path,
                &index,
                &config.analysis,
                &output_config,
            )
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
