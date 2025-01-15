mod cli;
mod config;
mod processing;

use clap::Parser;
use timsquery::models::aggregators::MultiCMGStatsFactory;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsseek::errors::TimsSeekError;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;
use timsseek::utils::tdf::get_ms1_frame_times_ms;
use tracing::level_filters::LevelFilter;
use tracing_subscriber::EnvFilter;

use cli::Cli;
use config::{
    Config,
    InputConfig,
    OutputConfig,
};

fn main() -> std::result::Result<(), TimsSeekError> {
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
    let conf = match std::fs::File::open(args.config) {
        Ok(x) => x,
        Err(e) => {
            return Err(TimsSeekError::Io {
                source: e,
                path: None,
            });
        }
    };
    let config: Result<Config, _> = serde_json::from_reader(conf);
    let mut config = match config {
        Ok(x) => x,
        Err(e) => {
            return Err(TimsSeekError::ParseError { msg: e.to_string() });
        }
    };

    // Override config with command line arguments if provided
    if let Some(dotd_file) = args.dotd_file {
        config.analysis.dotd_file = Some(dotd_file);
    }
    if let Some(speclib_file) = args.speclib_file {
        config.input = InputConfig::Speclib { path: speclib_file };
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
            return Err(TimsSeekError::Io {
                source: e,
                path: None,
            });
        }
    };

    let dotd_file_location = &config.analysis.dotd_file;
    let index = QuadSplittedTransposedIndex::from_path_centroided(
        dotd_file_location
            .clone()
            .unwrap() // TODO: Error handling
            .to_str()
            .expect("Path is not convertable to string"),
    )?;

    let tdf_path = &dotd_file_location.clone().unwrap().join("analysis.tdf");
    let ref_time_ms = get_ms1_frame_times_ms(tdf_path.to_str().unwrap()).unwrap();

    let factory = MultiCMGStatsFactory {
        converters: (index.mz_converter, index.im_converter),
        _phantom: std::marker::PhantomData::<SafePosition>,
    };

    // Process based on input type
    match config.input {
        InputConfig::Fasta { path, digestion } => {
            processing::process_fasta(
                path,
                &index,
                ref_time_ms.clone(),
                &factory,
                digestion,
                &config.analysis,
                &output_config,
            )?;
        }
        InputConfig::Speclib { path } => {
            processing::process_speclib(
                path,
                &index,
                ref_time_ms.clone(),
                &factory,
                &config.analysis,
                &output_config,
            )?;
        }
    }

    Ok(())
}
