use clap::Parser;
use timsquery::models::aggregators::MultiCMGStatsFactory;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsseek::cli::Cli;
use timsseek::config::{
    Config,
    InputConfig,
};
use timsseek::errors::TimsSeekError;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;
use timsseek::processing;

fn main() -> std::result::Result<(), TimsSeekError> {
    // Initialize logging
    env_logger::init();

    // Parse command line arguments
    let args = Cli::parse();

    // Load and parse configuration
    let config: Result<Config, _> = serde_json::from_reader(std::fs::File::open(args.config)?);
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
        config.output.directory = output_dir;
    }

    println!("{:?}", config);

    // Create output directory
    std::fs::create_dir_all(&config.output.directory)?;

    let dotd_file_location = &config.analysis.dotd_file;
    let index = QuadSplittedTransposedIndex::from_path_centroided(
        dotd_file_location
            .clone()
            .unwrap() // TODO: Error handling
            .to_str()
            .expect("Path is not convertable to string"),
    )?;

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
                &factory,
                digestion,
                &config.analysis,
                &config.output,
            )?;
        }
        InputConfig::Speclib { path } => {
            processing::process_speclib(path, &index, &factory, &config.analysis, &config.output)?;
        }
    }

    Ok(())
}
