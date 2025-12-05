/// This file is mainly used for testing and demonstration purposes.
///
/// In particular is used to check that the code to generate
/// spectral libs from python works correctly in Rust.
use clap::{
    Parser,
    Subcommand,
};
// use serde_json::to_string_pretty;
use timsseek::Speclib;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: SubCommands,
}

#[derive(Subcommand)]
enum SubCommands {
    // Sample,
    Parse {
        #[arg(short, long)]
        speclib_file: String,
    },
}

// fn generate_sample_speclib() {
//     // This function generates a sample speclib and prints it to stdout.
//     let speclib = Speclib::sample();
//     println!("{}", to_string_pretty(&speclib).unwrap());
// }

fn parse_speclib(speclib_file: &str) -> Result<(), timsseek::errors::LibraryReadingError> {
    // This function parses a speclib file and prints the parsed content to stdout.
    let speclib = Speclib::from_file(std::path::Path::new(speclib_file))?;
    println!("{:#?}", speclib);
    Ok(())
}

fn main() -> Result<(), timsseek::errors::LibraryReadingError> {
    // This is a placeholder for the main function.
    // The actual implementation would depend on the specific requirements of the application.
    let cli = Cli::parse();

    match &cli.command {
        // SubCommands::Sample => {
        //     // Generate a sample speclib
        //     generate_sample_speclib();
        //     Ok(())
        // }
        SubCommands::Parse { speclib_file } => {
            // Parse the provided speclib file
            parse_speclib(speclib_file)?;
            Ok(())
        }
    }
}
