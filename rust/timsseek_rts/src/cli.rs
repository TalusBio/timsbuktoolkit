use clap::Parser;
use std::path::PathBuf;
use timsquery::Tolerance;
use timsseek::errors::{
    Result,
    TimsSeekError,
};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Path to the .d file (will over-write the config file)
    #[arg(short, long)]
    pub dotd_file: PathBuf,

    /// Path to the speclib file
    #[arg(short, long)]
    pub config: PathBuf,

    /// Path to the output directory
    #[arg(short, long)]
    #[clap(default_value("127.0.0.1:3724"))]
    pub address: String,
}

impl Cli {
    pub fn read_config(&self) -> Result<Tolerance> {
        let conf = match std::fs::File::open(&self.config) {
            Ok(x) => x,
            Err(e) => {
                return Err(TimsSeekError::Io {
                    source: e,
                    path: None,
                });
            }
        };
        let config: core::result::Result<Tolerance, _> = serde_json::from_reader(conf);
        let config = match config {
            Ok(x) => x,
            Err(e) => {
                return Err(TimsSeekError::ParseError { msg: e.to_string() });
            }
        };
        Ok(config)
    }
}
