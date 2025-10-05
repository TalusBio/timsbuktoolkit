use serde::{
    Deserialize,
    Serialize,
};
use std::path::PathBuf;
use timsquery::Tolerance;

use crate::cli::Cli;
use crate::errors;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Config {
    pub input: Option<InputConfig>,
    pub analysis: AnalysisConfig,
    pub output: Option<OutputConfig>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(tag = "type")]
pub enum InputConfig {
    // TODO: Implement in-rust speclib generation to take a fasta file.
    // #[serde(rename = "fasta")]
    // Fasta {
    //     path: PathBuf,
    //     digestion: DigestionConfig,
    // },
    #[serde(rename = "speclib")]
    Speclib { path: PathBuf },
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AnalysisConfig {
    pub dotd_file: Option<PathBuf>,
    pub chunk_size: usize,
    pub tolerance: Tolerance,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct OutputConfig {
    pub directory: PathBuf,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DigestionConfig {
    pub min_length: u32,
    pub max_length: u32,
    pub max_missed_cleavages: u32,
    pub build_decoys: bool,
}

impl Default for DigestionConfig {
    fn default() -> Self {
        Self {
            min_length: 7,
            max_length: 30,
            max_missed_cleavages: 2,
            build_decoys: true,
        }
    }
}

impl Config {
    #![allow(dead_code)]
    pub fn with_cli_args(_config: Cli) -> Result<Self, errors::CliError> {
        unimplemented!()
    }
}
