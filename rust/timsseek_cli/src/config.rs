use serde::{
    Deserialize,
    Serialize,
};
use std::path::PathBuf;
use timsquery::Tolerance;

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
