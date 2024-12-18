use serde::{
    Deserialize,
    Serialize,
};
use std::path::PathBuf;
use timsquery::traits::tolerance::DefaultTolerance;

#[derive(Debug, Serialize, Deserialize)]
pub struct Config {
    /// Input configuration
    pub input: InputConfig,

    /// Analysis parameters
    pub analysis: AnalysisConfig,

    /// Output configuration
    pub output: OutputConfig,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum InputConfig {
    #[serde(rename = "fasta")]
    Fasta {
        path: PathBuf,
        digestion: DigestionConfig,
    },
    #[serde(rename = "speclib")]
    Speclib { path: PathBuf },
}

#[derive(Debug, Serialize, Deserialize)]
pub struct AnalysisConfig {
    /// Path to the .d file
    pub dotd_file: Option<PathBuf>,

    /// Processing parameters
    pub chunk_size: usize,

    /// Tolerance settings
    pub tolerance: DefaultTolerance,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct OutputConfig {
    /// Directory for results
    pub directory: PathBuf,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct DigestionConfig {
    pub min_length: u32,
    pub max_length: u32,
    pub max_missed_cleavages: u32,
    pub build_decoys: bool,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ToleranceConfig {
    pub ms_ppm: (f64, f64),
    pub mobility_pct: (f64, f64),
    pub quad_absolute: (f64, f64),
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

impl Default for ToleranceConfig {
    fn default() -> Self {
        Self {
            ms_ppm: (10.0, 20.0),
            mobility_pct: (0.015, 0.025),
            quad_absolute: (0.05, 0.1),
        }
    }
}
