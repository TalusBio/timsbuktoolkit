use serde::{
    Deserialize,
    Serialize,
};
use std::path::PathBuf;
use timsquery::Tolerance;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Config {
    pub input: Option<InputConfig>,
    pub analysis: AnalysisConfig,
    pub output: Option<OutputConfig>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(tag = "type")]
pub enum InputConfig {
    #[serde(rename = "speclib")]
    Speclib { path: PathBuf },
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AnalysisConfig {
    pub dotd_files: Option<Vec<PathBuf>>,
    pub chunk_size: usize,
    pub tolerance: Tolerance,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct OutputConfig {
    pub directory: PathBuf,
}

impl Config {
    /// Creates a default configuration with sensible defaults:
    /// - MS tolerance: 15 ppm
    /// - Mobility tolerance: 5%
    /// - Quad tolerance: 0.1 absolute
    /// - RT tolerance: unrestricted
    /// - Chunk size: 20000
    pub fn default_config() -> Self {
        Config {
            input: None,
            analysis: AnalysisConfig {
                dotd_files: None,
                chunk_size: 20000,
                tolerance: Tolerance {
                    ms: MzTolerance::Ppm((15.0, 15.0)),
                    mobility: MobilityTolerance::Pct((5.0, 5.0)),
                    quad: QuadTolerance::Absolute((0.1, 0.1)),
                    rt: RtTolerance::Unrestricted,
                },
            },
            output: None,
        }
    }
}
