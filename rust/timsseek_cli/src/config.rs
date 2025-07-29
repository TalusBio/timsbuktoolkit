use serde::{
    Deserialize,
    Serialize,
};
use timsseek::IonAnnot;
use std::path::PathBuf;
use timsquery::Tolerance;
use timsquery::models::indices::{
    ExpandedRawFrameIndex,
    QuadSplittedTransposedIndex,
};
use timsquery::GenerallyQueriable;
use std::sync::Arc;
use timsquery::IncludedRange;

use crate::cli::Cli;
use crate::errors;

#[derive(Debug, Serialize, Deserialize, Clone, Default)]
pub enum IndexType {
    #[default]
    Transposed,
    Expanded,
}

pub struct IndexElements {
    pub index: Box<dyn GenerallyQueriable<IonAnnot>>,
    pub index_cycle_rt_ms: Arc<[u32]>,
    pub fragmented_range: IncludedRange<f64>,
    
}

impl IndexType {
    pub fn build_index(
        &self,
        raw_file_path: &str,
    ) -> IndexElements {
        match self {
            IndexType::Expanded => {
                // Throughput seems to be ~ 30% better if I use the centroided version
                // But I like the idea of having the full resolution data available +
                // for the high throughput use case, we have the transposed index.
                let tmp = ExpandedRawFrameIndex::from_path(raw_file_path).unwrap();
                let rts = tmp.cycle_rt_ms.clone();
                let fragmented_range = tmp.fragmented_range();
                IndexElements {
                    index: Box::new(tmp),
                    index_cycle_rt_ms: rts,
                    fragmented_range,
                }
            }
            IndexType::Transposed => {
                let tmp = QuadSplittedTransposedIndex::from_path_centroided(raw_file_path).unwrap();
                let rts = tmp.cycle_rt_ms.clone();
                let fragmented_range = tmp.fragmented_range();
                IndexElements {
                    index: Box::new(tmp),
                    index_cycle_rt_ms: rts,
                    fragmented_range,
                }
            }
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Config {
    pub input: Option<InputConfig>,
    pub analysis: AnalysisConfig,
    pub output: Option<OutputConfig>,
    #[serde(default)]
    pub index_type: IndexType,
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

impl Config {
    #![allow(dead_code)]
    pub fn with_cli_args(_config: Cli) -> Result<Self, errors::CliError> {
        unimplemented!()
    }
}
