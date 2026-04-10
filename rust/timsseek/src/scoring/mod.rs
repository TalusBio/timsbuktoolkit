mod accumulator;
pub mod apex_finding;
pub mod full_results;
mod offsets;
pub mod parquet_writer;
pub mod pipeline;
pub mod results;
mod scores;
pub mod search_results;
pub mod timings;

pub use scores::hyperscore;

pub use pipeline::{
    CalibrantCandidate,
    CalibrantHeap,
    CalibrationConfig,
    Scorer,
};
pub use results::{ScoredCandidate, CompetedCandidate, FinalResult, ScoringFields};
pub use search_results::IonSearchResults;
pub use timings::{PipelineTimings, ScoreTimings};

pub const NUM_MS2_IONS: usize = 7;
pub const NUM_MS1_IONS: usize = 3;
