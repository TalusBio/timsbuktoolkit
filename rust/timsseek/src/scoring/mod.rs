mod accumulator;
pub mod apex_finding;
pub mod full_results;
mod offsets;
pub mod pipeline;
mod scores;
pub mod search_results;
pub mod timings;

// RN I am not the biggest fan of exposig this
pub use scores::{
    coelution,
    hyperscore,
};

pub use pipeline::{
    ScoringPipeline,
    ToleranceHierarchy,
};
pub use search_results::IonSearchResults;
pub use timings::ScoreTimings;

pub const NUM_MS2_IONS: usize = 7;
pub const NUM_MS1_IONS: usize = 3;
pub const COELUTION_WINDOW_WIDTH: usize = 7;
