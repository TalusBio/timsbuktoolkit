pub mod calculate_scores;
pub mod full_results;
mod offsets;
pub mod scorer;
mod scores;
pub mod search_results;

// RN I am not the biggest fan of exposig this
pub use scores::{
    coelution,
    hyperscore,
};

pub use scorer::Scorer;
pub use search_results::IonSearchResults;

const NUM_MS2_IONS: usize = 7;
const NUM_MS1_IONS: usize = 3;
const COELUTION_WINDOW_WIDTH: usize = 7;
