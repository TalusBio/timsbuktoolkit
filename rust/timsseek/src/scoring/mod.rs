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
