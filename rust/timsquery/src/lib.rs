// Re-export main structures
pub use crate::models::elution_group::ElutionGroup;
pub use crate::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
pub use crate::models::Tolerance;

// Re-export traits
pub use crate::traits::queriable_data::{QueriableData, GenerallyQueriable};
pub use crate::traits::KeyLike;

// Declare modules
pub mod errors;
pub mod models;
pub mod traits;
pub mod utils;

// Re-export errors
pub use crate::errors::{
    DataProcessingError,
    DataReadingError,
    TimsqueryError,
};
