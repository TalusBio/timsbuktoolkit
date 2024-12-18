// Re-export main structures
pub use crate::models::elution_group::ElutionGroup;
pub use crate::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;

// Re-export traits
pub use crate::traits::aggregator::Aggregator;
pub use crate::traits::queriable_data::QueriableData;
pub use crate::traits::tolerance::{
    Tolerance,
    ToleranceAdapter,
};

// Declare modules
pub mod errors;
pub mod models;
pub mod queriable_tims_data;
pub mod traits;
pub mod utils;

// Re-export errors
pub use crate::errors::{
    DataProcessingError,
    DataReadingError,
    TimsqueryError,
};
