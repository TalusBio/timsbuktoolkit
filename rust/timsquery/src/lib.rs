// Re-export main structures
pub use crate::models::Tolerance;
pub use crate::models::aggregators::{
    EGCAggregator,
    EGSAggregator,
};
pub use crate::models::base::{
    Array2D,
    MzMajorIntensityArray,
    RTMajorIntensityArray,
};
pub use crate::models::elution_group::ElutionGroup;
pub use crate::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;

// Re-export traits
pub use crate::traits::KeyLike;
pub use crate::traits::queriable_data::{
    GenerallyQueriable,
    QueriableData,
};

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
