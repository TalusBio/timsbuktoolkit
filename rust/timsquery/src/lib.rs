// Re-export main structures
pub use crate::models::aggregators::{
    ChromatogramCollector,
    MzMobilityStatsCollector,
    SpectralCollector,
};
pub use crate::models::base::{
    Array2D,
    MzMajorIntensityArray,
    RTMajorIntensityArray,
};
pub use crate::models::elution_group::ElutionGroup;
pub use crate::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
pub use crate::models::{
    OptionallyRestricted,
    Tolerance,
};

// Re-export traits
pub use crate::traits::queriable_data::{
    GenerallyQueriable,
    QueriableData,
};
pub use crate::traits::{
    KeyLike,
    PeakAddable,
    ValueLike,
};

// Declare modules
pub mod errors;
pub mod models;
pub mod traits;
pub mod utils;

pub use utils::tolerance_ranges::IncludedRange;

// Re-export errors
pub use crate::errors::{
    DataProcessingError,
    DataReadingError,
    TimsqueryError,
};
