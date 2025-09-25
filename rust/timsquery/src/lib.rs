// Re-export main structures
pub use crate::models::Tolerance;
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

// Re-export traits
pub use crate::models::PeakAddable;
pub use crate::traits::queriable_data::{
    GenerallyQueriable,
    QueriableData,
};
pub use crate::traits::{
    KeyLike,
    ValueLike,
};
pub use timscentroid::utils::OptionallyRestricted;
pub use timscentroid::{
    CentroidingConfig,
    IndexedTimstofPeaks,
    TimsTofPath,
};

// Declare modules
pub mod errors;
pub mod models;
pub mod traits;
pub mod utils;
pub use crate::utils::TupleRange;

// Re-export errors
pub use crate::errors::{
    DataProcessingError,
    DataReadingError,
    TimsqueryError,
};
