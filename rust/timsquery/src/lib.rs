#![doc = include_str!("../README.md")]

// Re export stuff from other crates ...
pub use {
    micromzpaf,
    tinyvec,
};

// Re-export main structures
pub use crate::models::Tolerance;
pub use crate::models::aggregators::{
    ChromatogramCollector,
    MzMobilityStatsCollector,
    PointIntensityAggregator,
    SpectralCollector,
};
pub use crate::models::base::{
    Array2D,
    MzMajorIntensityArray,
    RTMajorIntensityArray,
};
pub use crate::models::elution_group::TimsElutionGroup;

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
pub mod serde;
pub mod traits;
pub mod utils;
pub use crate::utils::TupleRange;

// Re-export errors
pub use crate::errors::{
    DataProcessingError,
    DataReadingError,
    TimsqueryError,
};

pub mod ion {
    pub use micromzpaf::{
        IonAnnot,
        IonParsingError,
        IonSeriesTerminality,
    };
}
