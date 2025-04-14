pub mod multi_chromatogram_agg;
pub mod point_agg;

pub use multi_chromatogram_agg::MultiCMGStatsAgg;
// TODO: reorganize this so I donr use direcly from `base``
pub use multi_chromatogram_agg::MultiCMGStatsFactory;
pub use point_agg::{
    RawPeakIntensityAggregator,
    RawPeakVectorAggregator,
};
