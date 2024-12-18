pub mod aggregator;
pub mod arrays;
#[allow(clippy::module_inception)]
pub mod multi_chromatogram_agg;
pub mod scored_arrays;

pub use multi_chromatogram_agg::{
    MultiCMGStatsAgg,
    MultiCMGStatsFactory,
};
