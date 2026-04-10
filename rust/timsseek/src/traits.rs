//! Trait definitions for timsseek query patterns.

use crate::IonAnnot;
use timsquery::{
    ChromatogramCollector,
    MzMobilityStatsCollector,
    QueriableData,
    SpectralCollector,
};

/// Trait for indexed data that supports the aggregators needed by [`crate::scoring::pipeline::Scorer`].
///
/// This trait is a convenience bound that documents exactly what query capabilities
/// the scoring engine requires. It's more specific than [`timsquery::GenerallyQueriable`]
/// which includes all standard aggregators.
///
/// # Required Aggregators
///
/// The scorer needs three specific aggregator types:
///
/// - **`ChromatogramCollector<IonAnnot, f32>`**: For prescoring and RT profile extraction
/// - **`SpectralCollector<IonAnnot, MzMobilityStatsCollector>`**: For secondary queries with statistics
/// - **`SpectralCollector<IonAnnot, f32>`**: For isotope pattern scoring
///
/// # Why This Trait?
///
/// Instead of using the general `GenerallyQueriable<IonAnnot>` trait (which includes
/// `PointIntensityAggregator` that Scorer doesn't use), this trait:
///
/// - Documents exactly what Scorer needs
/// - Makes function signatures more explicit
/// - Allows future pipeline variants to have different requirements
///
/// # Example
///
/// ```ignore
/// use timsseek::{Scorer, ScorerQueriable};
/// use timscentroid::IndexedTimstofPeaks;
/// use timsquery::Tolerance;
///
/// # let peaks: IndexedTimstofPeaks = unimplemented!();
/// # let broad_tol = Tolerance::default();
/// # let secondary_tol = Tolerance::default();
/// # let fragmented_range = (400.0, 1200.0).try_into().unwrap();
///
/// // IndexedTimstofPeaks implements ScorerQueriable
/// let scorer = Scorer {
///     index: peaks,
///     broad_tolerance: broad_tol,
///     secondary_tolerance: secondary_tol,
///     fragmented_range,
/// };
///
/// // Scorer can now use ChromatogramCollector, SpectralCollector variants
/// ```
///
/// # Implementation Note
///
/// This uses a blanket implementation - any type implementing the three required
/// `QueriableData` traits automatically implements `ScorerQueriable`.
pub trait ScorerQueriable:
    QueriableData<ChromatogramCollector<IonAnnot, f32>>
    + QueriableData<SpectralCollector<IonAnnot, MzMobilityStatsCollector>>
    + QueriableData<SpectralCollector<IonAnnot, f32>>
    + MappableRTCycles
{
}

// Blanket implementation: any type with the required QueriableData impls gets ScorerQueriable
impl<I> ScorerQueriable for I where
    I: QueriableData<ChromatogramCollector<IonAnnot, f32>>
        + QueriableData<SpectralCollector<IonAnnot, MzMobilityStatsCollector>>
        + QueriableData<SpectralCollector<IonAnnot, f32>>
        + MappableRTCycles
{
}

pub trait MappableRTCycles {
    fn ms1_cycle_mapping(
        &self,
    ) -> &timscentroid::rt_mapping::CycleToRTMapping<timscentroid::rt_mapping::MS1CycleIndex>;
}

impl MappableRTCycles for timscentroid::IndexedTimstofPeaks {
    fn ms1_cycle_mapping(
        &self,
    ) -> &timscentroid::rt_mapping::CycleToRTMapping<timscentroid::rt_mapping::MS1CycleIndex> {
        self.ms1_cycle_mapping()
    }
}
