use super::key_like::ValueLike;
use crate::models::aggregators::{
    ChromatogramCollector,
    MzMobilityStatsCollector,
    PointIntensityAggregator,
    SpectralCollector,
};
use crate::models::elution_group::TimsElutionGroup;
use crate::{
    KeyLike,
    Tolerance,
};
use rayon::prelude::*;
use std::ops::AddAssign;
use timscentroid::indexing::IndexedPeak;
use timscentroid::rt_mapping::RTIndex;

/// Data needed to build a query against the index.
///
/// Exposes exactly what `QueryRanges::from_query_data` and per-aggregator
/// `add_query` impls read — scalars for range building + iterators for
/// precursor and fragment m/z windows. A fully immutable view, implemented
/// by `TimsElutionGroup` as well as by aggregators that carry their own
/// query scalars.
pub trait HasQueryData<FH: KeyLike> {
    fn id(&self) -> u64;
    fn precursor_mz_limits(&self) -> (f64, f64);
    fn mobility_ook0(&self) -> f32;
    fn rt_seconds(&self) -> f32;
    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> + '_;
    /// Borrowed label + copied mz. Named `iter_fragments` (not `_refs`)
    /// since the mz is Copy — the `_refs` name would mislead readers.
    fn iter_fragments<'a>(&'a self) -> impl Iterator<Item = (&'a FH, f64)> + 'a
    where
        FH: 'a;
}

impl<FH: KeyLike> HasQueryData<FH> for TimsElutionGroup<FH> {
    fn id(&self) -> u64 {
        TimsElutionGroup::id(self)
    }

    fn precursor_mz_limits(&self) -> (f64, f64) {
        self.get_precursor_mz_limits()
    }

    fn mobility_ook0(&self) -> f32 {
        TimsElutionGroup::mobility_ook0(self)
    }

    fn rt_seconds(&self) -> f32 {
        TimsElutionGroup::rt_seconds(self)
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> + '_ {
        TimsElutionGroup::iter_precursors(self)
    }

    fn iter_fragments<'a>(&'a self) -> impl Iterator<Item = (&'a FH, f64)> + 'a
    where
        FH: 'a,
    {
        self.iter_fragments_refs().map(|(k, mz)| (k, *mz))
    }
}

/// Trait indicating that indexed data can be queried with a specific aggregator type.
///
/// This trait connects indexed peak data (typically `IndexedTimstofPeaks`) with
/// aggregators that accumulate query results. Each aggregator type has its own
/// implementation with optimized query patterns.
///
/// # Type Parameters
///
/// - `QA`: The aggregator type (e.g., [`ChromatogramCollector`], [`SpectralCollector`])
///
/// # Query Pattern
///
/// The typical usage flow:
///
/// 1. Create an aggregator with an `ElutionGroup` defining target m/z values
/// 2. Create a `Tolerance` defining search windows (m/z, RT, mobility, quad)
/// 3. Call `add_query()` to extract and accumulate matching peaks
/// 4. Repeat for multiple queries or use `par_add_query_multi()` for batch queries
///
/// # Example
///
/// ```ignore
/// use timsquery::{IndexedTimstofPeaks, ChromatogramCollector, ElutionGroup, Tolerance, QueriableData};
/// use std::sync::Arc;
///
/// let peaks: IndexedTimstofPeaks = unimplemented!();
/// let elution_group = Arc::new(ElutionGroup::default());
/// let ref_rt = Arc::new(vec![0u32; 100]);
/// let mut aggregator = ChromatogramCollector::new(elution_group, ref_rt)?;
/// let tolerance = Tolerance::default();
///
/// // Query peaks and accumulate into chromatogram
/// peaks.add_query(&mut aggregator, &tolerance);
///
/// // The aggregator is modified in place with results
/// ```
///
/// # Implementation Notes
///
/// Each aggregator has a tailored implementation because:
/// - Different aggregators need different binning strategies
/// - Query patterns can be optimized per aggregator type
/// - Some aggregators need extra metadata (e.g., reference RT/m/z arrays)
///
/// See [`crate::models::indexed_data`] for implementation details.
pub trait QueriableData<QA>
where
    QA: Send + Sync,
    Self: Send + Sync,
{
    /// Execute a query and accumulate results into the aggregator.
    ///
    /// This method:
    /// 1. Applies tolerance ranges to filter peaks
    /// 2. Iterates over matching peaks
    /// 3. Accumulates peak data into the aggregator
    ///
    /// # Arguments
    ///
    /// - `queriable_aggregator`: Mutable aggregator to accumulate results into
    /// - `tolerance`: Defines search windows for m/z, RT, mobility, and quadrupole
    fn add_query(&self, queriable_aggregator: &mut QA, tolerance: &Tolerance);

    /// Execute multiple queries in parallel, one per aggregator.
    ///
    /// Zips aggregators with tolerances and processes each pair in parallel via
    /// Rayon. Both parameters accept any `IntoParallelIterator`, so you can pass:
    ///
    /// - `&mut [QA]` for aggregators (the common case)
    /// - `&[Tolerance]` or `&Vec<Tolerance>` for per-query tolerances
    /// - `rayon::iter::repeatn(&tol, n)` for a single shared tolerance
    ///
    /// # Arguments
    ///
    /// - `queriable_aggregators`: Parallel iterator of mutable aggregator references
    /// - `tolerances`: Parallel iterator of tolerance references (one per aggregator)
    fn par_add_query_multi<'a, A, T>(&self, queriable_aggregators: A, tolerances: T)
    where
        QA: 'a,
        A: IntoParallelIterator<Item = &'a mut QA>,
        A::Iter: IndexedParallelIterator,
        T: IntoParallelIterator<Item = &'a Tolerance>,
        T::Iter: IndexedParallelIterator,
    {
        queriable_aggregators
            .into_par_iter()
            .zip(tolerances)
            .for_each(|(agg, tol)| self.add_query(agg, tol));
    }
}

/// Marker trait for value types that can accumulate indexed peaks.
///
/// This trait is used by [`SpectralCollector`] to constrain value types that
/// can directly accumulate peak data. Types must:
///
/// - Implement `AddAssign<IndexedPeak>` for peak accumulation
/// - Implement [`ValueLike`] (copyable, serializable, thread-safe)
/// - Implement `Default` for initialization
///
/// Common implementations: `f32`, custom statistics collectors.
pub trait PeakAddable<T: RTIndex>: AddAssign<IndexedPeak<T>> + ValueLike + Default {}

/// Convenience trait for indexed data that supports all common aggregator types.
///
/// This trait indicates that an index (like `IndexedTimstofPeaks`) can be queried
/// with the standard set of aggregators. It's a shorthand for requiring multiple
/// [`QueriableData`] implementations.
///
/// # Standard Aggregators
///
/// This trait requires implementations for:
/// - `SpectralCollector<T, f32>` - m/z spectra with float intensities
/// - `SpectralCollector<T, MzMobilityStatsCollector>` - spectra with statistics
/// - `ChromatogramCollector<T, f32>` - retention time profiles
/// - `PointIntensityAggregator<T>` - total intensity sums
///
/// # Why These Four?
///
/// These are ALL the aggregators currently exposed by timsquery:
/// - **Spectra**: Analyze mass distributions (with float or stats collectors)
/// - **Chromatograms**: Analyze elution profiles
/// - **Point intensities**: Simple intensity extraction
///
/// This trait guarantees support for the complete set of timsquery aggregators.
///
/// # Design Note
///
/// This trait uses a blanket implementation - any type implementing the four
/// required `QueriableData` traits automatically implements `GenerallyQueriable`.
///
/// **Important**: When adding new aggregators to timsquery, this trait MUST be updated
/// to include them. This is intentional - `GenerallyQueriable` represents "full support"
/// for the public timsquery API.
///
/// For domain-specific trait bounds that only require a subset of aggregators,
/// define custom traits (e.g., `timsseek::ScorerQueriable`).
///
/// # Example
///
/// ```ignore
/// use timsquery::{
///     IndexedTimstofPeaks, GenerallyQueriable, ChromatogramCollector,
///     SpectralCollector, PointIntensityAggregator, QueriableData
/// };
///
/// // IndexedTimstofPeaks implements GenerallyQueriable<usize>
/// fn process_data<T: GenerallyQueriable<usize>>(data: &T) {
///     // Can use any of the four standard aggregators:
///     // - ChromatogramCollector<usize, f32>
///     // - SpectralCollector<usize, f32>
///     // - SpectralCollector<usize, MzMobilityStatsCollector>
///     // - PointIntensityAggregator<usize>
/// }
///
/// # let peaks: IndexedTimstofPeaks = unimplemented!();
/// process_data(&peaks);  // ✓ compiles
/// ```
pub trait GenerallyQueriable<T: KeyLike>:
    QueriableData<SpectralCollector<T, f32>>
    + QueriableData<SpectralCollector<T, MzMobilityStatsCollector>>
    + QueriableData<ChromatogramCollector<T, f32>>
    + QueriableData<PointIntensityAggregator<T>>
{
}
impl<
    T: KeyLike,
    I: QueriableData<SpectralCollector<T, f32>>
        + QueriableData<ChromatogramCollector<T, f32>>
        + QueriableData<PointIntensityAggregator<T>>
        + QueriableData<SpectralCollector<T, MzMobilityStatsCollector>>,
> GenerallyQueriable<T> for I
{
}
