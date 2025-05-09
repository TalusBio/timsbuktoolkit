use crate::models::aggregators::{
    ChromatogramCollector,
    PointIntensityAggregator,
    SpectralCollector,
};
use crate::{
    KeyLike,
    Tolerance,
};
use rayon::prelude::*;

/// Trait meant to signal that some form of data is
/// queriable using a specific aggregator.
///
/// Since each aggregator can optimized in very specific ways
/// and mean very different things it makes sense to have
/// an implementation for each aggregator + index combination.
///
/// Note that its necessary to implement either `add_query`
/// (or both that and `add_query_multi_group`)
pub trait QueriableData<QA>
where
    QA: Send + Sync,
    Self: Send + Sync,
{
    fn add_query(&self, queriable_aggregator: &mut QA, tolerance: &Tolerance);

    fn par_add_query_multi(&self, queriable_aggregators: &mut [QA], tolerance: &Tolerance) {
        queriable_aggregators
            .par_iter_mut()
            .for_each(|queriable_aggregator| self.add_query(queriable_aggregator, tolerance));
    }
}

// Blanket trait implementation meaning that the index can be queried with any aggregator.
pub trait GenerallyQueriable<T: KeyLike>:
    QueriableData<SpectralCollector<T, f32>>
    + QueriableData<ChromatogramCollector<T>>
    + QueriableData<PointIntensityAggregator<T>>
{
}
impl<
    T: KeyLike,
    I: QueriableData<SpectralCollector<T, f32>>
        + QueriableData<ChromatogramCollector<T>>
        + QueriableData<PointIntensityAggregator<T>>,
> GenerallyQueriable<T> for I
{
}
