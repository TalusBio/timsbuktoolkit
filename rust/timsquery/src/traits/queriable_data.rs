use rayon::prelude::*;
use crate::Tolerance;


/// Trait meant to signal that some form of data is
/// queriable using a specific aggregator.
///
/// Since each aggregator can optimized in very specific ways
/// and mean very different things it makes sense to have
/// an implementation for each aggregator + index combination.
///
/// Note that its necessary to implement either `add_query` or
/// `add_query_multi_group` (or both), since the default implementation
/// the default implementations rely on each other. (and not implementing
/// any will result in infinite recursion)
///
pub trait QueriableData<QA>
where
    QA: Send + Sync,
    Self: Send + Sync,
{
    fn add_query(&self, queriable_aggregator: &mut QA, tolerance: &Tolerance) {
        self.add_query_multi_group(&mut [queriable_aggregator], tolerance);
    }

    fn par_add_query_multi(&self, queriable_aggregators: &mut [QA], tolerance: &Tolerance)
    {
        queriable_aggregators
            .par_iter_mut()
            .map(|queriable_aggregator| self.add_query(queriable_aggregator, tolerance));
    }
}
