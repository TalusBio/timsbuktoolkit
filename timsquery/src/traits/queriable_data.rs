use crate::traits::aggregator::{
    Aggregator,
    ProvidesContext,
};
use rayon::prelude::*;

pub trait QueriableData<QF, I, C>
where
    QF: Send + Sync + ProvidesContext<Context = C>,
    I: Send + Sync + Clone + Copy,
    Self: Send + Sync,
{
    fn query(&self, fragment_query: &QF) -> Vec<I>;
    fn add_query<A, O, AG, C2>(&self, fragment_query: &QF, aggregator: &mut AG)
    where
        A: From<I> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O, Context = C2>,
        C: Into<C2>;

    fn add_query_multi_group<A, O, AG, C2>(&self, fragment_queries: &[QF], aggregator: &mut [AG])
    where
        A: From<I> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O, Context = C2>,
        C: Into<C2>,
    {
        fragment_queries
            .par_iter()
            .zip(aggregator.par_iter_mut())
            .for_each(|(query, agg)| self.add_query(query, agg));
    }
}
