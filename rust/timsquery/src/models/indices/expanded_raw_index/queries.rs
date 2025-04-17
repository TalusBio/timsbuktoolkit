
use crate::{KeyLike, Tolerance};
use crate::traits::QueriableData;
use crate::models::aggregators::{EGCAggregator, EGSAggregator, RawPeakIntensityAggregator};
use super::model::ExpandedRawFrameIndex;

impl<FH: KeyLike> QueriableData<EGCAggregator<FH>>
    for ExpandedRawFrameIndex
{
    fn add_query(&self, aggregator: &mut EGCAggregator<FH>, tolerance: &Tolerance) {
        let quad_range = tolerance.quad_range(aggregator.eg.get_precursor_mz_limits());
        todo!();
    }
    // fn query(&self, fragment_query: &FragmentGroupIndexQuery<FH>) -> Vec<RawPeak> {
    //     let precursor_mz_range = (
    //         fragment_query.precursor_query.isolation_mz_range.0 as f64,
    //         fragment_query.precursor_query.isolation_mz_range.0 as f64,
    //     )
    //         .into();
    //     let scan_range = Some(fragment_query.precursor_query.mobility_index_range);

    //     fragment_query
    //         .mz_index_ranges
    //         .iter()
    //         .flat_map(|(_fh, tof_range)| {
    //             let mut local_vec: Vec<RawPeak> = vec![];
    //             self.query_peaks(
    //                 *tof_range,
    //                 precursor_mz_range,
    //                 scan_range,
    //                 fragment_query.precursor_query.frame_index_range,
    //                 &mut |x| local_vec.push(RawPeak::from(x)),
    //             );

    //             local_vec
    //         })
    //         .collect()
    // }

    // fn add_query<A, O, AG, C2>(
    //     &self,
    //     fragment_query: &FragmentGroupIndexQuery<FH>,
    //     aggregator: &mut AG,
    // ) where
    //     A: From<RawPeak> + Send + Sync + Clone + Copy,
    //     AG: Aggregator<Item = A, Output = O, Context = C2>,
    //     MsLevelContext<usize, FH>: Into<C2>,
    // {
    //     let precursor_mz_range = (
    //         fragment_query.precursor_query.isolation_mz_range.0 as f64,
    //         fragment_query.precursor_query.isolation_mz_range.0 as f64,
    //     )
    //         .into();
    //     let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
    //     fragment_query
    //         .mz_index_ranges
    //         .iter()
    //         .for_each(|(_fh, tof_range)| {
    //             self.query_peaks(
    //                 *tof_range,
    //                 precursor_mz_range,
    //                 scan_range,
    //                 fragment_query.precursor_query.frame_index_range,
    //                 &mut |peak| aggregator.add(RawPeak::from(peak)),
    //             );
    //         })
    // }

    // fn add_query_multi_group<A, O, AG, C2>(
    //     &self,
    //     fragment_queries: &[FragmentGroupIndexQuery<FH>],
    //     aggregator: &mut [AG],
    // ) where
    //     A: From<RawPeak> + Send + Sync + Clone + Copy,
    //     AG: Aggregator<Item = A, Output = O, Context = C2>,
    //     MsLevelContext<usize, FH>: Into<C2>,
    // {
    //     let prec_mz_ranges = fragment_queries
    //         .iter()
    //         .map(|x| {
    //             (
    //                 x.precursor_query.isolation_mz_range.start() as f64,
    //                 x.precursor_query.isolation_mz_range.end() as f64,
    //             )
    //                 .into()
    //         })
    //         .collect::<Vec<_>>();

    //     let scan_ranges = fragment_queries
    //         .iter()
    //         .map(|x| Some(x.precursor_query.mobility_index_range))
    //         .collect::<Vec<_>>();

    //     let frame_index_ranges = fragment_queries
    //         .iter()
    //         .map(|x| x.precursor_query.frame_index_range)
    //         .collect::<Vec<_>>();

    //     // Query the ms1 mzs first.
    //     aggregator.par_iter_mut().enumerate().for_each(|(i, agg)| {
    //         fragment_queries[i]
    //             .iter_ms1_mzs()
    //             .for_each(|(fh, mz_range)| {
    //                 if agg.supports_context() {
    //                     agg.set_context(fh.into());
    //                 }
    //                 self.query_ms1_peaks(
    //                     mz_range,
    //                     scan_ranges[i],
    //                     frame_index_ranges[i],
    //                     &mut |x| agg.add(RawPeak::from(x)),
    //                 );
    //             });
    //     });

    //     for quad_setting in self.flat_quad_settings.iter() {
    //         let local_index = quad_setting.index;

    //         let tqi = self
    //             .bundled_frames
    //             .get(&local_index)
    //             .expect("Only existing quads should be queried.");

    //         aggregator.par_iter_mut().enumerate().for_each(|(i, agg)| {
    //             if !matches_quad_settings(quad_setting, prec_mz_ranges[i], scan_ranges[i]) {
    //                 return;
    //             }

    //             for (fh, tof_range) in fragment_queries[i].iter_ms2_mzs() {
    //                 if agg.supports_context() {
    //                     agg.set_context(fh.into());
    //                 }
    //                 tqi.query_peaks(tof_range, scan_ranges[i], frame_index_ranges[i], &mut |x| {
    //                     agg.add(RawPeak::from(x))
    //                 });
    //             }
    //         });
    //     }
    // }
}


impl<FH: KeyLike> QueriableData<RawPeakIntensityAggregator<FH>>
    for ExpandedRawFrameIndex
{
    fn add_query(&self, aggregator: &mut RawPeakIntensityAggregator<FH>, tolerance: &Tolerance) {
        todo!();
    }
}

impl<FH: KeyLike> QueriableData<EGSAggregator<FH>>
    for ExpandedRawFrameIndex
{
    fn add_query(&self, aggregator: &mut EGSAggregator<FH>, tolerance: &Tolerance) {
        let quad_range = tolerance.quad_range(aggregator.eg.get_precursor_mz_limits());
        todo!();
    }
}
