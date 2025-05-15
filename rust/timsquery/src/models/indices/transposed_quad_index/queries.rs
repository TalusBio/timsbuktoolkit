use super::QuadSplittedTransposedIndex;
use crate::OptionallyRestricted::{
    Restricted,
    Unrestricted,
};
use crate::models::aggregators::{
    ChromatogramCollector,
    PointIntensityAggregator,
    SpectralCollector,
};
use crate::models::frames::peak_in_quad::PeakInQuad;
use crate::traits::QueriableData;
use crate::{
    KeyLike,
    PeakAddable,
    Tolerance,
};

impl<FH: KeyLike> QueriableData<PointIntensityAggregator<FH>> for QuadSplittedTransposedIndex {
    fn add_query(&self, aggregator: &mut PointIntensityAggregator<FH>, tolerance: &Tolerance) {
        let quad_range = tolerance.quad_range(aggregator.eg.get_precursor_mz_limits());
        let scan_range =
            tolerance.indexed_scan_range(aggregator.eg.mobility as f64, &self.im_converter);
        let rt_range_ms = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds);
        aggregator.eg.precursors.iter().for_each(|(_idx, mz)| {
            let mz_range = tolerance.indexed_tof_range(*mz, &self.mz_converter);
            let mut clsr = |peak: PeakInQuad| {
                aggregator.intensity += peak.corrected_intensity as f64;
            };
            self.query_ms1_peaks(mz_range, scan_range, rt_range_ms, &mut clsr);
        });

        self.iter_msn_ranges()
            .filter(|(_qi, quad)| {
                quad.quad_settings
                    .expect("all ms2 indices should have quad info")
                    .matches_quad_settings(quad_range, scan_range)
            })
            .for_each(|(_ii, quad_idx)| {
                aggregator.eg.fragments.iter().for_each(|(_idx, mz)| {
                    // there is a bit of repeated work here ...
                    // calculating the mz range for each quad instead of once per ion ...
                    // todo: fix later ...
                    let mz_range = tolerance.indexed_tof_range(*mz, &self.mz_converter);
                    quad_idx
                        .query_peaks(mz_range, scan_range, rt_range_ms)
                        .for_each(|x| {
                            aggregator.intensity += x.corrected_intensity as f64;
                        });
                });
            });
    }
}

impl<FH: KeyLike> QueriableData<ChromatogramCollector<FH, f32>> for QuadSplittedTransposedIndex {
    fn add_query(&self, aggregator: &mut ChromatogramCollector<FH, f32>, tolerance: &Tolerance) {
        let quad_range = tolerance.quad_range(aggregator.eg.get_precursor_mz_limits());
        let scan_range =
            tolerance.indexed_scan_range(aggregator.eg.mobility as f64, &self.im_converter);
        let rt_range_ms = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds);
        // TODO: This section is kind of ugly ... I can probably move some of the logic around ...
        let query_rt_ms = aggregator.rt_range();
        let rt_range_ms = match tolerance.rt_range_as_milis(aggregator.eg.rt_seconds) {
            Unrestricted => rt_range_ms,
            Restricted(tol_range_rt) => {
                let tmp = query_rt_ms.intersection(tol_range_rt);
                if tmp.is_none() {
                    return;
                }
                Restricted(tmp.unwrap())
            }
        };

        aggregator
            .iter_mut_precursors()
            .for_each(|((_idx, mz), mut chr)| {
                let mz_range = tolerance.indexed_tof_range(*mz, &self.mz_converter);
                let mut clsr = |peak: PeakInQuad| {
                    chr.add_at_close_rt(peak.retention_time_ms, peak.corrected_intensity)
                };
                self.query_ms1_peaks(mz_range, scan_range, rt_range_ms, &mut clsr);
            });

        self.iter_msn_ranges()
            .filter(|(_qi, quad)| {
                quad.quad_settings
                    .expect("all ms2 indices should have quad info")
                    .matches_quad_settings(quad_range, scan_range)
            })
            .for_each(|(_ii, quad_idx)| {
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), mut chr)| {
                        // there is a bit of repeated work here ...
                        // calculating the mz range for each quad instead of once per ion ...
                        // todo: fix later ...
                        let mz_range = tolerance.indexed_tof_range(*mz, &self.mz_converter);
                        quad_idx
                            .query_peaks(mz_range, scan_range, rt_range_ms)
                            .for_each(|x| {
                                chr.add_at_close_rt(x.retention_time_ms, x.corrected_intensity);
                            });
                    });
            });
    }
}

impl<FH: KeyLike> QueriableData<SpectralCollector<FH, f32>> for QuadSplittedTransposedIndex {
    fn add_query(&self, aggregator: &mut SpectralCollector<FH, f32>, tolerance: &Tolerance) {
        let quad_range = tolerance.quad_range(aggregator.eg.get_precursor_mz_limits());
        let scan_range =
            tolerance.indexed_scan_range(aggregator.eg.mobility as f64, &self.im_converter);
        let rt_range_ms = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds);
        aggregator
            .iter_mut_precursors()
            .for_each(|((_idx, mz), ion)| {
                let mz_range = tolerance.indexed_tof_range(*mz, &self.mz_converter);
                let mut clsr = |peak: PeakInQuad| {
                    *ion += peak.corrected_intensity;
                };
                self.query_ms1_peaks(mz_range, scan_range, rt_range_ms, &mut clsr);
            });

        self.iter_msn_ranges()
            .filter(|(_qi, quad)| {
                quad.quad_settings
                    .expect("All MS2 indices should have quad info")
                    .matches_quad_settings(quad_range, scan_range)
            })
            .for_each(|(_ii, quad_idx)| {
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), ion)| {
                        // There is a bit of repeated work here ...
                        // calculating the mz range for each quad instead of once per ion ...
                        // TODO: Fix later ...
                        let mz_range = tolerance.indexed_tof_range(*mz, &self.mz_converter);
                        quad_idx
                            .query_peaks(mz_range, scan_range, rt_range_ms)
                            .for_each(|x| {
                                *ion += x.corrected_intensity;
                            });
                    });
            });
    }
}

// This is some ugly copy-pasted code but I will fix it later if I have to ...
impl<FH: KeyLike, V: PeakAddable> QueriableData<SpectralCollector<FH, V>>
    for QuadSplittedTransposedIndex
{
    fn add_query(&self, aggregator: &mut SpectralCollector<FH, V>, tolerance: &Tolerance) {
        let quad_range = tolerance.quad_range(aggregator.eg.get_precursor_mz_limits());
        let scan_range =
            tolerance.indexed_scan_range(aggregator.eg.mobility as f64, &self.im_converter);
        let rt_range_ms = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds);
        aggregator
            .iter_mut_precursors()
            .for_each(|((_idx, mz), ion)| {
                let mz_range = tolerance.indexed_tof_range(*mz, &self.mz_converter);
                let mut clsr = |peak: PeakInQuad| {
                    let tmp = peak.resolve(&self.im_converter, &self.mz_converter);
                    *ion += tmp;
                };
                self.query_ms1_peaks(mz_range, scan_range, rt_range_ms, &mut clsr);
            });

        self.iter_msn_ranges()
            .filter(|(_qi, quad)| {
                quad.quad_settings
                    .expect("All MS2 indices should have quad info")
                    .matches_quad_settings(quad_range, scan_range)
            })
            .for_each(|(_ii, quad_idx)| {
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), ion)| {
                        // There is a bit of repeated work here ...
                        // calculating the mz range for each quad instead of once per ion ...
                        // TODO: Fix later ...
                        let mz_range = tolerance.indexed_tof_range(*mz, &self.mz_converter);
                        quad_idx
                            .query_peaks(mz_range, scan_range, rt_range_ms)
                            .for_each(|peak| {
                                let tmp = peak.resolve(&self.im_converter, &self.mz_converter);
                                *ion += tmp;
                            });
                    });
            });
    }
}
