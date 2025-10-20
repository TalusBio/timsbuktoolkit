use crate::models::aggregators::{
    ChromatogramCollector,
    PointIntensityAggregator,
    SpectralCollector,
};
use crate::traits::{
    PeakAddable,
    QueriableData,
};
use crate::{
    KeyLike,
    Tolerance,
};
use half::f16;
use timscentroid::IndexedTimstofPeaks;
use timscentroid::utils::OptionallyRestricted::{
    Restricted,
    Unrestricted,
};
use timscentroid::utils::TupleRange;

impl<FH: KeyLike> QueriableData<PointIntensityAggregator<FH>> for IndexedTimstofPeaks {
    fn add_query(&self, aggregator: &mut PointIntensityAggregator<FH>, tolerance: &Tolerance) {
        let prec_mz_limits: TupleRange<f64> =
            aggregator.eg.get_precursor_mz_limits().try_into().unwrap();
        let quad_range =
            tolerance.quad_range_f32((prec_mz_limits.start() as f32, prec_mz_limits.end() as f32));
        let im_range = tolerance.mobility_range_f16(aggregator.eg.mobility);
        let rt_range_milliseconds = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds);
        let cycle_range = match rt_range_milliseconds {
            Restricted(x) => Restricted(
                TupleRange::try_new(
                    self.rt_ms_to_cycle_index(x.start()),
                    self.rt_ms_to_cycle_index(x.end()),
                )
                .unwrap(),
            ),
            Unrestricted => Unrestricted,
        };
        aggregator.eg.precursors.iter().for_each(|(_idx, mz)| {
            let mz_range = tolerance.mz_range_f32(*mz as f32);
            self.query_peaks_ms1(mz_range, cycle_range, im_range)
                .for_each(|peak| {
                    aggregator.intensity += peak.intensity as f64;
                });
        });

        self.filter_precursor_ranges(quad_range, im_range)
            .for_each(|(quad_info, peaks)| {
                let constrained_im = im_range.map(|x| {
                    let tmp = quad_info.intersects_ranges(
                        (quad_range.start() as f64, quad_range.end() as f64),
                        (x.start().to_f64(), x.end().to_f64()),
                    );
                    let tmp = tmp.expect(
                        "Since we filtered based on the mz+ims there should always be a match",
                    );
                    (f16::from_f32(tmp.0 as f32), f16::from_f32(tmp.1 as f32))
                        .try_into()
                        .unwrap()
                });
                aggregator.eg.fragments.iter().for_each(|(_idx, mz)| {
                    let mz_range = tolerance.mz_range_f32(*mz as f32);
                    peaks
                        .query_peaks(mz_range, cycle_range, constrained_im)
                        .for_each(|peak| {
                            aggregator.intensity += peak.intensity as f64;
                        });
                });
            })
    }
}

impl<FH: KeyLike> QueriableData<ChromatogramCollector<FH, f32>> for IndexedTimstofPeaks {
    fn add_query(&self, aggregator: &mut ChromatogramCollector<FH, f32>, tolerance: &Tolerance) {
        let agg_rt_limits_seconds = aggregator.rt_range_milis();
        let prec_mz_limits = aggregator.eg.get_precursor_mz_limits();
        let quad_range =
            tolerance.quad_range_f32((prec_mz_limits.0 as f32, prec_mz_limits.1 as f32));
        let im_range = tolerance.mobility_range_f16(aggregator.eg.mobility);
        let rt_range_milliseconds = match tolerance
            .rt_range_as_milis(aggregator.eg.rt_seconds)
            .map(|x| x.try_intercept(agg_rt_limits_seconds))
        {
            Restricted(Some(x)) => Restricted(x),
            Restricted(None) => return,
            Unrestricted => Unrestricted,
        };
        let cycle_range = match rt_range_milliseconds {
            Restricted(x) => Restricted(
                TupleRange::try_new(
                    self.rt_ms_to_cycle_index(x.start()),
                    self.rt_ms_to_cycle_index(x.end()),
                )
                .unwrap(),
            ),
            Unrestricted => Unrestricted,
        };
        aggregator
            .iter_mut_precursors()
            .for_each(|((_idx, mz), mut chr)| {
                let mz_range = tolerance.mz_range_f32(*mz as f32);
                self.query_peaks_ms1(mz_range, cycle_range, im_range)
                    .for_each(|peak| chr.add_at_index(peak.cycle_index, peak.intensity));
            });

        self.filter_precursor_ranges(quad_range, im_range)
            .for_each(|(quad_info, peaks)| {
                let constrained_im = im_range.map(|x| {
                    let tmp = quad_info.intersects_ranges(
                        (quad_range.start() as f64, quad_range.end() as f64),
                        (x.start().to_f64(), x.end().to_f64()),
                    );
                    let tmp = tmp.expect(
                        "Since we filtered based on the mz+ims there should always be a match",
                    );
                    (f16::from_f32(tmp.0 as f32), f16::from_f32(tmp.1 as f32))
                        .try_into()
                        .unwrap()
                });
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), mut chr)| {
                        // there is a bit of repeated work here ...
                        // calculating the mz range for each quad instead of once per ion ...
                        // TODO: fix later ...

                        let mz_range = tolerance.mz_range_f32(*mz as f32);
                        peaks
                            .query_peaks(mz_range, cycle_range, constrained_im)
                            .for_each(|x| {
                                chr.add_at_index(x.cycle_index, x.intensity);
                            });
                    });
            })
    }
}

impl<FH: KeyLike> QueriableData<SpectralCollector<FH, f32>> for IndexedTimstofPeaks {
    fn add_query(&self, aggregator: &mut SpectralCollector<FH, f32>, tolerance: &Tolerance) {
        let prec_mz_limits = aggregator.eg.get_precursor_mz_limits();
        let quad_range =
            tolerance.quad_range_f32((prec_mz_limits.0 as f32, prec_mz_limits.1 as f32));
        let im_range = tolerance.mobility_range_f16(aggregator.eg.mobility);
        let rt_range_milliseconds = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds);
        let cycle_range = match rt_range_milliseconds {
            Restricted(x) => Restricted(
                TupleRange::try_new(
                    self.rt_ms_to_cycle_index(x.start()),
                    self.rt_ms_to_cycle_index(x.end()),
                )
                .unwrap(),
            ),
            Unrestricted => Unrestricted,
        };
        aggregator
            .iter_mut_precursors()
            .for_each(|((_idx, mz), ion)| {
                let mz_range = tolerance.mz_range_f32(*mz as f32);
                self.query_peaks_ms1(mz_range, cycle_range, im_range)
                    .for_each(|peak| {
                        *ion += peak.intensity;
                    });
            });

        self.filter_precursor_ranges(quad_range, im_range)
            .for_each(|(quad_info, peaks)| {
                let constrained_im = im_range.map(|x| {
                    let tmp = quad_info.intersects_ranges(
                        (quad_range.start() as f64, quad_range.end() as f64),
                        (x.start().to_f64(), x.end().to_f64()),
                    );
                    let tmp = tmp.expect(
                        "Since we filtered based on the mz+ims there should always be a match",
                    );
                    (f16::from_f32(tmp.0 as f32), f16::from_f32(tmp.1 as f32))
                        .try_into()
                        .unwrap()
                });
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), ion)| {
                        // There is a bit of repeated work here ...
                        // calculating the mz range for each quad instead of once per ion ...
                        // TODO: Fix later ...
                        let mz_range = tolerance.mz_range_f32(*mz as f32);
                        peaks
                            .query_peaks(mz_range, cycle_range, constrained_im)
                            .for_each(|x| {
                                *ion += x.intensity;
                            });
                    });
            });
    }
}

// This is some ugly copy-pasted code but I will fix it later if I have to ...
impl<FH: KeyLike, V: PeakAddable> QueriableData<SpectralCollector<FH, V>> for IndexedTimstofPeaks {
    fn add_query(&self, aggregator: &mut SpectralCollector<FH, V>, tolerance: &Tolerance) {
        let prec_mz_limits = aggregator.eg.get_precursor_mz_limits();
        let quad_range =
            tolerance.quad_range_f32((prec_mz_limits.0 as f32, prec_mz_limits.1 as f32));
        let im_range = tolerance.mobility_range_f16(aggregator.eg.mobility);
        let rt_range_milliseconds = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds);
        let cycle_range = match rt_range_milliseconds {
            Restricted(x) => Restricted(
                TupleRange::try_new(
                    self.rt_ms_to_cycle_index(x.start()),
                    self.rt_ms_to_cycle_index(x.end()),
                )
                .unwrap(),
            ),
            Unrestricted => Unrestricted,
        };
        aggregator
            .iter_mut_precursors()
            .for_each(|((_idx, mz), ion)| {
                let mz_range = tolerance.mz_range_f32(*mz as f32);
                self.query_peaks_ms1(mz_range, cycle_range, im_range)
                    .for_each(|peak| {
                        *ion += *peak;
                    });
            });

        self.filter_precursor_ranges(quad_range, im_range)
            .for_each(|(quad_info, peaks)| {
                let constrained_im = im_range.map(|x| {
                    let tmp = quad_info.intersects_ranges(
                        (quad_range.start() as f64, quad_range.end() as f64),
                        (x.start().to_f64(), x.end().to_f64()),
                    );
                    let tmp = tmp.expect(
                        "Since we filtered based on the mz+ims there should always be a match",
                    );
                    (f16::from_f32(tmp.0 as f32), f16::from_f32(tmp.1 as f32))
                        .try_into()
                        .unwrap()
                });
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), ion)| {
                        // There is a bit of repeated work here ...
                        // calculating the mz range for each quad instead of once per ion ...
                        // TODO: Fix later ...
                        let mz_range = tolerance.mz_range_f32(*mz as f32);
                        peaks
                            .query_peaks(mz_range, cycle_range, constrained_im)
                            .for_each(|x| {
                                *ion += *x;
                            });
                    });
            });
    }
}
