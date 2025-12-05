//! Query execution implementations for indexed timsTOF data.
//!
//! This module implements the [`QueriableData`] trait for [`IndexedTimstofPeaks`],
//! defining how different aggregators extract and accumulate peak data.
//!
//! # Query Flow
//!
//! The query pattern follows these steps:
//!
//! 1. **Tolerance Application**: The [`Tolerance`] struct defines search windows
//!    for m/z, retention time, ion mobility, and quadrupole isolation
//! 2. **Range Calculation**: For each dimension, calculate the query range based
//!    on the target value and tolerance (e.g., `mz ± ppm` or `rt ± minutes`)
//! 3. **Peak Filtering**: Query the indexed data structure for peaks within all ranges
//! 4. **Aggregation**: Accumulate filtered peaks into the aggregator (sum intensities,
//!    build profiles, compute statistics, etc.)
//!
//! # Aggregator-Specific Behavior
//!
//! ## PointIntensityAggregator
//!
//! Sums total intensity across all matching peaks (both MS1 and MS2).
//!
//! - **Precursors**: Queries MS1 data without quadrupole filtering
//! - **Fragments**: Queries MS2 data with quadrupole isolation constraints
//! - **Result**: Single f64 intensity sum
//!
//! ## ChromatogramCollector
//!
//! Builds retention time profiles by binning intensities across RT dimension.
//!
//! - **Binning**: Maps each peak to RT bin based on reference RT array
//! - **Precursors**: MS1 chromatogram trace
//! - **Fragments**: MS2 chromatogram trace per fragment ion
//! - **Result**: Vector of intensities indexed by RT bin
//!
//! ## SpectralCollector
//!
//! Builds m/z spectra by accumulating intensities across mass dimension.
//!
//! - **Binning**: Maps each peak to m/z bin based on reference m/z array
//! - **Precursors**: MS1 spectrum
//! - **Fragments**: MS2 spectrum per fragment ion
//! - **Result**: 2D array (precursor/fragment × m/z bins)
//!
//! # Quadrupole Filtering
//!
//! For MS2 data (fragments), peaks are filtered by quadrupole isolation window:
//!
//! 1. Get quadrupole isolation range from tolerance
//! 2. Filter precursor frames that overlap this range
//! 3. For overlapping frames, constrain ion mobility range to actual isolation window
//! 4. Query peaks within constrained ranges
//!
//! This ensures fragment queries only include peaks from appropriate precursor isolations.
//!
//! # OptionallyRestricted Ranges
//!
//! Some tolerance dimensions use `OptionallyRestricted<TupleRange>` (from `timscentroid::utils`):
//!
//! - `Restricted(range)`: Normal bounded range (e.g., RT window)
//! - `Unrestricted`: No filtering on this dimension (match all values)
//!
//! This allows queries like "all retention times" while still filtering by m/z and mobility.

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
        let im_range = tolerance.mobility_range_f16(aggregator.eg.mobility_ook0());
        let rt_range_milliseconds = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds());
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
        aggregator.eg.iter_precursors().for_each(|(_idx, mz)| {
            let mz_range = tolerance.mz_range_f32(mz as f32);
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
                aggregator.eg.iter_fragments_refs().for_each(|(_idx, mz)| {
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
        let im_range: timscentroid::utils::OptionallyRestricted<TupleRange<f16>> =
            tolerance.mobility_range_f16(aggregator.eg.mobility_ook0());
        let rt_range_milliseconds = match tolerance
            .rt_range_as_milis(aggregator.eg.rt_seconds())
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
        let im_range = tolerance.mobility_range_f16(aggregator.eg.mobility_ook0());
        let rt_range_milliseconds = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds());
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
                let mz_range = tolerance.mz_range_f32(mz as f32);
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
        let im_range = tolerance.mobility_range_f16(aggregator.eg.mobility_ook0());
        let rt_range_milliseconds = tolerance.rt_range_as_milis(aggregator.eg.rt_seconds());
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
                let mz_range = tolerance.mz_range_f32(mz as f32);
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
