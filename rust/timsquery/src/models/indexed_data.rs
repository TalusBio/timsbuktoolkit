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
use crate::serde::IndexedPeaksHandle;
use crate::traits::{
    PeakAddable,
    QueriableData,
};
use crate::{
    KeyLike,
    OptionallyRestricted,
    Tolerance,
};
use half::f16;
use timscentroid::IndexedTimstofPeaks;
use timscentroid::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
    WindowCycleIndex,
};
use timscentroid::utils::OptionallyRestricted::{
    Restricted,
    Unrestricted,
};
use timscentroid::utils::TupleRange;

/// Encapsulates the query ranges computed from an elution group and tolerance.
///
/// This struct centralizes the tolerance-to-range conversion logic that's used
/// across all aggregator implementations, reducing code duplication.
struct QueryRanges {
    quad_range: TupleRange<f32>,
    im_range: OptionallyRestricted<TupleRange<f16>>,
    ms1_cycle_range: OptionallyRestricted<TupleRange<MS1CycleIndex>>,
    ms2_cycle_range: OptionallyRestricted<TupleRange<WindowCycleIndex>>,
}

impl QueryRanges {
    /// Compute query ranges from an elution group and tolerance.
    ///
    /// This is the standard conversion used by most aggregators.
    fn from_elution_group<FH: KeyLike, A>(
        aggregator: &A,
        tolerance: &Tolerance,
        rt_ms_to_cycle: impl Fn(u32) -> MS1CycleIndex,
    ) -> Self
    where
        A: HasElutionGroup<FH>,
    {
        let prec_mz_limits = aggregator.elution_group().get_precursor_mz_limits();
        let quad_range =
            tolerance.quad_range_f32((prec_mz_limits.0 as f32, prec_mz_limits.1 as f32));
        let im_range = tolerance.mobility_range_f16(aggregator.elution_group().mobility_ook0());
        let rt_range_milliseconds =
            tolerance.rt_range_as_milis(aggregator.elution_group().rt_seconds());
        let ms1_cycle_range = match rt_range_milliseconds {
            Restricted(x) => Restricted(
                TupleRange::try_new(rt_ms_to_cycle(x.start()), rt_ms_to_cycle(x.end())).unwrap(),
            ),
            Unrestricted => Unrestricted,
        };
        let ms2_cycle_range = ms1_to_ms2_cycle_range(&ms1_cycle_range);

        Self {
            quad_range,
            im_range,
            ms1_cycle_range,
            ms2_cycle_range,
        }
    }

    /// Compute query ranges with RT intersection for ChromatogramCollector.
    ///
    /// Returns None if the RT ranges don't intersect (early termination case).
    fn from_elution_group_with_rt_intersection<FH: KeyLike, A>(
        aggregator: &A,
        tolerance: &Tolerance,
        rt_limits_milis: TupleRange<u32>,
        rt_ms_to_cycle: impl Fn(u32) -> MS1CycleIndex,
    ) -> Option<Self>
    where
        A: HasElutionGroup<FH>,
    {
        let prec_mz_limits = aggregator.elution_group().get_precursor_mz_limits();
        let quad_range =
            tolerance.quad_range_f32((prec_mz_limits.0 as f32, prec_mz_limits.1 as f32));
        let im_range = tolerance.mobility_range_f16(aggregator.elution_group().mobility_ook0());

        let rt_range_milliseconds = match tolerance
            .rt_range_as_milis(aggregator.elution_group().rt_seconds())
            .map(|x| x.try_intercept(rt_limits_milis))
        {
            Restricted(Some(x)) => Restricted(x),
            Restricted(None) => return None,
            Unrestricted => Unrestricted,
        };

        let ms1_cycle_range = match rt_range_milliseconds {
            Restricted(x) => Restricted(
                TupleRange::try_new(rt_ms_to_cycle(x.start()), rt_ms_to_cycle(x.end())).unwrap(),
            ),
            Unrestricted => Unrestricted,
        };
        let ms2_cycle_range = ms1_to_ms2_cycle_range(&ms1_cycle_range);

        Some(Self {
            quad_range,
            im_range,
            ms1_cycle_range,
            ms2_cycle_range,
        })
    }

    /// Constrain the IM range based on quadrupole isolation geometry.
    ///
    /// For MS2 queries, the ion mobility range must be further constrained to the
    /// intersection of the original IM tolerance and the quadrupole isolation window.
    /// This is the repeated pattern in all MS2 query implementations.
    fn constrain_im_for_quadrupole(
        &self,
        quad_info: &timscentroid::geometry::QuadrupoleIsolationScheme,
    ) -> OptionallyRestricted<TupleRange<f16>> {
        self.im_range.map(|x| {
            let tmp = quad_info.intersects_ranges(
                (self.quad_range.start() as f64, self.quad_range.end() as f64),
                (x.start().to_f64(), x.end().to_f64()),
            );
            let tmp =
                tmp.expect("Since we filtered based on the mz+ims there should always be a match");
            (f16::from_f32(tmp.0 as f32), f16::from_f32(tmp.1 as f32))
                .try_into()
                .unwrap()
        })
    }
}

/// Trait to abstract access to the elution group from different aggregator types.
trait HasElutionGroup<FH: KeyLike> {
    fn elution_group(&self) -> &crate::models::elution_group::TimsElutionGroup<FH>;
}

impl<FH: KeyLike> HasElutionGroup<FH> for PointIntensityAggregator<FH> {
    fn elution_group(&self) -> &crate::models::elution_group::TimsElutionGroup<FH> {
        &self.eg
    }
}

impl<FH: KeyLike> HasElutionGroup<FH> for ChromatogramCollector<FH, f32> {
    fn elution_group(&self) -> &crate::models::elution_group::TimsElutionGroup<FH> {
        &self.eg
    }
}

impl<FH: KeyLike, V: Default + crate::traits::ValueLike> HasElutionGroup<FH>
    for SpectralCollector<FH, V>
{
    fn elution_group(&self) -> &crate::models::elution_group::TimsElutionGroup<FH> {
        &self.eg
    }
}

impl<FH: KeyLike> QueriableData<PointIntensityAggregator<FH>> for IndexedTimstofPeaks {
    fn add_query(&self, aggregator: &mut PointIntensityAggregator<FH>, tolerance: &Tolerance) {
        let ranges = QueryRanges::from_elution_group(aggregator, tolerance, |rt| {
            self.rt_ms_to_cycle_index(rt)
        });

        aggregator.eg.iter_precursors().for_each(|(_idx, mz)| {
            let mz_range = tolerance.mz_range_f32(mz as f32);
            self.query_peaks_ms1(mz_range, ranges.ms1_cycle_range, ranges.im_range)
                .for_each(|peak| {
                    aggregator.intensity += peak.intensity as f64;
                });
        });

        self.filter_precursor_ranges(ranges.quad_range, ranges.im_range)
            .for_each(|(quad_info, peaks)| {
                let constrained_im = ranges.constrain_im_for_quadrupole(quad_info);
                aggregator.eg.iter_fragments_refs().for_each(|(_idx, mz)| {
                    let mz_range = tolerance.mz_range_f32(*mz as f32);
                    peaks
                        .query_peaks(mz_range, ranges.ms2_cycle_range, constrained_im)
                        .for_each(|peak| {
                            aggregator.intensity += peak.intensity as f64;
                        });
                });
            })
    }
}

impl<FH: KeyLike> QueriableData<ChromatogramCollector<FH, f32>> for IndexedTimstofPeaks {
    fn add_query(&self, aggregator: &mut ChromatogramCollector<FH, f32>, tolerance: &Tolerance) {
        let Some(ranges) = QueryRanges::from_elution_group_with_rt_intersection(
            aggregator,
            tolerance,
            aggregator.rt_range_milis(),
            |rt| self.rt_ms_to_cycle_index(rt),
        ) else {
            return; // No RT intersection, early exit
        };

        aggregator
            .iter_mut_precursors()
            .for_each(|((_idx, mz), mut chr)| {
                let mz_range = tolerance.mz_range_f32(*mz as f32);
                self.query_peaks_ms1(mz_range, ranges.ms1_cycle_range, ranges.im_range)
                    .for_each(|peak| chr.add_at_index(peak.cycle_index.as_u32(), peak.intensity));
            });

        self.filter_precursor_ranges(ranges.quad_range, ranges.im_range)
            .for_each(|(quad_info, peaks)| {
                let constrained_im = ranges.constrain_im_for_quadrupole(quad_info);
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), mut chr)| {
                        // there is a bit of repeated work here ...
                        // calculating the mz range for each quad instead of once per ion ...
                        // TODO: fix later ...

                        let mz_range = tolerance.mz_range_f32(*mz as f32);
                        peaks
                            .query_peaks(mz_range, ranges.ms2_cycle_range, constrained_im)
                            .for_each(|x| {
                                chr.add_at_index(x.cycle_index.as_u32(), x.intensity);
                            });
                    });
            })
    }
}

impl<FH: KeyLike> QueriableData<SpectralCollector<FH, f32>> for IndexedTimstofPeaks {
    fn add_query(&self, aggregator: &mut SpectralCollector<FH, f32>, tolerance: &Tolerance) {
        let ranges = QueryRanges::from_elution_group(aggregator, tolerance, |rt| {
            self.rt_ms_to_cycle_index(rt)
        });

        aggregator
            .iter_mut_precursors()
            .for_each(|((_idx, mz), ion)| {
                let mz_range = tolerance.mz_range_f32(mz as f32);
                self.query_peaks_ms1(mz_range, ranges.ms1_cycle_range, ranges.im_range)
                    .for_each(|peak| {
                        *ion += peak.intensity;
                    });
            });

        self.filter_precursor_ranges(ranges.quad_range, ranges.im_range)
            .for_each(|(quad_info, peaks)| {
                let constrained_im = ranges.constrain_im_for_quadrupole(quad_info);
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), ion)| {
                        // There is a bit of repeated work here ...
                        // calculating the mz range for each quad instead of once per ion ...
                        // TODO: Fix later ...
                        let mz_range = tolerance.mz_range_f32(*mz as f32);
                        peaks
                            .query_peaks(mz_range, ranges.ms2_cycle_range, constrained_im)
                            .for_each(|x| {
                                *ion += x.intensity;
                            });
                    });
            });
    }
}

impl<FH: KeyLike, V: PeakAddable<MS1CycleIndex> + PeakAddable<WindowCycleIndex>>
    QueriableData<SpectralCollector<FH, V>> for IndexedTimstofPeaks
{
    fn add_query(&self, aggregator: &mut SpectralCollector<FH, V>, tolerance: &Tolerance) {
        let ranges = QueryRanges::from_elution_group(aggregator, tolerance, |rt| {
            self.rt_ms_to_cycle_index(rt)
        });

        aggregator
            .iter_mut_precursors()
            .for_each(|((_idx, mz), ion)| {
                let mz_range = tolerance.mz_range_f32(mz as f32);
                self.query_peaks_ms1(mz_range, ranges.ms1_cycle_range, ranges.im_range)
                    .for_each(|peak| {
                        *ion += *peak;
                    });
            });

        self.filter_precursor_ranges(ranges.quad_range, ranges.im_range)
            .for_each(|(quad_info, peaks)| {
                let constrained_im = ranges.constrain_im_for_quadrupole(quad_info);
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), ion)| {
                        // There is a bit of repeated work here ...
                        // calculating the mz range for each quad instead of once per ion ...
                        // TODO: Fix later ...
                        let mz_range = tolerance.mz_range_f32(*mz as f32);
                        peaks
                            .query_peaks(mz_range, ranges.ms2_cycle_range, constrained_im)
                            .for_each(|x| {
                                *ion += *x;
                            });
                    });
            });
    }
}

fn ms1_to_ms2_cycle_range(
    cycle_range: &OptionallyRestricted<TupleRange<MS1CycleIndex>>,
) -> OptionallyRestricted<TupleRange<WindowCycleIndex>> {
    // FOR NOW we are assuming that the cycle for MS1 will be the same as MS2 ...
    // This is not true in some acquisition schemes but we can fix it later ...
    // This expression is pretty verbose ...
    match cycle_range {
        Restricted(x) => Restricted(
            TupleRange::try_new(
                WindowCycleIndex::new(x.start().as_u32()),
                WindowCycleIndex::new(x.end().as_u32()),
            )
            .unwrap(),
        ),
        Unrestricted => Unrestricted,
    }
}

impl<FH: KeyLike> QueriableData<ChromatogramCollector<FH, f32>> for IndexedPeaksHandle {
    fn add_query(&self, aggregator: &mut ChromatogramCollector<FH, f32>, tolerance: &Tolerance) {
        match self {
            IndexedPeaksHandle::Eager(eager) => {
                // Delegate to existing implementation
                eager.add_query(aggregator, tolerance)
            }
            IndexedPeaksHandle::Lazy(lazy) => {
                // Use QueryRanges to compute common tolerances, then convert to u32 for lazy API
                let Some(ranges) = QueryRanges::from_elution_group_with_rt_intersection(
                    aggregator,
                    tolerance,
                    aggregator.rt_range_milis(),
                    |rt| lazy.rt_ms_to_cycle_index(rt),
                ) else {
                    return; // No RT intersection, early exit
                };

                // Convert MS1CycleIndex ranges to u32 for lazy API
                let cycle_range_u32 = match ranges.ms1_cycle_range {
                    Restricted(x) => Restricted(
                        TupleRange::try_new(x.start().as_u32(), x.end().as_u32()).unwrap(),
                    ),
                    Unrestricted => Unrestricted,
                };

                // Query MS1 precursors
                aggregator
                    .iter_mut_precursors()
                    .for_each(|((_idx, mz), mut chr)| {
                        let mz_range = tolerance.mz_range_f32(*mz as f32);
                        lazy.query_peaks_ms1(mz_range, cycle_range_u32, ranges.im_range)
                            .for_each(|peak| {
                                chr.add_at_index(peak.cycle_index.as_u32(), peak.intensity)
                            });
                    });

                // Query MS2 fragments
                // LazyIndexedTimstofPeaks::query_peaks_ms2 already handles precursor filtering
                aggregator
                    .iter_mut_fragments()
                    .for_each(|((_idx, mz), mut chr)| {
                        let mz_range = tolerance.mz_range_f32(*mz as f32);
                        let results = lazy.query_peaks_ms2(
                            ranges.quad_range,
                            mz_range,
                            cycle_range_u32,
                            ranges.im_range,
                        );

                        // Iterate through all matching window groups and their peaks
                        for (_isolation_scheme, peaks) in results {
                            for peak in peaks {
                                chr.add_at_index(peak.cycle_index.as_u32(), peak.intensity);
                            }
                        }
                    });
            }
        }
    }
}
