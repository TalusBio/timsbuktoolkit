use super::arrays::PartitionedCMGArrayStats;
use crate::errors::{
    DataProcessingError,
    Result,
};
use crate::models::aggregators::rolling_calculators::{
    calculate_lazy_hyperscore,
    calculate_value_vs_baseline,
};
use crate::models::aggregators::streaming_aggregator::RunningStatsCalculator;
use crate::sort_vecs_by_first;
use crate::utils::dtw::dtw_max;
use crate::utils::math::lnfact_float;
use serde::Serialize;
use std::collections::HashMap;
use std::f64;
use std::hash::Hash;
use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};
use tracing::debug;

// fn check_nonnan(x: &[f64]) -> Result<()> {
//     if x.iter().any(|y| y.is_infinite() || y.is_nan()) {
//         let mut out = Vec::with_capacity(x.len());
//         for (i, val) in x.iter().enumerate() {
//             if !(val.is_finite() || val.is_nan()) {
//                 out.push((i, *val));
//             }
//         }
//         Err(crate::errors::DataProcessingError::UnexpectedInfiniteErrors(out).into())
//     } else {
//         Ok(())
//     }
// }

// This name is starting to get really long ...
#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGScoredStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub retention_time_miliseconds: Vec<u32>,
    pub average_mobility: Vec<f64>,
    pub summed_intensity: Vec<u64>,
    pub scores: ChromatogramScores,
    pub transition_mobilities: HashMap<FH, Vec<f64>>,
    pub transition_mzs: HashMap<FH, Vec<f64>>,
    pub expected_mzs: HashMap<FH, f64>,
    pub expected_mobility: f64,
    pub transition_intensities: HashMap<FH, Vec<u64>>,
    pub apex_primary_score_index: usize,
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct ScoresAtTime {
    pub lazyerscore: f64,
    pub lazyerscore_vs_baseline: f64,
    pub norm_lazyerscore_vs_baseline: f64,
    pub cosine_similarity: f64,
    pub npeaks: u8,
    pub retention_time_miliseconds: u32,
    pub average_mobility: f64,
    pub summed_intensity: u64,
    pub mz_errors: Vec<f64>,
    pub mobility_errors: Vec<f64>,
    pub transition_intensities: Vec<u64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ChromatogramScores {
    pub lazy_hyperscore: Vec<f64>,
    pub lazyerscore: Vec<f64>,
    pub lazy_hyperscore_vs_baseline: Vec<f64>,
    pub lazyerscore_vs_baseline: Vec<f64>,
    pub norm_hyperscore_vs_baseline: Vec<f64>,
    pub norm_lazyerscore_vs_baseline: Vec<f64>,
    pub npeaks: Vec<u8>,
    pub cosine_similarity: Option<Vec<f64>>,
}

impl ScoresAtTime {
    pub fn sorted_transition_intensities(&self) -> (Vec<u64>, Vec<f64>, Vec<f64>) {
        let sorted_vecs = sort_vecs_by_first!(
            &self.transition_intensities,
            &self.mobility_errors,
            &self.mz_errors
        );

        let partition_point = sorted_vecs.0.partition_point(|&x| x < 10);
        let mut transition_intensities = sorted_vecs.0[partition_point..].to_vec();
        transition_intensities.reverse();
        let mut mobility_errors = sorted_vecs.2[partition_point..].to_vec();
        mobility_errors.reverse();
        let mut mz_errors = sorted_vecs.1[partition_point..].to_vec();
        mz_errors.reverse();

        (transition_intensities, mobility_errors, mz_errors)
    }
}

impl ChromatogramScores {
    pub fn new<T: Clone + Eq + Serialize + Hash + Send + Sync>(
        arrays: &PartitionedCMGArrayStats<T>,
    ) -> Result<Self> {
        let lazy_hyperscore = calculate_lazy_hyperscore(&arrays.npeaks, &arrays.summed_intensity);
        let basline_window_len = 1 + (arrays.retention_time_miliseconds.len() / 20);
        let lazy_hyperscore_vs_baseline =
            calculate_value_vs_baseline(&lazy_hyperscore, basline_window_len);
        let lazyerscore: Vec<f64> = arrays
            .log_intensity_products
            .iter()
            .map(|x| lnfact_float(*x))
            .collect();
        let lazyerscore_vs_baseline = calculate_value_vs_baseline(&lazyerscore, basline_window_len);

        // Set 0 the NANs
        // Q: Can I do this in-place? Will the comnpiler do it for me?
        let lazy_hyperscore_vs_baseline: Vec<f64> = lazy_hyperscore_vs_baseline
            .into_iter()
            .map(|x| if x.is_nan() { 0.0 } else { x })
            .collect();
        let lazyerscore_vs_baseline: Vec<f64> = lazyerscore_vs_baseline
            .into_iter()
            .map(|x| if x.is_nan() { 0.0 } else { x })
            .collect();

        // Calculate the standard deviation of the lazyscores v baseline
        let mut sd_calculator_hyperscore = RunningStatsCalculator::new(1, 0.0);
        let mut sd_calculator_lazyscore = RunningStatsCalculator::new(1, 0.0);

        (0..lazyerscore_vs_baseline.len()).for_each(|i| {
            sd_calculator_hyperscore.add(lazy_hyperscore_vs_baseline[i], 1);
            sd_calculator_lazyscore.add(lazyerscore_vs_baseline[i], 1);
        });
        let sd_hyperscore = sd_calculator_hyperscore.standard_deviation().unwrap();
        let sd_lazyscore = sd_calculator_lazyscore.standard_deviation().unwrap();

        // Calculate the normalized scores
        let norm_hyperscore_vs_baseline: Vec<f64> = lazy_hyperscore_vs_baseline
            .iter()
            .map(|x| x / sd_hyperscore)
            .collect();
        let norm_lazyerscore_vs_baseline: Vec<f64> = lazyerscore_vs_baseline
            .iter()
            .map(|x| x / sd_lazyscore)
            .collect();

        // Make sure all arrays are the same length
        assert!(
            lazy_hyperscore_vs_baseline.len() == arrays.retention_time_miliseconds.len(),
            "Failed sanity check"
        );
        assert!(
            lazyerscore_vs_baseline.len() == arrays.retention_time_miliseconds.len(),
            "Failed sanity check"
        );

        // check_nonnan(&lazy_hyperscore)?;
        // check_nonnan(&lazyerscore)?;
        // check_nonnan(&lazy_hyperscore_vs_baseline)?;
        // check_nonnan(&lazyerscore_vs_baseline)?;
        // check_nonnan(&norm_hyperscore_vs_baseline)?;
        // check_nonnan(&norm_lazyerscore_vs_baseline)?;

        Ok(ChromatogramScores {
            lazy_hyperscore,
            lazyerscore,
            lazy_hyperscore_vs_baseline,
            lazyerscore_vs_baseline,
            norm_hyperscore_vs_baseline,
            norm_lazyerscore_vs_baseline,
            npeaks: arrays.npeaks.to_owned(),
            cosine_similarity: arrays.cosine_similarity.to_owned(),
        })
    }

    fn get_apex_score_index(&self) -> usize {
        let mut apex_primary_score_index = 0;
        let mut max_primary_score = 0.0f64;
        let primary_scores = &self.norm_lazyerscore_vs_baseline;
        for (i, val) in primary_scores.iter().enumerate() {
            if max_primary_score.is_nan() || *val > max_primary_score {
                max_primary_score = *val;
                apex_primary_score_index = i;
            }
        }
        apex_primary_score_index
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>
    PartitionedCMGScoredStatsArrays<FH>
{
    pub fn index_at_rt(&self, rt_miliseconds: u32) -> usize {
        let part_point = self
            .retention_time_miliseconds
            .partition_point(|x| x < &rt_miliseconds);
        if part_point == 0 { 0 } else { part_point - 1 }
    }

    pub fn scores_at_rt(&self, rt_miliseconds: u32) -> ScoresAtTime {
        let index = self.index_at_rt(rt_miliseconds);
        let lazyerscore = self.scores.lazyerscore[index];
        let lazyerscore_vs_baseline = self.scores.lazyerscore_vs_baseline[index];
        let norm_lazyerscore_vs_baseline = self.scores.norm_lazyerscore_vs_baseline[index];
        let cosine_similarity = match self.scores.cosine_similarity.as_ref() {
            Some(x) => x[index],
            None => 0.0,
        };
        let npeaks = self.scores.npeaks[index];
        let retention_time_miliseconds = self.retention_time_miliseconds[index];
        // TODO decide whether to report errors where the expected val is 0.

        let transition_msz: HashMap<FH, f64> = self
            .transition_mzs
            .iter()
            .map(|x| (x.0.clone(), x.1[index]))
            .collect();
        let transition_mobilities: HashMap<FH, f64> = self
            .transition_mobilities
            .iter()
            .map(|x| (x.0.clone(), x.1[index]))
            .collect();
        let transition_intensities_map: HashMap<FH, u64> = self
            .transition_intensities
            .iter()
            .map(|x| (x.0.clone(), x.1[index]))
            .collect();

        let key_order: Vec<FH> = self.transition_intensities.keys().cloned().collect();
        let mz_errors: Vec<f64> = key_order
            .iter()
            .map(|x| self.expected_mzs[x] - transition_msz[x])
            .collect();
        let mobility_errors: Vec<f64> = key_order
            .iter()
            .map(|x| self.expected_mobility - transition_mobilities[x])
            .collect();
        let transition_intensities: Vec<u64> = key_order
            .iter()
            .map(|x| transition_intensities_map[x])
            .collect();

        let summed_intensity = self.summed_intensity[index];
        let average_mobility = self.average_mobility[index];

        ScoresAtTime {
            lazyerscore,
            lazyerscore_vs_baseline,
            norm_lazyerscore_vs_baseline,
            cosine_similarity,
            npeaks,
            retention_time_miliseconds,
            mz_errors,
            mobility_errors,
            average_mobility,
            summed_intensity,
            transition_intensities,
        }
    }

    /// Calculates the cross score between the two partitioned arrays.
    ///
    /// USES THE OTHER ARRAYS LAZYERSCORE AND PROJECTS ONTO THE SELF SPACE.
    ///
    /// In brief .. it is
    ///     (self.cosine_similarity * other.cosine_similarity) * other.lazyerscore_vs_baseline
    ///
    /// But it projects the 'self' values to the 'other' values.
    pub fn cross_scores<T: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>(
        &self,
        other: &PartitionedCMGScoredStatsArrays<T>,
    ) -> Result<Vec<f64>> {
        let self_cosine = match self.scores.cosine_similarity.as_ref() {
            Some(x) => x,
            None => return Err(DataProcessingError::ExpectedNonEmptyData.into()),
        };
        let self_rt = self.retention_time_miliseconds.as_slice();
        let other_cosine = match other.scores.cosine_similarity.as_ref() {
            Some(x) => x,
            None => return Err(DataProcessingError::ExpectedNonEmptyData.into()),
        };
        let other_rt = other.retention_time_miliseconds.as_slice();
        let other_lzb = other.scores.lazyerscore_vs_baseline.as_slice();

        // Should I project self -> other or other -> self?
        let proj_other_cosine = dtw_max(self_rt, other_rt, other_cosine)?;
        let proj_other_lzb = dtw_max(self_rt, other_rt, other_lzb)?;

        assert!(proj_other_cosine.len() == proj_other_lzb.len());
        assert!(proj_other_cosine.len() == self_rt.len());

        let mut out = Vec::with_capacity(proj_other_cosine.len());
        for (self_cos, (other_cos, other_lzb)) in self_cosine
            .iter()
            .zip(proj_other_cosine.iter().zip(proj_other_lzb.iter()))
        {
            // TODO fix this a layer deeper ...
            // maybe two ... 1. Either prevent needing to wrap the extensions by having
            // The aggregator map to the closest cycle, rather than the frame.
            // 2. Fix the wrapping to not generate -inf values.
            let score = self_cos * other_cos * other_lzb;
            if score.is_infinite() || score.is_nan() {
                out.push(0.0);
                continue;
            }
            out.push(score);
        }

        // Check if there are infinite values
        if out.iter().any(|x| x.is_infinite()) {
            let infinite_indices = out
                .iter()
                .enumerate()
                .filter(|(_, x)| x.is_infinite())
                .collect::<Vec<_>>();
            panic!(
                "Cross scores contain infinite values at indices: \n\n>>> {:?}",
                infinite_indices,
            );
        }

        Ok(out)
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> PartitionedCMGScoredStatsArrays<FH> {
    pub fn new(
        other: PartitionedCMGArrayStats<FH>,
        mz_converter: &Tof2MzConverter,
        mobility_converter: &Scan2ImConverter,
    ) -> Result<Self> {
        let scores = ChromatogramScores::new(&other)?;
        let apex_primary_score_index = scores.get_apex_score_index();
        let average_mobility = other
            .weighted_scan_index_mean
            .into_iter()
            .map(|x| {
                let out = mobility_converter.convert(x);
                if !(0.5..=2.0).contains(&out) {
                    debug!("Bad mobility value: {:?}, input was {:?}", out, x);
                }
                out
            })
            .collect();

        let transition_mobilities = other
            .scan_index_means
            .into_iter()
            .map(|(k, v)| {
                let new_v = v
                    .into_iter()
                    .map(|x| mobility_converter.convert(x))
                    .collect();
                (k, new_v)
            })
            .collect();

        let transition_mzs = other
            .tof_index_means
            .into_iter()
            .map(|(k, v)| {
                let new_v = v.into_iter().map(|x| mz_converter.convert(x)).collect();
                (k, new_v)
            })
            .collect();

        let expected_mobility = mobility_converter.convert(other.expected_scan_index as f64);
        let expected_mzs = other
            .expected_tof_index
            .into_iter()
            .map(|x| (x.0, mz_converter.convert(x.1 as f64)))
            .collect();

        Ok(PartitionedCMGScoredStatsArrays {
            retention_time_miliseconds: other.retention_time_miliseconds,
            average_mobility,
            summed_intensity: other.summed_intensity,
            transition_mobilities,
            transition_mzs,
            transition_intensities: other.intensities,
            apex_primary_score_index,
            scores,
            expected_mzs,
            expected_mobility,
        })
    }
}
