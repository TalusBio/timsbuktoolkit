use super::cv::{
    CrossValidatedScorer,
    DataBuffer,
    FeatureLike,
    GBMConfig,
};
use super::{
    LabelledScore,
    TargetDecoy,
};
use rand::prelude::*;
use rayon::prelude::*;
use tracing::debug;

/// Assign q_values in place.
///
/// # Invariants
/// * `scores` must be sorted in descending order (e.g. best PSM is first)
///
/// Implementation derived from the Sage implementation of qval (Thanks Mike) github.com/lazear/sage
fn assign_qval<T: LabelledScore>(scores: &mut [T], key: impl Fn(&T) -> f32) {
    assert!(
        scores.windows(2).all(|w| key(&w[0]) >= key(&w[1])),
        "Expecting scores to be sorted in descending order",
    );

    let mut decoy = 1;
    let mut target = 0;

    for score in scores.iter_mut() {
        match score.get_label() {
            TargetDecoy::Decoy => decoy += 1,
            TargetDecoy::Target => target += 1,
        }
        score.assign_qval((decoy as f32) / (target as f32));
    }

    // Reverse slice, and calculate the cumulative minimum
    let mut q_min = 1.0f32;
    for score in scores.iter_mut().rev() {
        q_min = q_min.min(score.get_qval());
        score.assign_qval(q_min);
    }

    // We do a third pass to ensure that values with the same score have the same q-value
    let mut last_score = f32::NAN;
    let mut last_qval = 1.0f32;
    for score in scores.iter_mut() {
        let current_score = key(score);
        if current_score != last_score {
            last_score = current_score;
            last_qval = score.get_qval();
            continue;
        }
        score.assign_qval(last_qval);
    }
}

pub fn report_qvalues_at_thresholds<T: LabelledScore + std::fmt::Debug>(
    scores: &[T],
    thresholds: &[f32],
) -> Vec<(f32, usize, usize, usize)> {
    let mut out = Vec::new();

    for &thresh in thresholds {
        let n_below_thresh = scores.iter().filter(|s| s.get_qval() <= thresh).count();
        let n_targets = scores
            .iter()
            .filter(|s| s.get_qval() <= thresh && matches!(s.get_label(), TargetDecoy::Target))
            .count();
        let n_decoys = scores
            .iter()
            .filter(|s| s.get_qval() <= thresh && matches!(s.get_label(), TargetDecoy::Decoy))
            .count();
        out.push((thresh, n_below_thresh, n_targets, n_decoys));
    }

    out
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn rescore<T: LabelledScore + FeatureLike + Send + Sync + std::fmt::Debug>(
    mut data: Vec<T>,
) -> Vec<T> {
    let config = GBMConfig::default();

    data.shuffle(&mut rand::rng());

    let mut scorer = CrossValidatedScorer::new_from_shuffled(5, data, config);
    scorer
        .fit(&mut DataBuffer::default(), &mut DataBuffer::default())
        .unwrap();

    let mut out = scorer.score();
    // Sort by score descending
    out.par_sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    assign_qval(&mut out, |x| T::get_score(x) as f32);
    debug!("Best:\n{:#?}", out.first());
    debug!("Worst:\n{:#?}", out.last());
    out
}

use crate::IonSearchResults;

impl FeatureLike for IonSearchResults {
    fn as_feature(&self) -> impl IntoIterator<Item = f64> + '_ {
        let Self {
            sequence: _,
            library_id: _,
            decoy_group_id: _,
            precursor_mz,
            precursor_charge,
            precursor_mobility_query,
            precursor_rt_query_seconds,
            recalibrated_query_rt,
            nqueries,
            is_target: _,

            // Combined
            main_score,
            delta_next,
            delta_second_next,
            obs_rt_seconds,
            obs_mobility,
            delta_theo_rt,
            sq_delta_theo_rt,
            calibrated_sq_delta_theo_rt,
            delta_ms1_ms2_mobility,
            // ms1_ms2_correlation,
            sq_delta_ms1_ms2_mobility,
            raising_cycles,
            falling_cycles,

            // MS2
            npeaks,
            apex_lazyerscore,
            apex_lazyerscore_vs_baseline,
            apex_norm_lazyerscore_vs_baseline,
            ms2_cosine_ref_similarity,
            ms2_coelution_score,
            ms2_corr_v_gauss,
            ms2_summed_transition_intensity,
            ms2_lazyerscore,
            ms2_isotope_lazyerscore,
            ms2_isotope_lazyerscore_ratio,

            // MS2 - Split
            // Flattening manually bc serde(flatten)
            // is not supported by csv ...
            // https,
            // Q,
            // A,
            ms2_mz_error_0,
            ms2_mz_error_1,
            ms2_mz_error_2,
            ms2_mz_error_3,
            ms2_mz_error_4,
            ms2_mz_error_5,
            ms2_mz_error_6,
            ms2_mobility_error_0,
            ms2_mobility_error_1,
            ms2_mobility_error_2,
            ms2_mobility_error_3,
            ms2_mobility_error_4,
            ms2_mobility_error_5,
            ms2_mobility_error_6,

            // MS1
            ms1_cosine_ref_similarity,
            ms1_coelution_score,
            ms1_summed_precursor_intensity,
            ms1_corr_v_gauss,

            // MS1 Split
            ms1_mz_error_0,
            ms1_mz_error_1,
            ms1_mz_error_2,
            ms1_mobility_error_0,
            ms1_mobility_error_1,
            ms1_mobility_error_2,

            // Relative Intensities
            ms1_inten_ratio_0,
            ms1_inten_ratio_1,
            ms1_inten_ratio_2,

            ms2_inten_ratio_0,
            ms2_inten_ratio_1,
            ms2_inten_ratio_2,
            ms2_inten_ratio_3,
            ms2_inten_ratio_4,
            ms2_inten_ratio_5,
            ms2_inten_ratio_6,

            discriminant_score: _,
            qvalue: _,
            delta_group,
            delta_group_ratio,
        } = *self;

        vec![
            (precursor_mz / 5.0).round(),
            precursor_charge as f64,
            precursor_mobility_query as f64,
            precursor_rt_query_seconds.round() as f64,
            nqueries as f64,
            // Combined
            main_score as f64,
            (main_score / delta_next) as f64,
            delta_next as f64,
            delta_second_next as f64,
            obs_rt_seconds as f64,
            obs_mobility as f64,
            delta_theo_rt as f64,
            sq_delta_theo_rt as f64,
            delta_ms1_ms2_mobility as f64,
            sq_delta_ms1_ms2_mobility as f64,
            raising_cycles as f64,
            falling_cycles as f64,
            // MS2
            npeaks as f64,
            apex_lazyerscore as f64,
            apex_lazyerscore_vs_baseline as f64,
            apex_norm_lazyerscore_vs_baseline as f64,
            ms2_cosine_ref_similarity as f64,
            ms2_coelution_score as f64,
            ms2_corr_v_gauss as f64,
            (ms2_summed_transition_intensity as f64).ln_1p(),
            ms2_lazyerscore as f64,
            ms2_isotope_lazyerscore as f64,
            ms2_isotope_lazyerscore_ratio as f64,
            // MS2 - Split
            ms2_mz_error_0 as f64,
            ms2_mz_error_1 as f64,
            ms2_mz_error_2 as f64,
            ms2_mz_error_3 as f64,
            ms2_mz_error_4 as f64,
            ms2_mz_error_5 as f64,
            ms2_mz_error_6 as f64,
            ms2_mobility_error_0 as f64,
            ms2_mobility_error_1 as f64,
            ms2_mobility_error_2 as f64,
            ms2_mobility_error_3 as f64,
            ms2_mobility_error_4 as f64,
            ms2_mobility_error_5 as f64,
            ms2_mobility_error_6 as f64,
            // MS as f641
            ms1_cosine_ref_similarity as f64,
            ms1_coelution_score as f64,
            (ms1_summed_precursor_intensity as f64).ln_1p(),
            ms1_corr_v_gauss as f64,
            // MS1 Spli as f64t
            ms1_mz_error_0 as f64,
            ms1_mz_error_1 as f64,
            ms1_mz_error_2 as f64,
            ms1_mobility_error_0 as f64,
            ms1_mobility_error_1 as f64,
            ms1_mobility_error_2 as f64,
            // Relative Intensities
            ms1_inten_ratio_0 as f64,
            ms1_inten_ratio_1 as f64,
            ms1_inten_ratio_2 as f64,
            ms2_inten_ratio_0 as f64,
            ms2_inten_ratio_1 as f64,
            ms2_inten_ratio_2 as f64,
            ms2_inten_ratio_3 as f64,
            ms2_inten_ratio_4 as f64,
            ms2_inten_ratio_5 as f64,
            ms2_inten_ratio_6 as f64,
            delta_group as f64,
            delta_group_ratio as f64,
            recalibrated_query_rt as f64,
            calibrated_sq_delta_theo_rt as f64,
        ]
    }

    fn get_y(&self) -> f64 {
        if self.is_target { 1.0 } else { 0.0 }
    }

    fn assign_score(&mut self, score: f64) {
        self.discriminant_score = score as f32;
    }

    fn get_score(&self) -> f64 {
        self.discriminant_score as f64
    }
}

impl LabelledScore for IonSearchResults {
    fn get_label(&self) -> TargetDecoy {
        if self.is_target {
            TargetDecoy::Target
        } else {
            TargetDecoy::Decoy
        }
    }

    fn assign_qval(&mut self, qval: f32) {
        self.qvalue = qval;
    }

    fn get_qval(&self) -> f32 {
        self.qvalue
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // scores = np.array([10, 10, 9, 8, 7, 7, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1])
    // target = np.array([1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0])
    // qvals = np.array([
    //     1 / 4,
    //     1 / 4,
    //     1 / 4,
    //     1 / 4,
    //     2 / 6,
    //     2 / 6,
    //     2 / 6,
    //     3 / 7,
    //     3 / 7,
    //     4 / 7,
    //     5 / 8,
    //     5 / 8,
    //     1,
    //     1,
    //     1,
    //     1,
    // ])

    #[test]
    fn test_assign_qval() {
        let scores = vec![10, 10, 9, 8, 7, 7, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1];
        let target = vec![1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0];
        let qvals = vec![
            1. / 4.,
            1. / 4.,
            1. / 4.,
            1. / 4.,
            2. / 6.,
            2. / 6.,
            2. / 6.,
            3. / 7.,
            3. / 7.,
            4. / 8.,   // 4. / 7.,
            4. / 8.,   // 5. / 8.,
            4. / 8.,   // 5. / 8.,
            6.0 / 8.0, // 1.,
            6.0 / 8.0, // 1.,
            6.0 / 8.0, // 1.,
            6.0 / 8.0, // 1.,
        ];
        struct TestScore {
            score: f64,
            label: TargetDecoy,
            qval: f32,
        }

        impl LabelledScore for TestScore {
            fn get_label(&self) -> TargetDecoy {
                self.label
            }

            fn assign_qval(&mut self, q: f32) {
                self.qval = q
            }

            fn get_qval(&self) -> f32 {
                self.qval
            }
        }

        let mut data = scores
            .iter()
            .zip(target.iter())
            .map(|(&s, &l)| TestScore {
                score: s as f64,
                label: if l == 1 {
                    TargetDecoy::Target
                } else {
                    TargetDecoy::Decoy
                },
                qval: 0.0,
            })
            .collect::<Vec<_>>();

        assign_qval(&mut data, |x| x.score as f32);

        for i in 0..qvals.len() {
            let real = qvals[i];
            let model = data[i].get_qval();

            assert!(
                (real - model).abs() < 1e-6,
                "At index {}: expected {}, got {}; REAL: {:?} MODEL {:?}",
                i,
                real,
                model,
                qvals,
                data.iter().map(|x| x.get_qval()).collect::<Vec<_>>(),
            );
        }
    }
}
