pub use forust_ml::constraints::{
    Constraint,
    ConstraintMap,
};
pub use forust_ml::errors::ForustError;
pub use forust_ml::gradientbooster::{
    GrowPolicy,
    MissingNodeTreatment,
};
pub use forust_ml::metric::{
    EvaluationMetric,
    Metric,
};
pub use forust_ml::objective::ObjectiveType;
pub use forust_ml::sampler::SampleMethod;
pub use forust_ml::{
    GradientBooster,
    Matrix,
};
pub use rayon::prelude::*;
pub use std::collections::HashSet;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TargetDecoy {
    Target,
    Decoy,
}

pub trait LabelledScore {
    fn get_score(&self) -> f64;
    fn get_label(&self) -> TargetDecoy;
    fn assign_qval(&mut self, q: f32);
    fn get_qval(&self) -> f32;
}

impl LabelledScore for (f64, TargetDecoy, f32) {
    fn get_score(&self) -> f64 {
        self.0
    }

    fn get_label(&self) -> TargetDecoy {
        self.1
    }

    fn assign_qval(&mut self, q: f32) {
        self.2 = q
    }

    fn get_qval(&self) -> f32 {
        self.2
    }
}

/// Assign q_values in place.
///
/// # Invariants
/// * `scores` must be sorted in descending order (e.g. best PSM is first)
///
/// Implementation derived from the Sage implementation of qval (Thanks Mike) github.com/lazear/sage
fn assign_qval<T: LabelledScore>(scores: &mut [T]) {
    // FDR Calculation:
    // * Sort by score, descending
    // * Estimate FDR
    // * Calculate q-value
    //
    let first_score = scores.first().unwrap();
    let last_score = scores.last().unwrap();
    assert!(first_score.get_score() >= last_score.get_score());

    let mut decoy = 1;
    let mut target = 0;

    for score in scores.iter_mut() {
        match score.get_label() {
            TargetDecoy::Decoy => decoy += 1,
            TargetDecoy::Target => target += 1,
        }
        score.assign_qval(decoy as f32 / target as f32);
    }

    // Reverse slice, and calculate the cumulative minimum
    let mut q_min = 1.0f32;
    for score in scores.iter_mut().rev() {
        q_min = q_min.min(score.get_qval());
        score.assign_qval(q_min);
    }
}

pub struct GBMConfig {
    iterations: usize,
    learning_rate: f32,
    max_depth: usize,
    max_leaves: usize,
    l1: f32,
    l2: f32,
    gamma: f32,
    max_delta_step: f32,
    min_leaf_weight: f32,
    base_score: f64,
    nbins: u16,
    parallel: bool,
    allow_missing_splits: bool,
    monotone_constraints: Option<ConstraintMap>,
    subsample: f32,
    top_rate: f64,
    other_rate: f64,
    colsample_bytree: f64,
    seed: u64,
    missing: f64,
    create_missing_branch: bool,
    sample_method: SampleMethod,
    grow_policy: GrowPolicy,
    evaluation_metric: Option<Metric>,
    early_stopping_rounds: Option<usize>,
    initialize_base_score: bool,
    terminate_missing_features: HashSet<usize>,
    missing_node_treatment: MissingNodeTreatment,
    log_iterations: usize,
    force_children_to_bound_parent: bool,
}

impl Clone for GBMConfig {
    fn clone(&self) -> Self {
        let Self {
            iterations,
            learning_rate,
            max_depth,
            max_leaves,
            l1,
            l2,
            gamma,
            max_delta_step,
            min_leaf_weight,
            base_score,
            nbins,
            parallel,
            allow_missing_splits,
            monotone_constraints,
            subsample,
            top_rate,
            other_rate,
            colsample_bytree,
            seed,
            missing,
            create_missing_branch,
            sample_method,
            grow_policy,
            evaluation_metric,
            early_stopping_rounds,
            initialize_base_score,
            terminate_missing_features,
            missing_node_treatment,
            log_iterations,
            force_children_to_bound_parent,
        } = self;

        Self {
            iterations: *iterations,
            learning_rate: *learning_rate,
            max_depth: *max_depth,
            max_leaves: *max_leaves,
            l1: *l1,
            l2: *l2,
            gamma: *gamma,
            max_delta_step: *max_delta_step,
            min_leaf_weight: *min_leaf_weight,
            base_score: *base_score,
            nbins: *nbins,
            parallel: *parallel,
            allow_missing_splits: *allow_missing_splits,
            monotone_constraints: monotone_constraints.clone(),
            subsample: *subsample,
            top_rate: *top_rate,
            other_rate: *other_rate,
            colsample_bytree: *colsample_bytree,
            seed: *seed,
            missing: *missing,
            create_missing_branch: *create_missing_branch,
            sample_method: match sample_method {
                SampleMethod::None => SampleMethod::None,
                SampleMethod::Random => SampleMethod::Random,
                SampleMethod::Goss => SampleMethod::Goss,
            },
            grow_policy: match grow_policy {
                GrowPolicy::DepthWise => GrowPolicy::DepthWise,
                GrowPolicy::LossGuide => GrowPolicy::LossGuide,
            },
            evaluation_metric: *evaluation_metric,
            early_stopping_rounds: *early_stopping_rounds,
            initialize_base_score: *initialize_base_score,
            terminate_missing_features: terminate_missing_features.clone(),
            missing_node_treatment: *missing_node_treatment,
            log_iterations: *log_iterations,
            force_children_to_bound_parent: *force_children_to_bound_parent,
        }
    }
}

impl Default for GBMConfig {
    fn default() -> Self {
        GBMConfig {
            iterations: 500,
            learning_rate: 0.1,
            max_depth: 5,
            max_leaves: usize::MAX,
            l1: 0.,
            l2: 1.,
            gamma: 0.,
            max_delta_step: 0.,
            min_leaf_weight: 5.,
            base_score: 0.5,
            nbins: 256,
            parallel: true,
            allow_missing_splits: true,
            monotone_constraints: None,
            subsample: 0.75,
            top_rate: 0.1,
            other_rate: 0.2,
            colsample_bytree: 1.0,
            seed: 0,
            missing: f64::NAN,
            create_missing_branch: false,
            sample_method: SampleMethod::None,
            grow_policy: GrowPolicy::DepthWise,
            evaluation_metric: None,
            early_stopping_rounds: Some(10),
            initialize_base_score: true,
            terminate_missing_features: HashSet::new(),
            missing_node_treatment: MissingNodeTreatment::AssignToParent,
            log_iterations: 10,
            force_children_to_bound_parent: false,
        }
    }
}

impl GBMConfig {
    fn try_build(&self) -> Result<GradientBooster, ForustError> {
        let Self {
            iterations,
            learning_rate,
            max_depth,
            max_leaves,
            l1,
            l2,
            gamma,
            max_delta_step,
            min_leaf_weight,
            base_score,
            nbins,
            parallel,
            allow_missing_splits,
            monotone_constraints,
            subsample,
            top_rate,
            other_rate,
            colsample_bytree,
            seed,
            missing,
            create_missing_branch,
            sample_method,
            grow_policy,
            evaluation_metric,
            early_stopping_rounds,
            initialize_base_score,
            terminate_missing_features,
            missing_node_treatment,
            log_iterations,
            force_children_to_bound_parent,
        } = self;

        GradientBooster::new(
            ObjectiveType::LogLoss,
            *iterations,
            *learning_rate,
            *max_depth,
            *max_leaves,
            *l1,
            *l2,
            *gamma,
            *max_delta_step,
            *min_leaf_weight,
            *base_score,
            *nbins,
            *parallel,
            *allow_missing_splits,
            monotone_constraints.clone(),
            *subsample,
            *top_rate,
            *other_rate,
            *colsample_bytree,
            *seed,
            *missing,
            *create_missing_branch,
            match sample_method {
                SampleMethod::None => SampleMethod::None,
                SampleMethod::Random => SampleMethod::Random,
                SampleMethod::Goss => SampleMethod::Goss,
            },
            match grow_policy {
                GrowPolicy::DepthWise => GrowPolicy::DepthWise,
                GrowPolicy::LossGuide => GrowPolicy::LossGuide,
            },
            *evaluation_metric,
            *early_stopping_rounds,
            *initialize_base_score,
            terminate_missing_features.clone(),
            *missing_node_treatment,
            *log_iterations,
            *force_children_to_bound_parent,
        )
    }
}

pub trait FeatureLike<const N: usize> {
    fn as_feature(&self) -> [f64; N];
    fn get_y(&self) -> f64;
}

#[derive(Default)]
pub struct DataBuffer<const N: usize> {
    fold_buffer: Vec<f64>,
    accum_buffer: Vec<[f64; N]>,
    response_buffer: Vec<f64>,
}

impl<const N: usize> DataBuffer<N> {
    fn fill_buffer(&mut self, assigned_fold: &[u8], data: &[impl FeatureLike<N>], fold: u8) {
        self.fold_buffer.clear();
        self.accum_buffer.clear();
        self.response_buffer.clear();

        for (elem_fold, elem) in assigned_fold.iter().zip(data.iter()) {
            if fold == *elem_fold {
                self.accum_buffer.push(elem.as_feature());
                self.response_buffer.push(elem.get_y());
            }
        }

        for feat_idx in (0..N).into_iter() {
            for e in self.accum_buffer.iter() {
                // Maybe its faster to make the vec with zeros and assign to the
                // positions ... so we iterate only once per elemnent instead
                // of once per element*feature
                let value = e[feat_idx];
                self.fold_buffer.push(value);
            }
        }

        assert_eq!(self.accum_buffer.len(), self.response_buffer.len());
        assert!(self.accum_buffer.len() > 0, "No data for fold {}", fold);
    }

    fn as_matrix(&self) -> (Matrix<'_, f64>, &'_ [f64]) {
        let nrows = self.accum_buffer.len();
        (
            Matrix::new(self.fold_buffer.as_slice(), nrows, N),
            self.response_buffer.as_slice(),
        )
    }
}

/// This is meant to IN ESSENCE ...
///
/// Provided we have a number of splits k >= 3.
/// We would train k classifiers, for for classifier n,
/// we would use as train data fold n, as early stop data
/// fold n + 1 and thus we can safely use as test/scoring
/// data any other fold.
///
/// So the score for any point in the data is the average of
/// the results for all classifiers that didint use it
/// for either training or early_stopping_rounds.
pub struct CrossValidatedScorer<const N: usize, T: FeatureLike<N>> {
    n_folds: u8,
    data: Vec<T>,
    assigned_fold: Vec<u8>,
    fold_classifiers: Vec<Option<GradientBooster>>,
    // I tried this but makes no difference ...
    // fold_classifiers: Vec<Option<SelfSupervisedBooster>>,
    config: GBMConfig,
}

impl<const N: usize, T: FeatureLike<N>> CrossValidatedScorer<N, T> {
    /// Create a new CrossValidatedScorer
    ///
    /// NOTE: THIS ASSUMES YOUR DATA IS ALREADY SHUFFLED
    /// FOLDS WILL BE ASSIGNED IN ORDER (0, 1, 2, ..., n_folds-1, 0, 1, ...)
    /// IF YOUR DATA IS ORDERED IN ANY WAY, COULD LEAD TO BIASED RESULTS.
    pub fn new(n_folds: u8, data: Vec<T>, config: GBMConfig) -> Self {
        let assigned_fold: Vec<u8> = (0..data.len())
            .into_iter()
            .map(|x| (x % n_folds as usize).try_into().unwrap())
            .collect();
        Self {
            n_folds,
            data,
            assigned_fold,
            fold_classifiers: Vec::new(),
            config,
        }
    }

    pub fn fit<'a>(
        &mut self,
        train_buffer: &'a mut DataBuffer<N>,
        val_buffer: &'a mut DataBuffer<N>,
    ) -> Result<(), ForustError> {
        self.fold_classifiers.clear();
        // 3 folds == [0, 1, 2]
        (0..self.n_folds)
            .into_iter()
            .for_each(|_| self.fold_classifiers.push(None));
        for fold in (0..self.n_folds).into_iter() {
            self.fit_fold(fold, train_buffer, val_buffer)?
        }
        Ok(())
    }

    pub fn score(&self) -> Vec<f64> {
        let mut scores = vec![0.0; self.data.len()];
        let mut buffer = DataBuffer::default();

        for train_i in (0..self.n_folds).into_iter() {
            let train_i = train_i as u8;
            let early_stop_i = self.next_fold(train_i);

            for inference_i in (0..self.n_folds).into_iter() {
                let inference_i = inference_i as u8;
                if inference_i == train_i {
                    continue;
                };
                if inference_i == early_stop_i {
                    continue;
                };
                let scorer = self.fold_classifiers[train_i as usize].as_ref().unwrap();
                let (matrix, _res) = self.fold_to_matrix(inference_i, &mut buffer);
                let preds = scorer.predict(&matrix, true);
                // Now we need to add the predictions to the scores
                let mut score_i = 0;
                for pred in preds.into_iter() {
                    while self.assigned_fold[score_i] != inference_i {
                        score_i += 1;
                    }
                    scores[score_i] += pred;
                    score_i += 1;
                }
            }
        }

        let div_factor = (self.n_folds - 2) as f64;
        scores.iter_mut().for_each(|x| {
            *x /= div_factor;
        });

        scores
    }

    fn next_fold(&self, fold: u8) -> u8 {
        let mut maybe_next = fold + 1;
        if maybe_next >= self.n_folds {
            maybe_next -= self.n_folds;
        }
        assert!(maybe_next < self.n_folds);
        maybe_next
    }

    fn fold_to_matrices<'a>(
        &'a self,
        fold: u8,
        train_buffer: &'a mut DataBuffer<N>,
        val_buffer: &'a mut DataBuffer<N>,
    ) -> ((Matrix<'a, f64>, &'a [f64]), (Matrix<'a, f64>, &'a [f64])) {
        let next_fold_id = self.next_fold(fold);
        (
            self.fold_to_matrix(fold, train_buffer),
            self.fold_to_matrix(next_fold_id, val_buffer),
        )
    }

    fn fold_to_matrix<'a>(
        &self,
        fold: u8,
        buffer: &'a mut DataBuffer<N>,
    ) -> (Matrix<'a, f64>, &'a [f64]) {
        buffer.fill_buffer(self.assigned_fold.as_slice(), self.data.as_slice(), fold);
        buffer.as_matrix()
    }

    fn fit_fold<'a>(
        &mut self,
        fold: u8,
        train_buffer: &'a mut DataBuffer<N>,
        val_buffer: &'a mut DataBuffer<N>,
    ) -> Result<(), ForustError> {
        // let mut model = SelfSupervisedBooster::try_new(&self.config, 10)?;
        let mut model = self.config.try_build()?;
        let ((matrix, response), (v_matrix, v_response)) =
            self.fold_to_matrices(fold, train_buffer, val_buffer);
        // let (val_matrix, val_response) = self.fold_to_matrix(fold);
        // model.fit_unweighted(matrix, response, evaluation_data)
        let evaluation_data = Some(vec![(v_matrix, v_response, &[1.0f64; 1][..])]);
        model.fit_unweighted(&matrix, response, evaluation_data)?;
        // model
        //     .fit_unweighted(matrix, response, v_matrix, v_response)
        //     .map_err(|x| x.0)?;
        self.fold_classifiers[fold as usize] = Some(model);
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::distr::{
        Distribution,
        Uniform,
    };

    use crate::ml::FeatureLike;

    struct MyFeature {
        vals: [f64; 5],
        class: f64,
    }

    impl FeatureLike<5> for MyFeature {
        fn as_feature(&self) -> [f64; 5] {
            self.vals
        }

        fn get_y(&self) -> f64 {
            self.class
        }
    }

    fn random_data(num_targets: usize, num_decoys: usize) -> Vec<MyFeature> {
        let between_unch = Uniform::try_from(1.0..10.0).unwrap();
        let between_targ = Uniform::try_from(1.0..20.0).unwrap();
        let mut rng = rand::rng();
        let mut out = Vec::new();

        for nt in 0..num_targets {
            let tmp = MyFeature {
                vals: [
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                    between_targ.sample(&mut rng),
                    between_targ.sample(&mut rng),
                ],
                class: 1.0,
            };
            out.push(tmp);
        }
        for nt in 0..num_decoys {
            let tmp = MyFeature {
                vals: [
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                    between_unch.sample(&mut rng),
                ],
                class: 0.0,
            };
            out.push(tmp);
        }
        out
    }

    #[test]
    fn test_cv_train() {
        let config = GBMConfig::default();
        let data = random_data(500, 500);
        let data_len = data.len();

        let mut scorer = CrossValidatedScorer::new(3, data, config);
        scorer
            .fit(&mut DataBuffer::default(), &mut DataBuffer::default())
            .unwrap();

        let out = scorer.score();
        let num_t_gt0 = out[..=500].iter().filter(|&&x| x > 0.0).count();
        let num_d_gt0 = out[500..].iter().filter(|&&x| x > 0.0).count();
        // There are 2 features that have 2x the uniform range.
        // Thus 75% chance of having at least one feature where the value is out of the possible
        // ranges in the decoys ... (500 * 0.75) = 375
        // I can crunch more formal stats bc ... prob distributions ...
        // but this should be a safe enough value to test.
        assert!(num_t_gt0 > 375, "num_t_gt0: {}", num_t_gt0);
        assert!(num_d_gt0 < 40, "num_d_gt0: {}", num_d_gt0);
        assert_eq!(out.len(), data_len);
    }
}
