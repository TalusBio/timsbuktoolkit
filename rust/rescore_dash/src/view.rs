//! The single type crossing from a rescoring run into the dashboard.
//!
//! Deliberately free of any pipeline type: anything that can fill a
//! [`RescoreView`] can drive the TUI, which is what keeps a future standalone
//! viewer a new constructor rather than a new dashboard.

use std::sync::Arc;

/// Per-fold feature importance and stats, as produced by the rescorer.
///
/// A local plain struct rather than `timsseek`'s `RescoreFeatureStats` so this
/// crate stays free of a `timsseek` dependency; the caller converts.
pub struct FoldImportance {
    pub fold: u8,
    /// `(feature name, GBM gain)`
    pub gain: Vec<(Arc<str>, f32)>,
    /// `(feature name, mean, NaN ratio)`
    pub stats: Vec<(Arc<str>, f32, f32)>,
}

/// One rescoring run's rows: the model's feature matrix plus labels, scores and
/// q-values.
///
/// `features` is row-major and borrowed — the caller owns the matrix (it can be
/// hundreds of MB) and lends it for the dashboard's lifetime.
pub struct RescoreView<'a> {
    /// ALL-lane feature names, in model column order.
    pub feature_names: Vec<Arc<str>>,
    /// Row-major, `n_rows * n_features` long.
    pub features: &'a [f64],
    pub is_target: Vec<bool>,
    /// `discriminant_score`, row-aligned.
    pub score: Vec<f32>,
    pub qvalue: Vec<f32>,
    pub importance: Vec<FoldImportance>,
}

#[derive(Debug)]
pub enum ViewError {
    /// `features.len()` is not `n_rows * n_features`.
    MatrixLen { expected: usize, got: usize },
    /// `score`/`qvalue` disagree with `is_target` on the row count.
    RowLen {
        field: &'static str,
        expected: usize,
        got: usize,
    },
    /// Nothing to show.
    Empty,
}

impl std::fmt::Display for ViewError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MatrixLen { expected, got } => {
                write!(f, "feature matrix has {got} values, expected {expected}")
            }
            Self::RowLen {
                field,
                expected,
                got,
            } => {
                write!(f, "`{field}` has {got} rows, expected {expected}")
            }
            Self::Empty => write!(f, "no rows to display"),
        }
    }
}

impl std::error::Error for ViewError {}

impl<'a> RescoreView<'a> {
    pub fn n_rows(&self) -> usize {
        self.is_target.len()
    }

    pub fn n_features(&self) -> usize {
        self.feature_names.len()
    }

    /// Values of feature `j` in row order. Panics only on an out-of-range `j`,
    /// which is a programming error, not user data.
    pub fn feature_column(&self, j: usize) -> impl Iterator<Item = f64> + '_ {
        assert!(j < self.n_features(), "feature index out of range");
        self.features
            .iter()
            .skip(j)
            .step_by(self.n_features())
            .copied()
    }

    /// Mean GBM gain for `name` across folds; `0.0` for a feature no fold
    /// reported (the LDA path reports coefficients for the linear lane only).
    pub fn mean_gain(&self, name: &str) -> f32 {
        let mut sum = 0.0f32;
        let mut n = 0u32;
        for fold in &self.importance {
            for (nm, gain) in &fold.gain {
                if &**nm == name {
                    sum += gain;
                    n += 1;
                }
            }
        }
        if n == 0 { 0.0 } else { sum / n as f32 }
    }

    pub fn validate(&self) -> Result<(), ViewError> {
        let rows = self.n_rows();
        if rows == 0 || self.n_features() == 0 {
            return Err(ViewError::Empty);
        }
        for (field, got) in [("score", self.score.len()), ("qvalue", self.qvalue.len())] {
            if got != rows {
                return Err(ViewError::RowLen {
                    field,
                    expected: rows,
                    got,
                });
            }
        }
        let expected = rows * self.n_features();
        if self.features.len() != expected {
            return Err(ViewError::MatrixLen {
                expected,
                got: self.features.len(),
            });
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(n: usize) -> Vec<Arc<str>> {
        (0..n)
            .map(|i| Arc::from(format!("f{i}").as_str()))
            .collect()
    }

    #[test]
    fn feature_column_walks_the_right_stride() {
        // 3 rows x 2 features, row-major.
        let matrix = [1.0, 10.0, 2.0, 20.0, 3.0, 30.0];
        let view = RescoreView {
            feature_names: names(2),
            features: &matrix,
            is_target: vec![true, false, true],
            score: vec![0.0; 3],
            qvalue: vec![1.0; 3],
            importance: Vec::new(),
        };
        assert_eq!(view.n_rows(), 3);
        assert_eq!(view.n_features(), 2);
        assert_eq!(
            view.feature_column(0).collect::<Vec<_>>(),
            vec![1.0, 2.0, 3.0]
        );
        assert_eq!(
            view.feature_column(1).collect::<Vec<_>>(),
            vec![10.0, 20.0, 30.0]
        );
    }

    #[test]
    fn validate_rejects_a_ragged_matrix() {
        let matrix = [1.0, 2.0, 3.0];
        let view = RescoreView {
            feature_names: names(2),
            features: &matrix,
            is_target: vec![true, false],
            score: vec![0.0; 2],
            qvalue: vec![1.0; 2],
            importance: Vec::new(),
        };
        assert!(matches!(view.validate(), Err(ViewError::MatrixLen { .. })));
    }

    #[test]
    fn validate_rejects_mismatched_row_vectors() {
        let matrix = [1.0, 2.0, 3.0, 4.0];
        let view = RescoreView {
            feature_names: names(2),
            features: &matrix,
            is_target: vec![true, false],
            score: vec![0.0; 1],
            qvalue: vec![1.0; 2],
            importance: Vec::new(),
        };
        assert!(matches!(view.validate(), Err(ViewError::RowLen { .. })));
    }

    #[test]
    fn mean_gain_averages_folds_and_defaults_to_zero() {
        let matrix = [1.0, 2.0];
        let view = RescoreView {
            feature_names: names(2),
            features: &matrix,
            is_target: vec![true],
            score: vec![0.0],
            qvalue: vec![1.0],
            importance: vec![
                FoldImportance {
                    fold: 0,
                    gain: vec![(Arc::from("f0"), 2.0)],
                    stats: Vec::new(),
                },
                FoldImportance {
                    fold: 1,
                    gain: vec![(Arc::from("f0"), 4.0)],
                    stats: Vec::new(),
                },
            ],
        };
        assert_eq!(view.mean_gain("f0"), 3.0);
        assert_eq!(view.mean_gain("f1"), 0.0);
    }
}
