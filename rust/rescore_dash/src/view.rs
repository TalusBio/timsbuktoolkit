//! The single type crossing from a rescoring run into the dashboard.
//!
//! Deliberately free of any pipeline type: anything that can fill a
//! [`RescoreView`] can drive the TUI, which is what keeps a future standalone
//! viewer a new constructor rather than a new dashboard.

use std::sync::Arc;

/// How many rows pass at one q-value cutoff.
///
/// An unlabeled row is not counted at all, rather than counted into a total
/// alone, which is why there is no separate total to pass in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThresholdRow {
    pub q: f32,
    pub targets: usize,
    pub decoys: usize,
}

impl ThresholdRow {
    pub(crate) fn total(&self) -> usize {
        self.targets + self.decoys
    }
}

/// One rescoring run's rows: the model's feature matrix plus labels, scores,
/// q-values and per-feature gain.
///
/// Borrowed slices and no owned state: the caller keeps the data and lends it
/// for the duration of [`crate::precompute::Dashboard`] construction. After
/// that the view is dead and the matrix can go.
pub struct RescoreView<'a> {
    /// ALL-lane feature names, in model column order.
    pub feature_names: &'a [Arc<str>],
    /// Row-major, `n_rows * n_features` long.
    pub features: &'a [f64],
    pub is_target: &'a [bool],
    /// `discriminant_score`, row-aligned.
    pub score: &'a [f32],
    pub qvalue: &'a [f32],
    /// Rows passing at each q-value cutoff, tightest first.
    ///
    /// Passed in rather than counted here: these are the numbers the caller
    /// already put in its run log, and a panel that disagreed with the log
    /// about the headline ID count would be worse than no panel.
    pub thresholds: &'a [ThresholdRow],
    /// Fold-averaged GBM gain, aligned to `feature_names`.
    ///
    /// Averaged by the caller: the dashboard wants one number per feature, and
    /// column-aligned makes reading it an array index rather than a search.
    pub gain: &'a [f32],
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
    /// `gain` is not aligned to `feature_names`.
    GainLen { expected: usize, got: usize },
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
            Self::GainLen { expected, got } => {
                write!(f, "`gain` has {got} entries, expected {expected} features")
            }
            Self::Empty => write!(f, "no rows to display"),
        }
    }
}

impl std::error::Error for ViewError {}

impl RescoreView<'_> {
    pub fn n_rows(&self) -> usize {
        self.is_target.len()
    }

    pub fn n_features(&self) -> usize {
        self.feature_names.len()
    }

    /// Row `i`'s feature values.
    pub fn row(&self, i: usize) -> &[f64] {
        let nf = self.n_features();
        &self.features[i * nf..(i + 1) * nf]
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
        if self.gain.len() != self.n_features() {
            return Err(ViewError::GainLen {
                expected: self.n_features(),
                got: self.gain.len(),
            });
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
    fn row_walks_the_right_stride() {
        // 3 rows x 2 features, row-major.
        let matrix = [1.0, 10.0, 2.0, 20.0, 3.0, 30.0];
        let names = names(2);
        let view = RescoreView {
            feature_names: &names,
            features: &matrix,
            is_target: &[true, false, true],
            score: &[0.0; 3],
            qvalue: &[1.0; 3],
            thresholds: &[],
            gain: &[0.0; 2],
        };
        assert_eq!(view.n_rows(), 3);
        assert_eq!(view.n_features(), 2);
        assert_eq!(view.row(0), &[1.0, 10.0]);
        assert_eq!(view.row(2), &[3.0, 30.0]);
        view.validate().expect("well-formed");
    }

    #[test]
    fn validate_rejects_a_ragged_matrix() {
        let names = names(2);
        let view = RescoreView {
            feature_names: &names,
            features: &[1.0, 2.0, 3.0],
            is_target: &[true, false],
            score: &[0.0; 2],
            qvalue: &[1.0; 2],
            thresholds: &[],
            gain: &[0.0; 2],
        };
        assert!(matches!(view.validate(), Err(ViewError::MatrixLen { .. })));
    }

    #[test]
    fn validate_rejects_mismatched_row_vectors() {
        let names = names(2);
        let view = RescoreView {
            feature_names: &names,
            features: &[1.0, 2.0, 3.0, 4.0],
            is_target: &[true, false],
            score: &[0.0; 1],
            qvalue: &[1.0; 2],
            thresholds: &[],
            gain: &[0.0; 2],
        };
        assert!(matches!(view.validate(), Err(ViewError::RowLen { .. })));
    }

    /// `gain` is indexed by feature column with no name lookup, so a
    /// misaligned slice would silently attribute one feature's importance to
    /// another rather than failing.
    #[test]
    fn validate_rejects_gain_that_is_not_column_aligned() {
        let names = names(2);
        let view = RescoreView {
            feature_names: &names,
            features: &[1.0, 2.0, 3.0, 4.0],
            is_target: &[true, false],
            score: &[0.0; 2],
            qvalue: &[1.0; 2],
            thresholds: &[],
            gain: &[0.0; 1],
        };
        assert!(matches!(view.validate(), Err(ViewError::GainLen { .. })));
    }

    #[test]
    fn validate_rejects_an_empty_view() {
        let view = RescoreView {
            feature_names: &[],
            features: &[],
            is_target: &[],
            score: &[],
            qvalue: &[],
            thresholds: &[],
            gain: &[],
        };
        assert!(matches!(view.validate(), Err(ViewError::Empty)));
    }
}
