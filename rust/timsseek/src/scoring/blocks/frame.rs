//! `FeatFrame` + `FrameSink` — name-bound, row-major-transposable feature
//! matrix.
//!
//! Additive primitive for the `#[derive(ScoreBlock)]` migration: a column
//! store (`Vec<Arc<str>>` names + `Vec<Vec<f64>>` data) that a [`FrameSink`]
//! fills row by row, binding each pushed value to its column name so the two
//! can never desync (unlike [`super::FeatSink`]/[`super::NameSink`], which are
//! built and walked separately and rely on the caller keeping order in sync).
//! Nothing consumes this yet.

use std::sync::Arc;

// ---------------------------------------------------------------------------
// FeatFrame — column-major feature store
// ---------------------------------------------------------------------------

/// A row-major-convertible, name-bound feature matrix: one `Arc<str>` name and
/// one `Vec<f64>` column per feature, all columns sharing `nrows`.
#[derive(Debug, Default, Clone)]
pub struct FeatFrame {
    names: Vec<Arc<str>>,
    cols: Vec<Vec<f64>>,
    nrows: usize,
}

impl FeatFrame {
    /// Pre-allocate `ncols_hint` columns and `nrows_hint` capacity in each.
    pub fn with_capacity(ncols_hint: usize, nrows_hint: usize) -> Self {
        Self {
            names: Vec::with_capacity(ncols_hint),
            cols: Vec::with_capacity(ncols_hint),
            nrows: nrows_hint,
        }
    }

    pub fn names(&self) -> &[Arc<str>] {
        &self.names
    }

    pub fn ncols(&self) -> usize {
        self.cols.len()
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }

    pub fn column(&self, j: usize) -> &[f64] {
        &self.cols[j]
    }

    /// Append one column. `values.len()` must equal `self.nrows`.
    pub fn push_column(&mut self, name: impl Into<Arc<str>>, values: Vec<f64>) {
        assert_eq!(
            values.len(),
            self.nrows,
            "column length must match frame nrows"
        );
        self.names.push(name.into());
        self.cols.push(values);
    }

    /// Append all of `other`'s columns + names. `other.nrows()` must equal
    /// `self.nrows()`.
    pub fn extend(&mut self, other: FeatFrame) {
        assert_eq!(
            other.nrows, self.nrows,
            "cannot extend frame with mismatched nrows"
        );
        self.names.extend(other.names);
        self.cols.extend(other.cols);
    }

    /// Transpose the column store into a row-major `Vec<f64>`, so
    /// `out[i*ncols + j]` is row `i`, column `j` — the layout LDA/forust want.
    pub fn row_major(&self) -> Vec<f64> {
        let ncols = self.ncols();
        let mut out = vec![0.0; self.nrows * ncols];
        for (j, col) in self.cols.iter().enumerate() {
            for (i, v) in col.iter().enumerate() {
                out[i * ncols + j] = *v;
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// FrameSink — row-oriented, name-bound accumulator for FeatFrame
// ---------------------------------------------------------------------------

/// Fills a [`FeatFrame`] row by row, binding each value to its column name at
/// push time. On the first row (`row == 0`) each `push*` call creates its
/// column (name + `Vec::with_capacity(nrows)`); on later rows it walks a
/// per-row cursor over the already-established columns,
/// `debug_assert_eq!`ing that the name at the cursor matches, and appends.
/// Mirrors [`super::ColSink::slot`]'s first-row-vs-append idiom.
pub struct FrameSink<'f> {
    frame: &'f mut FeatFrame,
    nrows: usize,
    cursor: usize,
}

impl<'f> FrameSink<'f> {
    pub fn new(frame: &'f mut FeatFrame, nrows: usize) -> Self {
        frame.nrows = nrows;
        Self {
            frame,
            nrows,
            cursor: 0,
        }
    }

    /// Start a new row: reset the per-row column cursor.
    pub fn begin_row(&mut self) {
        self.cursor = 0;
    }

    /// Locate the column at the cursor (creating it the first time the
    /// cursor reaches it — i.e. on the first row — with a
    /// `Vec::with_capacity(nrows)` buffer), advance the cursor, and append
    /// `v`. On later rows this instead checks the name at the cursor matches
    /// (order must be identical every row). Mirrors `ColSink::slot` exactly:
    /// `cursor == cols.len()` is "never seen this column before", with no
    /// separate row counter needed.
    fn slot(&mut self, name: String, v: f64) {
        let i = self.cursor;
        if i == self.frame.cols.len() {
            self.frame.names.push(Arc::from(name));
            let mut col = Vec::with_capacity(self.nrows);
            col.push(v);
            self.frame.cols.push(col);
        } else {
            debug_assert_eq!(
                self.frame.names[i].as_ref(),
                name.as_str(),
                "column order mismatch"
            );
            self.frame.cols[i].push(v);
        }
        self.cursor += 1;
    }

    pub fn push(&mut self, name: &str, v: f64) {
        self.slot(name.to_string(), v);
    }

    pub fn push_ln1p(&mut self, name: &str, x: f64) {
        self.slot(format!("{name}_ln1p"), x.ln_1p());
    }

    pub fn push_log2(&mut self, name: &str, x: f64) {
        self.slot(format!("{name}_log2"), x.log2());
    }

    pub fn push_round(&mut self, name: &str, x: f64) {
        self.slot(format!("{name}_round"), x.round());
    }

    /// Magnitude fold (`|x|`, NaN preserved).
    pub fn push_abs(&mut self, name: &str, x: f64) {
        self.slot(format!("{name}_abs"), x.abs());
    }

    /// Missingness indicator (`1.0` if `x` is non-finite, else `0.0`).
    pub fn push_isna(&mut self, name: &str, x: f64) {
        self.slot(
            format!("{name}_isna"),
            if x.is_finite() { 0.0 } else { 1.0 },
        );
    }

    /// One bare column per element, named `{prefix}_{i}`.
    pub fn push_slice(&mut self, prefix: &str, vals: &[f32]) {
        for (i, v) in vals.iter().enumerate() {
            self.slot(format!("{prefix}_{i}"), *v as f64);
        }
    }

    /// Consume the sink. The row cursor's final row count is tracked by the
    /// caller via `nrows` passed to [`FrameSink::new`]; this only exists to
    /// give callers an explicit place to end the borrow.
    pub fn finish(self) {}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_column_and_extend() {
        let mut f = FeatFrame::with_capacity(2, 3);
        f.push_column("a", vec![1.0, 2.0, 3.0]);
        f.push_column("b", vec![4.0, 5.0, 6.0]);
        assert_eq!(f.ncols(), 2);
        assert_eq!(f.nrows(), 3);
        assert_eq!(f.names(), &[Arc::from("a"), Arc::from("b")]);
        assert_eq!(f.column(1), &[4.0, 5.0, 6.0]);

        let mut g = FeatFrame::with_capacity(1, 3);
        g.push_column("c", vec![7.0, 8.0, 9.0]);
        f.extend(g);
        assert_eq!(f.ncols(), 3);
        assert_eq!(f.names()[2], Arc::from("c"));
    }

    #[test]
    fn row_major_transposes() {
        let mut f = FeatFrame::with_capacity(2, 2);
        f.push_column("a", vec![1.0, 3.0]);
        f.push_column("b", vec![2.0, 4.0]);
        // row-major: [row0: a,b][row1: a,b] = [1,2,3,4]
        assert_eq!(f.row_major(), vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn frame_sink_binds_name_and_value() {
        let mut f = FeatFrame::with_capacity(1, 2);
        {
            let mut s = FrameSink::new(&mut f, 2);
            s.begin_row();
            s.push("x", 10.0);
            s.push_ln1p("y", 0.0);
            s.begin_row();
            s.push("x", 20.0);
            s.push_ln1p("y", std::f64::consts::E - 1.0);
            s.finish();
        }
        assert_eq!(f.names(), &[Arc::from("x"), Arc::from("y_ln1p")]);
        assert_eq!(f.column(0), &[10.0, 20.0]);
        assert!((f.column(1)[1] - 1.0).abs() < 1e-12);
    }
}
