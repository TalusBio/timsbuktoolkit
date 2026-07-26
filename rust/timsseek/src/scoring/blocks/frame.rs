//! `FeatFrame` + `FrameSink` — name-bound, row-major-transposable feature
//! matrix.
//!
//! A column store (`Vec<Arc<str>>` names + `Vec<Vec<f64>>` data) that a
//! [`FrameSink`] fills row by row, binding each pushed value to its column name
//! so values and names cannot desync. This is the ML consumers' feature matrix;
//! [`super::NameSink`] builds the same names set-level, without any record.

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
/// Mirrors `super::ColSink::slot`'s first-row-vs-append idiom.
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
    /// (order must be identical every row) via a `&str` comparison, with no
    /// allocation. Mirrors `ColSink::slot` exactly: `cursor == cols.len()` is
    /// "never seen this column before", with no separate row counter needed;
    /// the owned `Arc<str>` is only created in the create branch, once per
    /// column, ever.
    ///
    /// Callers that already hold the final name use this; callers that would
    /// have to *build* one (the `{prefix}_{i}` array fan-out) use
    /// `FrameSink::slot_lazy` instead, which never builds it on rows 1..N.
    fn slot(&mut self, name: &str, v: f64) {
        let i = self.cursor;
        if i == self.frame.cols.len() {
            self.open_column(Arc::from(name), v);
        } else {
            debug_assert_eq!(&*self.frame.names[i], name, "column order mismatch");
            self.frame.cols[i].push(v);
        }
        self.cursor += 1;
    }

    /// [`FrameSink::slot`] for names that have to be *constructed*: `name` is
    /// invoked only when a column is actually created (the first row) and, in
    /// debug builds, for the column-order check. On rows 1..N of a release
    /// build the closure is never called, so the name costs nothing.
    ///
    /// `Fn`, not `FnOnce`, precisely so the debug order check can also call it.
    fn slot_lazy(&mut self, name: impl Fn() -> String, v: f64) {
        let i = self.cursor;
        if i == self.frame.cols.len() {
            self.open_column(Arc::from(name().as_str()), v);
        } else {
            #[cfg(debug_assertions)]
            debug_assert_eq!(
                &*self.frame.names[i],
                name().as_str(),
                "column order mismatch"
            );
            self.frame.cols[i].push(v);
        }
        self.cursor += 1;
    }

    /// Create the column at the cursor: record its name and open a
    /// `nrows`-sized buffer holding `v`. Does not touch the cursor — the
    /// `slot*` caller advances it.
    fn open_column(&mut self, name: Arc<str>, v: f64) {
        self.frame.names.push(name);
        let mut col = Vec::with_capacity(self.nrows);
        col.push(v);
        self.frame.cols.push(col);
    }

    pub fn push(&mut self, name: &str, v: f64) {
        self.slot(name, v);
    }

    /// `x.ln_1p()` under the already-suffixed `name`. Every `push_*` below
    /// takes the FINAL column name and only applies its arithmetic: the
    /// suffix table lives in exactly one place, the macro's
    /// `Generator::name_suffix`, which hands both the value walk and the name
    /// walk the same compile-time literal.
    pub fn push_ln1p(&mut self, name: &str, x: f64) {
        self.slot(name, x.ln_1p());
    }

    /// `x.log2()` under the already-suffixed `name` (see
    /// [`FrameSink::push_ln1p`]).
    pub fn push_log2(&mut self, name: &str, x: f64) {
        self.slot(name, x.log2());
    }

    /// `x.round()` under the already-suffixed `name` (see
    /// [`FrameSink::push_ln1p`]).
    pub fn push_round(&mut self, name: &str, x: f64) {
        self.slot(name, x.round());
    }

    /// Magnitude fold (`|x|`, NaN preserved) under the already-suffixed
    /// `name` (see [`FrameSink::push_ln1p`]).
    pub fn push_abs(&mut self, name: &str, x: f64) {
        self.slot(name, x.abs());
    }

    /// Missingness indicator (`1.0` if `x` is non-finite, else `0.0`) under
    /// the already-suffixed `name` (see [`FrameSink::push_ln1p`]).
    pub fn push_isna(&mut self, name: &str, x: f64) {
        self.slot(name, if x.is_finite() { 0.0 } else { 1.0 });
    }

    /// One bare column per element, named `{prefix}_{i}`. The index is a
    /// runtime value, so this is the one naming the macro cannot hand over as
    /// a literal; `FrameSink::slot_lazy` keeps it off the per-row path.
    pub fn push_slice(&mut self, prefix: &str, vals: &[f32]) {
        for (i, v) in vals.iter().enumerate() {
            self.slot_lazy(|| format!("{prefix}_{i}"), *v as f64);
        }
    }

    /// [`FrameSink::push_isna`] for each element of an array field: one column
    /// per element, named `{prefix}_{i}_isna` — again built lazily, see
    /// [`FrameSink::push_slice`].
    pub fn push_slice_isna(&mut self, prefix: &str, vals: &[f32]) {
        for (i, v) in vals.iter().enumerate() {
            self.slot_lazy(
                || format!("{prefix}_{i}_isna"),
                if v.is_finite() { 0.0 } else { 1.0 },
            );
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
            s.push_ln1p("y_ln1p", 0.0);
            s.begin_row();
            s.push("x", 20.0);
            s.push_ln1p("y_ln1p", std::f64::consts::E - 1.0);
            s.finish();
        }
        assert_eq!(f.names(), &[Arc::from("x"), Arc::from("y_ln1p")]);
        assert_eq!(f.column(0), &[10.0, 20.0]);
        assert!((f.column(1)[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn frame_sink_push_log2() {
        let mut f = FeatFrame::with_capacity(1, 1);
        {
            let mut s = FrameSink::new(&mut f, 1);
            s.begin_row();
            s.push_log2("x_log2", 8.0);
            s.finish();
        }
        assert_eq!(f.names(), &[Arc::from("x_log2")]);
        assert_eq!(f.column(0), &[3.0]);
    }

    #[test]
    fn frame_sink_push_round() {
        let mut f = FeatFrame::with_capacity(1, 1);
        {
            let mut s = FrameSink::new(&mut f, 1);
            s.begin_row();
            s.push_round("x_round", 2.6);
            s.finish();
        }
        assert_eq!(f.names(), &[Arc::from("x_round")]);
        assert_eq!(f.column(0), &[3.0]);
    }

    #[test]
    fn frame_sink_push_abs() {
        let mut f = FeatFrame::with_capacity(1, 2);
        {
            let mut s = FrameSink::new(&mut f, 2);
            s.begin_row();
            s.push_abs("x_abs", -4.5);
            s.begin_row();
            s.push_abs("x_abs", f64::NAN);
            s.finish();
        }
        assert_eq!(f.names(), &[Arc::from("x_abs")]);
        assert_eq!(f.column(0)[0], 4.5);
        assert!(f.column(0)[1].is_nan());
    }

    #[test]
    fn frame_sink_push_isna() {
        let mut f = FeatFrame::with_capacity(1, 2);
        {
            let mut s = FrameSink::new(&mut f, 2);
            s.begin_row();
            s.push_isna("x_isna", 1.0);
            s.begin_row();
            s.push_isna("x_isna", f64::NAN);
            s.finish();
        }
        assert_eq!(f.names(), &[Arc::from("x_isna")]);
        assert_eq!(f.column(0), &[0.0, 1.0]);
    }

    #[test]
    fn frame_sink_push_slice() {
        let mut f = FeatFrame::with_capacity(3, 1);
        {
            let mut s = FrameSink::new(&mut f, 1);
            s.begin_row();
            s.push_slice("v", &[1.0f32, 2.5f32, 3.0f32]);
            s.finish();
        }
        assert_eq!(
            f.names(),
            &[Arc::from("v_0"), Arc::from("v_1"), Arc::from("v_2")]
        );
        assert_eq!(f.column(0), &[1.0]);
        assert_eq!(f.column(1), &[2.5]);
        assert_eq!(f.column(2), &[3.0]);
    }

    #[test]
    fn frame_sink_push_slice_isna() {
        let mut f = FeatFrame::with_capacity(2, 1);
        {
            let mut s = FrameSink::new(&mut f, 1);
            s.begin_row();
            s.push_slice_isna("v", &[1.0f32, f32::NAN]);
            s.finish();
        }
        assert_eq!(f.names(), &[Arc::from("v_0_isna"), Arc::from("v_1_isna")]);
        assert_eq!(f.column(0), &[0.0]);
        assert_eq!(f.column(1), &[1.0]);
    }

    /// The array fan-out builds its `{prefix}_{i}` names through a closure that
    /// only the create branch (row 0) invokes, so rows 1..N never construct a
    /// name. Across several rows the names must still be exactly the row-0 ones
    /// and every value must land in its own column, in order.
    #[test]
    fn frame_sink_slice_names_come_from_first_row_only() {
        let rows = [
            [1.0f32, 2.0f32],
            [3.0f32, 4.0f32],
            [5.0f32, f32::NAN],
            [7.0f32, 8.0f32],
        ];
        let mut f = FeatFrame::with_capacity(4, rows.len());
        {
            let mut s = FrameSink::new(&mut f, rows.len());
            for r in &rows {
                s.begin_row();
                s.push_slice("v", r);
                s.push_slice_isna("v", r);
            }
            s.finish();
        }
        assert_eq!(
            f.names(),
            &[
                Arc::from("v_0"),
                Arc::from("v_1"),
                Arc::from("v_0_isna"),
                Arc::from("v_1_isna"),
            ]
        );
        assert_eq!(f.nrows(), 4);
        assert_eq!(f.column(0), &[1.0, 3.0, 5.0, 7.0]);
        assert_eq!(f.column(2), &[0.0, 0.0, 0.0, 0.0]);
        assert_eq!(f.column(3), &[0.0, 0.0, 1.0, 0.0]);
        // NaN does not compare equal, so column 1 is checked elementwise.
        assert_eq!(f.column(1)[0], 2.0);
        assert_eq!(f.column(1)[1], 4.0);
        assert!(f.column(1)[2].is_nan());
        assert_eq!(f.column(1)[3], 8.0);
    }

    /// Scalar `push_*` take the FINAL name (the macro already appended the
    /// suffix) and must use it verbatim, not append anything of their own.
    #[test]
    fn scalar_pushes_use_the_name_verbatim() {
        let mut f = FeatFrame::with_capacity(5, 3);
        {
            let mut s = FrameSink::new(&mut f, 3);
            for i in 0..3 {
                let x = (i + 1) as f64;
                s.begin_row();
                s.push("a", x);
                s.push_log2("b_log2", x);
                s.push_ln1p("c_ln1p", x);
                s.push_abs("d_abs", -x);
                s.push_isna("e_isna", x);
            }
            s.finish();
        }
        assert_eq!(
            f.names(),
            &[
                Arc::from("a"),
                Arc::from("b_log2"),
                Arc::from("c_ln1p"),
                Arc::from("d_abs"),
                Arc::from("e_isna"),
            ]
        );
        assert_eq!(f.column(0), &[1.0, 2.0, 3.0]);
        assert_eq!(f.column(1)[1], 1.0);
        assert_eq!(f.column(3), &[1.0, 2.0, 3.0]);
        assert_eq!(f.column(4), &[0.0, 0.0, 0.0]);
    }

    #[test]
    #[should_panic]
    fn push_column_wrong_length_panics() {
        let mut f = FeatFrame::with_capacity(1, 3);
        f.push_column("a", vec![1.0, 2.0]);
    }

    #[test]
    #[should_panic]
    fn extend_mismatched_nrows_panics() {
        let mut f = FeatFrame::with_capacity(1, 3);
        f.push_column("a", vec![1.0, 2.0, 3.0]);

        let mut g = FeatFrame::with_capacity(1, 2);
        g.push_column("b", vec![4.0, 5.0]);

        f.extend(g);
    }
}
