//! Seven-cycle quadratic Savitzky-Golay smoothing on the acquisition cycle axis.

use array2d::Array2D;
use timsquery::models::MzMajorIntensityArray;
use timsquery::traits::KeyLike;

use super::pipeline::TOP_N_FRAGMENTS;

const WINDOW: usize = 7;
const RADIUS: usize = 3;
const CHUNK: usize = 8;
const COSINE_FLOOR: f32 = 1e-3;

const _: () = assert!(WINDOW == 2 * RADIUS + 1);

// Numerators for the orthogonal basis (1, x, x² - 4) on x = -3..=3.
// Check all three reproduction moments at each evaluation point.
const NUMERATORS: [[i16; WINDOW]; WINDOW] = {
    let positions = [-3, -2, -1, 0, 1, 2, 3];
    let mut taps = [[0; WINDOW]; WINDOW];
    let mut sample = 0;
    while sample < WINDOW {
        let t = positions[sample];
        let mut moments = [0; 3];
        let mut j = 0;
        while j < WINDOW {
            let x = positions[j];
            let tap = 12 + 3 * t * x + (t * t - 4) * (x * x - 4);
            taps[sample][j] = tap;
            moments[0] += tap;
            moments[1] += tap * x;
            moments[2] += tap * x * x;
            j += 1;
        }
        assert!(moments[0] == 84);
        assert!(moments[1] == 84 * t);
        assert!(moments[2] == 84 * t * t);
        sample += 1;
    }
    taps
};
const _: () = assert!(NUMERATORS[RADIUS][0] == -8);
const _: () = assert!(NUMERATORS[RADIUS][1] == 12);
const _: () = assert!(NUMERATORS[RADIUS][2] == 24);

#[inline]
fn edge(window: &[f32; WINDOW], sample: usize) -> f32 {
    let baseline = window[sample];
    baseline
        + window
            .iter()
            .zip(NUMERATORS[sample])
            .map(|(&value, tap)| (value - baseline) * (f32::from(tap) / 84.0))
            .sum::<f32>()
}

#[inline]
fn center(window: &[f32], i: usize) -> f32 {
    let baseline = window[i + 3];
    let outer = (window[i] - baseline) + (window[i + 6] - baseline);
    let middle = (window[i + 1] - baseline) + (window[i + 5] - baseline);
    let inner = (window[i + 2] - baseline) + (window[i + 4] - baseline);
    baseline + (-2.0 * outer + 3.0 * middle + 6.0 * inner) * (1.0 / 21.0)
}

fn smooth(input: &[f32], output: &mut [f32]) {
    assert_eq!(output.len(), input.len());
    let Some(first) = input.first_chunk::<WINDOW>() else {
        output.copy_from_slice(input);
        return;
    };
    let last = input.last_chunk::<WINDOW>().expect("checked window length");
    let (left, rest) = output.split_at_mut(RADIUS);
    let (interior, right) = rest.split_at_mut(input.len() - 2 * RADIUS);
    assert_eq!(right.len(), RADIUS);
    for (i, dst) in left.iter_mut().enumerate() {
        *dst = edge(first, i);
    }
    let (chunks, tail) = interior.as_chunks_mut::<CHUNK>();
    for (chunk_index, dst) in chunks.iter_mut().enumerate() {
        let start = chunk_index * CHUNK;
        let window = input[start..]
            .first_chunk::<{ CHUNK + 2 * RADIUS }>()
            .expect("interior chunk has a full halo");
        for (i, dst) in dst.iter_mut().enumerate() {
            *dst = center(window, i);
        }
    }
    let window = &input[chunks.len() * CHUNK..];
    for (i, dst) in tail.iter_mut().enumerate() {
        *dst = center(window, i);
    }
    for (i, dst) in right.iter_mut().enumerate() {
        *dst = edge(last, RADIUS + 1 + i);
    }
}

#[derive(Clone, Copy)]
pub(super) struct ShapeTraces<'a> {
    pub(super) fragments: &'a Array2D<f32>,
    pub(super) precursor: &'a [f32],
    pub(super) cosine: &'a [f32],
}

#[derive(Debug)]
pub(super) struct ShapeScratch {
    fragments: Array2D<f32>,
    precursor: Vec<f32>,
    cosine: Vec<f32>,
    cycle_capacity: u32,
    row_capacity: u32,
    valid_offset: Option<u32>,
}

impl ShapeScratch {
    pub(super) fn new(cycle_capacity: usize, row_capacity: usize) -> Self {
        assert!(row_capacity <= TOP_N_FRAGMENTS);
        Self {
            fragments: Array2D::with_capacity(cycle_capacity, row_capacity),
            precursor: Vec::with_capacity(cycle_capacity),
            cosine: Vec::with_capacity(cycle_capacity),
            cycle_capacity: u32::try_from(cycle_capacity).expect("run cycle count fits u32"),
            row_capacity: u32::try_from(row_capacity).expect("retained row count fits u32"),
            valid_offset: None,
        }
    }

    pub(super) fn invalidate(&mut self) {
        self.valid_offset = None;
    }

    pub(super) fn view<T: KeyLike>(
        &self,
        raw: &MzMajorIntensityArray<T, f32>,
    ) -> Option<ShapeTraces<'_>> {
        let offset = usize::try_from(self.valid_offset?).expect("cycle offset fits usize");
        assert_eq!(offset, raw.cycle_offset);
        assert_eq!(self.fragments.nrows(), raw.num_ions());
        assert_eq!(self.fragments.ncols(), raw.num_cycles());
        assert_eq!(self.precursor.len(), raw.num_cycles());
        assert_eq!(self.cosine.len(), raw.num_cycles());
        Some(ShapeTraces {
            fragments: &self.fragments,
            precursor: &self.precursor,
            cosine: &self.cosine,
        })
    }

    pub(super) fn refill<T: KeyLike>(
        &mut self,
        raw: &MzMajorIntensityArray<T, f32>,
        precursor: &[f32],
        pred_norms: &[(u32, f32)],
        accumulators: [&mut [f32]; 3],
    ) {
        self.invalidate();
        let cycles = raw.num_cycles();
        let rows = raw.num_ions();
        assert!(cycles <= usize::try_from(self.cycle_capacity).expect("cycle capacity fits usize"));
        assert!(rows <= usize::try_from(self.row_capacity).expect("row capacity fits usize"));
        assert_eq!(raw.arr.nrows(), rows);
        assert_eq!(precursor.len(), cycles);
        if cycles < WINDOW {
            return;
        }

        self.fragments.reset_with_value(cycles, rows, 0.0);
        self.precursor.clear();
        self.precursor.resize(cycles, 0.0);
        self.cosine.clear();
        self.cosine.resize(cycles, 0.0);
        assert_eq!(self.fragments.nrows(), raw.arr.nrows());
        assert_eq!(self.fragments.ncols(), raw.arr.ncols());

        smooth(precursor, &mut self.precursor);
        if !self.precursor.iter().all(|value| value.is_finite()) {
            return;
        }
        for (input, dst) in raw.arr.iter_rows().zip(self.fragments.iter_mut_rows()) {
            smooth(input, dst);
            if !dst.iter().all(|value| value.is_finite()) {
                return;
            }
        }

        if fill_cosine(&self.fragments, pred_norms, accumulators, &mut self.cosine) {
            self.valid_offset =
                Some(u32::try_from(raw.cycle_offset).expect("cycle offset fits u32"));
        }
    }
}

fn fill_cosine(
    fragments: &Array2D<f32>,
    pred_norms: &[(u32, f32)],
    [dot, norm_sq, sum]: [&mut [f32]; 3],
    cosine: &mut [f32],
) -> bool {
    let cycles = fragments.ncols();
    assert_eq!(dot.len(), cycles);
    assert_eq!(norm_sq.len(), cycles);
    assert_eq!(sum.len(), cycles);
    assert_eq!(cosine.len(), cycles);
    assert!(pred_norms.len() <= fragments.nrows());
    assert!(pred_norms.windows(2).all(|pair| pair[0].0 < pair[1].0));
    dot.fill(0.0);
    norm_sq.fill(0.0);
    sum.fill(0.0);

    // The raw scribe weights are sqrt(prediction) divided by a common
    // factor. Renormalize that same sparse list to unit L2 norm.
    let pred_norm = pred_norms
        .iter()
        .map(|&(_, weight)| weight * weight)
        .sum::<f32>()
        .sqrt();
    if !pred_norms.is_empty() && (!pred_norm.is_finite() || pred_norm <= 0.0) {
        return false;
    }
    for &(row, weight) in pred_norms {
        let row = usize::try_from(row).expect("retained row index fits usize");
        let row = fragments
            .get_row(row)
            .expect("weight matches a retained row");
        accumulate_cosine(row, weight / pred_norm, dot, norm_sq, sum);
    }

    for (((dst, &dot), &norm_sq), &sum) in cosine
        .iter_mut()
        .zip(dot.iter())
        .zip(norm_sq.iter())
        .zip(sum.iter())
    {
        if !dot.is_finite() || !norm_sq.is_finite() || !sum.is_finite() {
            return false;
        }
        let cosine = if norm_sq > 0.0 {
            (dot / norm_sq.sqrt()).clamp(COSINE_FLOOR, 1.0)
        } else {
            COSINE_FLOOR
        };
        *dst = cosine * cosine * cosine * sum.ln_1p();
    }
    true
}

fn accumulate_cosine(
    row: &[f32],
    weight: f32,
    dot: &mut [f32],
    norm_sq: &mut [f32],
    sum: &mut [f32],
) {
    assert_eq!(row.len(), dot.len());
    assert_eq!(row.len(), norm_sq.len());
    assert_eq!(row.len(), sum.len());
    let (row, row_tail) = row.as_chunks::<CHUNK>();
    let (dot, dot_tail) = dot.as_chunks_mut::<CHUNK>();
    let (norm_sq, norm_tail) = norm_sq.as_chunks_mut::<CHUNK>();
    let (sum, sum_tail) = sum.as_chunks_mut::<CHUNK>();
    for (((row, dot), norm_sq), sum) in row.iter().zip(dot).zip(norm_sq).zip(sum) {
        for i in 0..CHUNK {
            let value = row[i].max(0.0);
            dot[i] += weight * value;
            norm_sq[i] += value * value;
            sum[i] += value;
        }
    }
    for (((&value, dot), norm_sq), sum) in
        row_tail.iter().zip(dot_tail).zip(norm_tail).zip(sum_tail)
    {
        let value = value.max(0.0);
        *dot += weight * value;
        *norm_sq += value * value;
        *sum += value;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproduces_quadratic_basis_at_every_sample() {
        // y = 3 + 2x + x^2 at x = -3..=3: a quadratic fit reproduces it exactly.
        let input = [6.0, 3.0, 2.0, 3.0, 6.0, 11.0, 18.0];

        let mut output = [0.0; WINDOW];
        smooth(&input, &mut output);
        for (actual, expected) in output.into_iter().zip(input) {
            assert!((actual - expected).abs() < 1e-6, "{actual} vs {expected}");
        }
        assert_eq!([output[0], output[3], output[6]], [6.0, 3.0, 18.0]);
    }

    #[test]
    fn preserves_signed_impulse_undershoot_but_clips_it_in_cosine() {
        let input = vec![0.0, 0.0, 0.0, 21.0, 0.0, 0.0, 0.0];
        let raw = MzMajorIntensityArray {
            arr: Array2D::from_flat_vector(input, 1, WINDOW).unwrap(),
            mz_order: vec![(0_u8, 250.0)],
            cycle_offset: 0,
        };
        let mut scratch = ShapeScratch::new(WINDOW, 1);
        let mut accumulators = [[0.0; WINDOW]; 3];
        let [dot, norm_sq, sum] = &mut accumulators;
        scratch.refill(&raw, &[0.0; WINDOW], &[(0, 1.0)], [dot, norm_sq, sum]);

        let shape = scratch.view(&raw).unwrap();
        let row = shape.fragments.get_row(0).unwrap();
        assert_eq!(row[0], -2.0);
        assert_eq!(row[3], 7.0);
        assert_eq!(row[6], -2.0);
        assert_eq!(shape.cosine[0], 0.0);
        assert_eq!(shape.cosine[3], 7.0_f32.ln_1p());
        assert_eq!(shape.cosine[6], 0.0);
    }

    #[test]
    fn leaves_six_cells_exactly_unchanged() {
        let input = [0.0, 3.0, 2.0, 1.0, -1.0, 0.0];
        let mut output = [0.0; 6];
        smooth(&input, &mut output);
        assert_eq!(output, input);
    }

    #[test]
    fn keeps_flat_rows_exact_through_edges_chunks_and_tail() {
        let input = [5.0; WINDOW];
        let mut output = [0.0; WINDOW];
        smooth(&input, &mut output);
        assert_eq!(output, input);

        let input = [f32::MAX; 23];
        let mut output = [0.0; 23];
        smooth(&input, &mut output);
        assert_eq!(output, input);
    }

    #[test]
    fn zero_prediction_excludes_row_from_norm_and_intensity_sum() {
        const CYCLES: usize = 9;
        let raw = MzMajorIntensityArray {
            arr: Array2D::from_flat_vector(
                [vec![3.0; CYCLES], vec![4.0; CYCLES], vec![50.0; CYCLES]].concat(),
                3,
                CYCLES,
            )
            .unwrap(),
            mz_order: vec![(0_u8, 100.0), (1, 200.0), (2, 300.0)],
            cycle_offset: 0,
        };
        let mut scratch = ShapeScratch::new(CYCLES, 3);
        let mut accumulators = [[-1.0; CYCLES]; 3];
        let [dot, norm_sq, sum] = &mut accumulators;
        // Predictions (1, 1, 0) produce these two nonzero raw scribe weights.
        scratch.refill(
            &raw,
            &[0.0; CYCLES],
            &[(0, 0.5), (1, 0.5)],
            [dot, norm_sq, sum],
        );

        let shape = scratch.view(&raw).unwrap();
        assert_eq!(shape.fragments.get_row(2).unwrap(), &[50.0; CYCLES]);
        let cosine = 7.0 / (5.0 * 2.0_f32.sqrt());
        let expected = cosine.powi(3) * 7.0_f32.ln_1p();
        for &value in shape.cosine {
            assert!((value - expected).abs() < 1e-6);
        }
        let [_, norm_sq, sum] = accumulators;
        assert_eq!(norm_sq, [25.0; CYCLES]);
        assert_eq!(sum, [7.0; CYCLES]);
    }
}
