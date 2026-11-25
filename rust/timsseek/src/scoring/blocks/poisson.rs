//! Fragment counts under a fitted signal-plus-background model.

use timsquery::IonAnnot;
use timsquery::models::MzMajorIntensityArray;
use timsseek_macros::ScoreBlock;

use crate::data_sources::reference_library::MAX_FRAGMENTS;
use crate::scoring::pipeline::TOP_N_FRAGMENTS;

const MAX_HALF_WIDTH: u32 = 9;
const MAX_ON_CYCLES: usize = (2 * MAX_HALF_WIDTH + 1) as usize;
const OFF_GUARD_CYCLES: u32 = 2;
const OFF_FLANK_CYCLES: u32 = 4;

#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct PoissonEvidence {
    #[feat(ln1p)]
    pub ms2_poisson_llr: f64,
    #[feat(ln1p)]
    pub ms2_poisson_signal_to_background: f64,
    #[feat(raw)]
    pub ms2_poisson_missing_fragments: f32,
}

impl PoissonEvidence {
    pub(crate) const MISSING: Self = Self {
        ms2_poisson_llr: f64::NAN,
        ms2_poisson_signal_to_background: f64::NAN,
        ms2_poisson_missing_fragments: f32::NAN,
    };
}

#[derive(Clone, Copy)]
struct Fragment<'a> {
    on_row: Option<&'a [f32]>,
    fraction: f64,
    background: f64,
    on_is_zero: bool,
}

pub(crate) fn compute_poisson_evidence(
    apex_profile: &[f32],
    apex: usize,
    rise_and_fall: (u8, u8),
    observed: &MzMajorIntensityArray<IonAnnot, f32>,
    predicted: impl Iterator<Item = (IonAnnot, f32)>,
    gain: f64,
) -> PoissonEvidence {
    assert!(gain.is_finite() && gain > 0.0);
    assert_eq!(apex_profile.len(), observed.num_cycles());
    assert!(observed.num_ions() <= TOP_N_FRAGMENTS);
    assert!(apex < apex_profile.len());
    let cycle_count = u32::try_from(apex_profile.len()).expect("cycle count fits u32");
    let apex = u32::try_from(apex).expect("apex index fits u32");
    let (rise, fall) = rise_and_fall;

    // Rise and fall both include the apex.
    let on_start = apex.saturating_sub(u32::from(rise.saturating_sub(1)).min(MAX_HALF_WIDTH));
    let on_end = apex
        .saturating_add(u32::from(fall.saturating_sub(1)).min(MAX_HALF_WIDTH) + 1)
        .min(cycle_count);
    let left_end = on_start.saturating_sub(OFF_GUARD_CYCLES);
    let left_start = left_end.saturating_sub(OFF_FLANK_CYCLES);
    let right_start = on_end.saturating_add(OFF_GUARD_CYCLES).min(cycle_count);
    let right_end = right_start
        .saturating_add(OFF_FLANK_CYCLES)
        .min(cycle_count);
    assert!(left_start <= left_end && left_end <= on_start && on_start <= apex);
    assert!(apex < on_end && on_end <= right_start && right_start <= right_end);
    assert!(right_end <= cycle_count);
    let on_len = on_end - on_start;
    let off_len = (left_end - left_start) + (right_end - right_start);
    assert!(on_len as usize <= MAX_ON_CYCLES && off_len <= 2 * OFF_FLANK_CYCLES);
    if off_len == 0 {
        return PoissonEvidence::MISSING;
    }

    let on = on_start as usize..on_end as usize;
    let left = left_start as usize..left_end as usize;
    let right = right_start as usize..right_end as usize;

    let top = top_predicted(predicted);
    let top = &top[..top.partition_point(|(_, intensity)| *intensity > 0.0)];
    let prediction_sum: f64 = top.iter().map(|(_, intensity)| f64::from(*intensity)).sum();
    let mut shape = [0.0; MAX_ON_CYCLES];
    let shape = &mut shape[..on_len as usize];
    for (weight, &value) in shape.iter_mut().zip(&apex_profile[on.clone()]) {
        *weight = f64::from(value.max(0.0));
    }
    let shape_sum: f64 = shape.iter().sum();
    if prediction_sum == 0.0 || shape_sum == 0.0 {
        return PoissonEvidence::MISSING;
    }
    assert!(shape_sum.is_finite());
    shape.iter_mut().for_each(|weight| *weight /= shape_sum);
    let shape_sum = shape.iter().sum::<f64>();
    assert!((shape_sum - 1.0).abs() <= 1e-5, "shape_sum={shape_sum}");

    let background_floor = 0.5 / f64::from(on_len.max(off_len));
    let mut fragments = [Fragment {
        on_row: None,
        fraction: 0.0,
        background: 0.0,
        on_is_zero: true,
    }; TOP_N_FRAGMENTS];
    let fragments = &mut fragments[..top.len()];
    let mut on_sum = 0.0;
    let mut background_sum = 0.0;
    for (fragment, &(key, intensity)) in fragments.iter_mut().zip(top) {
        // Top predictions absent after extraction pruning stay in the model as zero rows.
        let row = observed.get_row(&key);
        let off_sum = row.map_or(0.0, |row| {
            sum_intensity(&row[left.clone()]) + sum_intensity(&row[right.clone()])
        });
        let on_row = row.map(|row| &row[on.clone()]);
        let row_sum = on_row.map_or(0.0, sum_intensity);
        let background = (off_sum / gain / f64::from(off_len)).max(background_floor);
        assert!(background.is_finite() && background > 0.0);
        *fragment = Fragment {
            on_row,
            fraction: f64::from(intensity) / prediction_sum,
            background,
            on_is_zero: row_sum == 0.0,
        };
        on_sum += row_sum;
        background_sum += background;
    }
    let fraction_sum = fragments.iter().map(|row| row.fraction).sum::<f64>();
    assert!(
        (fraction_sum - 1.0).abs() <= 1e-5,
        "fraction_sum={fraction_sum}"
    );
    fit_signal(
        shape,
        fragments,
        on_sum / gain,
        background_sum,
        gain.recip(),
    )
}

fn top_predicted(
    predicted: impl Iterator<Item = (IonAnnot, f32)>,
) -> [(IonAnnot, f32); TOP_N_FRAGMENTS] {
    let mut top = [(IonAnnot::default(), 0.0); TOP_N_FRAGMENTS];
    for (index, (key, intensity)) in predicted.enumerate() {
        assert!(index < MAX_FRAGMENTS);
        if intensity.is_finite() && intensity > top[TOP_N_FRAGMENTS - 1].1 {
            let rank = top.partition_point(|(_, existing)| *existing >= intensity);
            top.copy_within(rank..TOP_N_FRAGMENTS - 1, rank + 1);
            top[rank] = (key, intensity);
        }
    }
    top
}

fn sum_intensity(row: &[f32]) -> f64 {
    let (chunks, remainder) = row.as_chunks::<4>();
    let mut sums = [0.0; 4];
    for chunk in chunks {
        for (sum, &value) in sums.iter_mut().zip(chunk) {
            *sum += f64::from(value.max(0.0));
        }
    }
    sums.into_iter().sum::<f64>()
        + remainder
            .iter()
            .map(|&value| f64::from(value.max(0.0)))
            .sum::<f64>()
}

fn fit_signal(
    shape: &[f64],
    fragments: &[Fragment<'_>],
    total_count: f64,
    background_sum: f64,
    inv_gain: f64,
) -> PoissonEvidence {
    assert!(total_count.is_finite() && total_count >= 0.0);
    assert!(background_sum.is_finite() && background_sum > 0.0);
    let on_len = shape.len() as f64;
    let seed = (total_count - on_len * background_sum).max(0.0);
    let at_seed = evaluate_signal::<true>(shape, fragments, inv_gain, seed);

    // The derivative is nonpositive above total_count, since each n*q <= n/S.
    let trial = if at_seed.information > 0.0 {
        (seed + at_seed.gradient / at_seed.information).clamp(0.0, total_count)
    } else {
        seed
    };
    let trial_llr = if trial == seed {
        at_seed.llr
    } else {
        evaluate_signal::<false>(shape, fragments, inv_gain, trial).llr
    };
    let (signal, llr) = if trial_llr > 0.0 {
        (trial, trial_llr)
    } else if at_seed.llr > 0.0 {
        (seed, at_seed.llr)
    } else {
        (0.0, 0.0)
    };

    let missing = fragments
        .iter()
        .filter(|row| row.on_is_zero && signal * row.fraction + on_len * row.background > 3.0)
        .count();
    let ratio = signal / background_sum;
    assert!((0.0..=total_count).contains(&signal));
    assert!(llr.is_finite() && llr >= 0.0 && ratio.is_finite() && ratio >= 0.0);
    assert!(missing <= fragments.len() && fragments.len() <= TOP_N_FRAGMENTS);
    PoissonEvidence {
        ms2_poisson_llr: llr,
        ms2_poisson_signal_to_background: ratio,
        ms2_poisson_missing_fragments: missing as f32,
    }
}

struct SignalEvaluation {
    llr: f64,
    gradient: f64,
    information: f64,
}

fn evaluate_signal<const NEWTON: bool>(
    shape: &[f64],
    fragments: &[Fragment<'_>],
    inv_gain: f64,
    signal: f64,
) -> SignalEvaluation {
    // These terms cover the change in mean across every cell, even zero rows.
    let mut half_llr = -signal;
    let mut gradient = -1.0;
    let mut information = 0.0;
    if NEWTON || signal > 0.0 {
        for fragment in fragments.iter().filter(|row| !row.on_is_zero) {
            let row = fragment
                .on_row
                .expect("positive count requires an observed row");
            let signal_over_bg = signal * fragment.fraction / fragment.background;
            let mut accumulate = |&value: &f32, &time_weight: &f64| {
                let count = f64::from(value.max(0.0)) * inv_gain;
                half_llr += count * (signal_over_bg * time_weight).ln_1p();
                if NEWTON {
                    let weight = time_weight * fragment.fraction;
                    let slope = weight / (fragment.background + signal * weight);
                    gradient += count * slope;
                    information += count * slope * slope;
                }
            };
            let (row_chunks, row_tail) = row.as_chunks::<4>();
            let (shape_chunks, shape_tail) = shape.as_chunks::<4>();
            for (values, weights) in row_chunks.iter().zip(shape_chunks) {
                for (value, weight) in values.iter().zip(weights) {
                    accumulate(value, weight);
                }
            }
            for (value, weight) in row_tail.iter().zip(shape_tail) {
                accumulate(value, weight);
            }
        }
    }
    SignalEvaluation {
        llr: 2.0 * half_llr,
        gradient,
        information,
    }
}

/// Reusable five-cycle window for gain measurements, in fragment-major order.
pub struct CalibrantGainScratch {
    cells: [[f32; 5]; MAX_FRAGMENTS],
}

impl CalibrantGainScratch {
    pub fn new() -> Self {
        Self {
            cells: [[0.0; 5]; MAX_FRAGMENTS],
        }
    }

    /// Pool the squared fourth differences and their mean-intensity weights
    /// across fragments. The contrast removes locally cubic mean profiles.
    ///
    /// # Panics
    /// Panics if the collector exceeds the validated library fragment limit.
    pub fn measure(
        &mut self,
        fragments: &MzMajorIntensityArray<IonAnnot, f32>,
        apex: usize,
    ) -> Option<f64> {
        let start = apex.checked_sub(2)?;
        let end = start.checked_add(5)?;
        if end > fragments.num_cycles() {
            return None;
        }
        assert!(fragments.num_ions() <= MAX_FRAGMENTS);
        let cells = &mut self.cells[..fragments.num_ions()];
        for (window, (_, row)) in cells.iter_mut().zip(fragments.iter_mzs()) {
            window.copy_from_slice(&row[start..end]);
        }

        let mut squared_difference = 0.0;
        let mut weighted_intensity = 0.0;
        for window in cells {
            if !window
                .iter()
                .all(|&value| value.is_finite() && value >= 0.0)
            {
                return None;
            }
            let [a, b, c, d, e] = window.map(f64::from);
            let difference = (a + e) - 4.0 * (b + d) + 6.0 * c;
            squared_difference += difference * difference;
            weighted_intensity += (a + e) + 16.0 * (b + d) + 36.0 * c;
        }
        let gain = squared_difference / weighted_intensity;
        (gain.is_finite() && gain >= 0.0).then_some(gain)
    }
}

impl Default for CalibrantGainScratch {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ion(index: u8) -> IonAnnot {
        IonAnnot::try_new('y', Some(index), 1, 0).unwrap()
    }

    fn observed(traces: &[&[f32]]) -> MzMajorIntensityArray<IonAnnot, f32> {
        let keys = (1..=traces.len())
            .map(|i| (ion(u8::try_from(i).unwrap()), 300.0 + i as f64))
            .collect();
        let mut observed = MzMajorIntensityArray::try_new_empty(keys, traces[0].len(), 0).unwrap();
        for (index, trace) in traces.iter().enumerate() {
            observed.arr.try_replace_row_with(index, trace).unwrap();
        }
        observed
    }

    fn assert_missing(score: PoissonEvidence) {
        assert!(score.ms2_poisson_llr.is_nan());
        assert!(score.ms2_poisson_signal_to_background.is_nan());
        assert!(score.ms2_poisson_missing_fragments.is_nan());
    }

    #[test]
    fn one_cycle_uses_off_counts_and_gain() {
        let traces = observed(&[&[4.0, 90.0, 90.0, 14.0, 90.0, 90.0, 4.0]]);
        let score = compute_poisson_evidence(
            &[0.0, 0.0, 0.0, 7.0, 0.0, 0.0, 0.0],
            3,
            (1, 1),
            &traces,
            [(ion(1), 12.0)].into_iter(),
            2.0,
        );

        let expected = 2.0 * (7.0 * (7.0_f64 / 2.0).ln() - 5.0);
        assert!((score.ms2_poisson_llr - expected).abs() < 1e-12);
        assert_eq!(score.ms2_poisson_signal_to_background, 5.0 / 2.0);
        assert_eq!(score.ms2_poisson_missing_fragments, 0.0);
    }

    #[test]
    fn zero_off_counts_and_absent_predicted_row() {
        let traces = observed(&[&[0.0, 90.0, 90.0, 8.0, 8.0, 90.0, 90.0, 0.0]]);
        let score = compute_poisson_evidence(
            &[0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0],
            3,
            (1, 2),
            &traces,
            [(ion(1), 1.0), (ion(2), 1.0)].into_iter(),
            1.0,
        );

        // b is 0.25 per row, S is 15, and the absent row expects 8 counts.
        let expected = 2.0 * (16.0 * 16.0_f64.ln() - 15.0);
        assert!((score.ms2_poisson_llr - expected).abs() < 1e-12);
        assert_eq!(score.ms2_poisson_signal_to_background, 30.0);
        assert_eq!(score.ms2_poisson_missing_fragments, 1.0);
    }

    #[test]
    fn no_off_cycles_returns_all_missing() {
        assert_missing(compute_poisson_evidence(
            &[1.0, 1.0],
            0,
            (1, 2),
            &observed(&[&[8.0, 8.0]]),
            [(ion(1), 1.0), (ion(2), 1.0)].into_iter(),
            1.0,
        ));
    }

    #[test]
    fn one_newton_step_uses_heterogeneous_background() {
        let traces = observed(&[&[1.0, 0.0, 0.0, 4.0], &[9.0, 0.0, 0.0, 0.0]]);
        let score = compute_poisson_evidence(
            &[0.0, 0.0, 0.0, 1.0],
            3,
            (1, 1),
            &traces,
            [(ion(1), 1.0), (ion(2), 1.0)].into_iter(),
            1.0,
        );

        // One step gives S=1. The exact optimum S=2 gives D=1.54518.
        let expected = 2.0 * (4.0 * 1.5_f64.ln() - 1.0);
        assert!((score.ms2_poisson_llr - expected).abs() < 1e-12);
        assert_eq!(score.ms2_poisson_signal_to_background, 1.0 / 10.0);
        assert_eq!(score.ms2_poisson_missing_fragments, 1.0);
    }

    #[test]
    fn gain_pools_five_cycle_contrasts() {
        let traces = observed(&[&[3.0, 3.0, 5.0, 3.0, 3.0], &[2.0; 5]]);
        let mut scratch = CalibrantGainScratch::new();
        assert_eq!(scratch.measure(&traces, 2), Some(144.0 / 422.0));
        assert_eq!(scratch.measure(&traces, 1), None);
        assert_eq!(scratch.measure(&traces, 3), None);
    }

    #[test]
    fn skips_nonfinite_predictions_and_clips_negative_intensity() {
        let traces = observed(&[
            &[-4.0, 90.0, 90.0, 8.0, -4.0, 90.0, 90.0, -8.0],
            &[90.0; 8],
            &[90.0; 8],
        ]);
        let score = compute_poisson_evidence(
            &[0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0],
            3,
            (1, 2),
            &traces,
            [(ion(2), f32::NAN), (ion(1), 1.0), (ion(3), f32::INFINITY)].into_iter(),
            1.0,
        );

        // The usable row has counts (8, 0), b=0.25, and equal time weights.
        // The seed is S=7.5 and is already optimal: each cell expects 4 counts.
        let expected = 2.0 * (8.0 * 16.0_f64.ln() - 7.5);
        assert!((score.ms2_poisson_llr - expected).abs() < 1e-12);
        assert_eq!(score.ms2_poisson_signal_to_background, 30.0);
        assert_eq!(score.ms2_poisson_missing_fragments, 0.0);
    }

    #[test]
    fn top_predictions_keep_first_ties_through_fragment_limit() {
        let top = top_predicted((1..=MAX_FRAGMENTS).map(|i| {
            let index = u8::try_from(i).unwrap();
            let intensity = if i == MAX_FRAGMENTS { 2.0 } else { 1.0 };
            (ion(index), intensity)
        }));
        assert_eq!(top[0], (ion(u8::MAX), 2.0));
        for (index, &(key, intensity)) in top[1..].iter().enumerate() {
            assert_eq!(key, ion(u8::try_from(index + 1).unwrap()));
            assert_eq!(intensity, 1.0);
        }
    }
}
