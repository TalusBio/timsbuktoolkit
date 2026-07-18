//! Synthetic extraction-product generator.
//!
//! Builds a [`timsseek::scoring::apex_finding::Extraction`] from tunable knobs
//! WITHOUT any real timsTOF index or `.d` data. The produced extraction is fed
//! verbatim to the real `TraceScorer`, so what the app scores is exactly what
//! the timsseek pipeline would score.
//!
//! Model: all "real" fragments co-elute at a single shared apex cycle (no RT
//! drift between real fragments -- coelution is always real). Realism knobs:
//! per-fragment noise, theoretical-vs-observed ratio mismatch (`obs_scale`),
//! absent/empty transitions (`obs_scale == 0`), randomly injected peak-like
//! objects, and a global noise floor.

use rand::{
    Rng,
    SeedableRng,
};
use rand_chacha::ChaCha8Rng;

use timsquery::{
    ChromatogramCollector,
    MzMajorIntensityArray,
    TupleRange,
};
use timsseek::ExpectedIntensities;
use timsseek::scoring::apex_finding::Extraction;
use timsseek::scoring::pipeline::{
    TOP_N_FRAGMENTS,
    filter_zero_intensity_ions,
    select_top_n_fragments,
};

/// A library fragment. Theoretical and observed intensities are decoupled:
/// `theo_intensity` is what the spectral library predicts (goes into
/// `ExpectedIntensities`), while the observed peak height is
/// `theo_intensity * obs_scale`. Real speclibs are often measured on a
/// different instrument, so observed ratios rarely match theory exactly.
///
/// - `obs_scale == 1.0` -> observed matches theory.
/// - `obs_scale != 1.0` -> instrument ratio mismatch.
/// - `obs_scale == 0.0` -> transition is EXPECTED but ABSENT: its trace is
///   pure noise (an "empty" transition), which stresses fragment coverage.
#[derive(Debug, Clone)]
pub struct FragmentSpec {
    pub label: String,
    /// Theoretical (library) relative intensity -> `ExpectedIntensities`.
    pub theo_intensity: f32,
    /// Observed height multiplier relative to theory (0 => empty/noise-only).
    pub obs_scale: f32,
    /// Per-fragment noise multiplier (some fragments noisier than others).
    pub noise_mult: f32,
}

/// Random "peak-like objects" injected onto existing transitions/precursors.
///
/// These are `count`
/// gaussian bumps scattered uniformly at random -- each lands on a randomly
/// chosen existing transition (and optionally precursor) row, at a uniform-
/// random cycle, with height/width drawn uniformly from the given ranges.
/// Because they add onto REAL rows, they perturb cosine/scribe/apex directly.
#[derive(Debug, Clone)]
pub struct RandomPeaks {
    pub enabled: bool,
    pub count: usize,
    /// Peak height as a fraction of the global `height`, drawn uniformly.
    pub height_frac_min: f32,
    pub height_frac_max: f32,
    /// If true, precursor rows are also eligible targets.
    pub hit_precursors: bool,
}

impl Default for RandomPeaks {
    fn default() -> Self {
        Self {
            enabled: true,
            count: 200,
            height_frac_min: 0.1,
            height_frac_max: 10.0,
            hit_precursors: true,
        }
    }
}

/// Full set of simulation knobs. Deterministic given `seed`.
#[derive(Debug, Clone)]
pub struct SimParams {
    pub n_cycles: usize,
    /// Shared apex position (cycle index) for all real fragments + precursor.
    pub apex_cycle: f32,
    /// Elution peak width (gaussian sigma, in cycles).
    pub width_sigma: f32,
    /// Global height scale applied on top of each fragment's rel_intensity.
    pub height: f32,

    pub real_fragments: Vec<FragmentSpec>,
    pub random_peaks: RandomPeaks,

    /// Number of precursor isotopes (keys 0..n, all positive => scored).
    pub n_isotopes: usize,
    pub precursor_intensity: f32,

    /// Global additive noise scale (uniform 0..1 * this, per cell).
    pub noise_floor: f32,

    pub seed: u64,

    // Cycle <-> ms mapping for the rt_mapper.
    pub rt_start_ms: u32,
    pub cycle_period_ms: u32,
}

impl Default for SimParams {
    fn default() -> Self {
        // (label, theo_intensity, obs_scale, noise_mult)
        // obs_scale defaults deviate from 1.0 => realistic cross-instrument
        // ratio mismatch out of the box (observed != theoretical).
        let real_fragments = vec![
            ("y3", 1.00, 0.9, 1.0),
            ("y4", 0.80, 1.1, 1.0),
            ("y5", 0.55, 0.6, 1.5),
            ("b3", 0.40, 1.3, 1.0),
            ("b4", 0.30, 0.8, 2.5),
            ("y6", 0.20, 1.0, 1.0),
        ]
        .into_iter()
        .map(|(l, theo, obs, n)| FragmentSpec {
            label: l.to_string(),
            theo_intensity: theo,
            obs_scale: obs,
            noise_mult: n,
        })
        .collect();

        Self {
            n_cycles: 60,
            apex_cycle: 30.0,
            width_sigma: 1.0,
            height: 1000.0,
            real_fragments,
            random_peaks: RandomPeaks::default(),
            n_isotopes: 3,
            precursor_intensity: 0.7,
            noise_floor: 0.02,
            seed: 42,
            rt_start_ms: 60_000,
            cycle_period_ms: 1_000,
        }
    }
}

impl SimParams {
    /// Linear cycle -> retention-time (ms) mapper, matching what the pipeline
    /// passes into `TraceScorer::suggest_apex` / `score_at`.
    pub fn rt_mapper(&self) -> impl Fn(usize) -> u32 + '_ {
        move |cycle: usize| self.rt_start_ms + (cycle as u32) * self.cycle_period_ms
    }
}

/// A single observed transition row, kept for display alongside the extraction.
#[derive(Debug, Clone)]
pub struct TransitionRow {
    pub label: String,
    /// True when the transition is expected but observed as noise only
    /// (`obs_scale == 0`). Drawn dashed so it reads as "missing".
    pub is_absent: bool,
    pub intensities: Vec<f32>,
}

/// The synthesized data: the `Extraction` fed to the scorer, plus the raw rows
/// (fragments + precursors) kept separately so the UI can plot them directly.
pub struct SimData {
    pub extraction: Extraction<String>,
    pub fragment_rows: Vec<TransitionRow>,
    pub precursor_rows: Vec<TransitionRow>,
}

/// Gaussian elution value at `cycle` for a peak centered at `center`.
fn gaussian(cycle: f32, center: f32, sigma: f32, height: f32) -> f32 {
    let z = (cycle - center) / sigma;
    height * (-0.5 * z * z).exp()
}

/// Build the full synthetic extraction + display rows from `params`.
pub fn build(params: &SimParams) -> SimData {
    let mut rng = ChaCha8Rng::seed_from_u64(params.seed);
    let n = params.n_cycles;
    let dummy_mz = 500.0_f64;

    // --- Expected intensities: THEORETICAL (library) ratios. ---
    // These drive cosine/scribe. Observed peaks below may deviate (obs_scale).
    let expected = ExpectedIntensities::try_from_pairs(
        params
            .real_fragments
            .iter()
            .map(|f| (f.label.clone(), f.theo_intensity)),
        (0..params.n_isotopes as i8).map(|iso| {
            // simple decaying isotope pattern
            (iso, params.precursor_intensity * 0.6_f32.powi(iso as i32))
        }),
    )
    .expect("labels are unique by construction");

    // --- Fragment rows: OBSERVED = theo * obs_scale (co-eluting). ---
    let mut frag_order: Vec<(String, f64)> = Vec::new();
    let mut frag_rows: Vec<TransitionRow> = Vec::new();

    for f in &params.real_fragments {
        // Observed height decoupled from theory. obs_scale == 0 => absent
        // transition: no peak, trace is pure noise.
        let peak = params.height * f.theo_intensity * f.obs_scale;
        let noise = params.noise_floor * params.height * f.noise_mult;
        let intensities = (0..n)
            .map(|c| {
                let signal = gaussian(c as f32, params.apex_cycle, params.width_sigma, peak);
                sample_cell(&mut rng, signal, noise)
            })
            .collect();
        frag_order.push((f.label.clone(), dummy_mz));
        frag_rows.push(TransitionRow {
            label: f.label.clone(),
            is_absent: f.obs_scale <= 0.0,
            intensities,
        });
    }

    // --- Precursor rows (isotopes, co-eluting at shared apex) ---
    let mut prec_order: Vec<(i8, f64)> = Vec::new();
    let mut prec_rows: Vec<TransitionRow> = Vec::new();
    for iso in 0..params.n_isotopes as i8 {
        let peak = params.height * params.precursor_intensity * 0.6_f32.powi(iso as i32);
        let noise = params.noise_floor * params.height;
        let intensities = (0..n)
            .map(|c| {
                let signal = gaussian(c as f32, params.apex_cycle, params.width_sigma, peak);
                sample_cell(&mut rng, signal, noise)
            })
            .collect();
        prec_order.push((iso, dummy_mz));
        prec_rows.push(TransitionRow {
            label: format!("M+{iso}"),
            is_absent: false,
            intensities,
        });
    }

    // --- Inject random peak-like objects onto existing rows ---
    inject_random_peaks(&mut rng, &mut frag_rows, &mut prec_rows, params);

    // --- Pack into MzMajorIntensityArrays (cycle_offset = 0) ---
    let mut fragments = MzMajorIntensityArray::<String, f32>::try_new_empty(frag_order, n, 0)
        .expect("non-empty fragments");
    fill_array(&mut fragments, &frag_rows);

    let mut precursors = MzMajorIntensityArray::<i8, f32>::try_new_empty(prec_order, n, 0)
        .expect("non-empty precursors");
    fill_array(&mut precursors, &prec_rows);

    let n_frag_peaks: u64 = frag_rows
        .iter()
        .map(|r| r.intensities.iter().filter(|v| **v > 0.0).count() as u64)
        .sum();
    let n_prec_peaks: u64 = prec_rows
        .iter()
        .map(|r| r.intensities.iter().filter(|v| **v > 0.0).count() as u64)
        .sum();

    let map = params.rt_mapper();
    let rt_range_ms =
        TupleRange::try_new(map(0), map(n - 1)).expect("start < end for positive period");

    let chromatograms = ChromatogramCollector::<String, f32> {
        id: 0,
        mobility_ook0: 1.0,
        rt_seconds: (map((params.apex_cycle as usize).min(n - 1)) as f32) / 1000.0,
        precursor_mono_mz: dummy_mz,
        precursor_charge: 2,
        precursor_mz_limits: (dummy_mz - 1.0, dummy_mz + 1.0),
        precursors,
        fragments,
        rt_range_ms,
        n_precursor_peaks_added: n_prec_peaks,
        n_fragment_peaks_added: n_frag_peaks,
        n_quad_windows_matched: 1,
    };

    let mut extraction = Extraction {
        expected_intensities: expected,
        chromatograms,
    };

    // Apply the SAME extraction-time filtering the production pipeline always
    // runs before scoring (broad + calibrated paths, extraction.rs): drop
    // zero-intensity ions, then keep the top-N fragments by predicted
    // intensity. Without this the scorer would see rows production discards
    // (e.g. an absent transition's all-noise row, or >TOP_N_FRAGMENTS ions),
    // making the sim diverge from what timsseek actually scores.
    filter_zero_intensity_ions(
        &mut extraction.chromatograms,
        &mut extraction.expected_intensities,
    );
    select_top_n_fragments(
        &mut extraction.chromatograms,
        &mut extraction.expected_intensities,
        TOP_N_FRAGMENTS,
    );

    // Display rows are kept RAW (everything simulated), so the plots show the
    // full input; the scored `extraction` reflects post-filter reality.
    SimData {
        extraction,
        fragment_rows: frag_rows,
        precursor_rows: prec_rows,
    }
}

/// Scatter `count` random gaussian bumps across the transition (and optionally
/// precursor) rows. Each bump: uniform-random target row, uniform-random center
/// cycle, uniform height (fraction of global height) and sigma.
fn inject_random_peaks(
    rng: &mut ChaCha8Rng,
    frag_rows: &mut [TransitionRow],
    prec_rows: &mut [TransitionRow],
    params: &SimParams,
) {
    let rp = &params.random_peaks;
    if !rp.enabled || rp.count == 0 {
        return;
    }
    let n = params.n_cycles;
    let n_frag = frag_rows.len();
    let n_targets = n_frag
        + if rp.hit_precursors {
            prec_rows.len()
        } else {
            0
        };
    if n_targets == 0 {
        return;
    }

    for _ in 0..rp.count {
        let target = rng.random_range(0..n_targets);
        let row = if target < n_frag {
            &mut frag_rows[target]
        } else {
            &mut prec_rows[target - n_frag]
        };
        let center = rng.random_range(0.0..n as f32);
        // Width inherited from the peptide's elution width (the fragments the
        // injected peak is contaminating).
        let sigma = params.width_sigma;
        let hfrac =
            rng.random_range(rp.height_frac_min..=rp.height_frac_max.max(rp.height_frac_min));
        let peak = params.height * hfrac;
        for (cy, v) in row.intensities.iter_mut().enumerate() {
            *v += gaussian(cy as f32, center, sigma, peak);
        }
    }
}

/// Sample one observed cell: signal + uniform noise floor, clamped >=0.
fn sample_cell(rng: &mut ChaCha8Rng, signal: f32, noise: f32) -> f32 {
    (signal + rng.random::<f32>() * noise).max(0.0)
}

/// Copy display rows into an [`MzMajorIntensityArray`] via `add_at_index`.
fn fill_array<K: timsquery::KeyLike>(
    arr: &mut MzMajorIntensityArray<K, f32>,
    rows: &[TransitionRow],
) {
    // iter_mut_mzs yields rows in mz_order order == `rows` order; zip directly.
    for (row, ((_key, _mz), mut chrom)) in rows.iter().zip(arr.iter_mut_mzs()) {
        for (cycle, &val) in row.intensities.iter().enumerate() {
            if val > 0.0 {
                chrom.add_at_index(cycle as u32, val);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clean_peak_apex_lands_at_configured_cycle() {
        let mut p = SimParams::default();
        p.noise_floor = 0.0;
        p.random_peaks.enabled = false;
        p.apex_cycle = 25.0;
        let data = build(&p);
        // The tallest real fragment row should peak at cycle 25.
        let row = &data.fragment_rows[0];
        let (max_idx, _) = row
            .intensities
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap();
        assert_eq!(max_idx, 25);
    }

    #[test]
    fn deterministic_given_seed() {
        let p = SimParams::default();
        let a = build(&p);
        let b = build(&p);
        assert_eq!(
            a.fragment_rows[0].intensities,
            b.fragment_rows[0].intensities
        );
    }

    #[test]
    fn empty_transition_at_zero_noise_is_filtered_out() {
        // At zero noise an absent (obs_scale=0) transition is an all-zero row.
        // Production's filter_zero_intensity_ions drops such rows AND their
        // expected entry, so the sim must too (matches timsseek exactly).
        let mut p = SimParams::default();
        p.noise_floor = 0.0;
        p.random_peaks.enabled = false;
        p.real_fragments[0].obs_scale = 0.0;
        let label = p.real_fragments[0].label.clone();
        let data = build(&p);
        // Display keeps the raw simulated row (all zeros, marked absent)...
        assert!(data.fragment_rows[0].intensities.iter().all(|&v| v == 0.0));
        assert!(data.fragment_rows[0].is_absent);
        // ...but it is dropped from the SCORED extraction's expected set.
        assert!(
            data.extraction
                .expected_intensities
                .get_fragment(&label)
                .is_none()
        );
    }

    #[test]
    fn top_n_fragments_capped_for_scoring() {
        // Add well beyond TOP_N; the scored extraction caps to TOP_N while
        // display keeps all raw rows.
        let mut p = SimParams::default();
        p.random_peaks.enabled = false;
        for i in 0..12 {
            p.real_fragments.push(FragmentSpec {
                label: format!("z{i}"),
                theo_intensity: 0.05,
                obs_scale: 1.0,
                noise_mult: 1.0,
            });
        }
        let data = build(&p);
        assert!(data.fragment_rows.len() > 8);
        assert!(data.extraction.expected_intensities.fragment_len() <= 8);
    }

    #[test]
    fn observed_decoupled_from_theoretical() {
        let mut p = SimParams::default();
        p.noise_floor = 0.0;
        p.random_peaks.enabled = false;
        p.apex_cycle = 30.0;
        p.real_fragments[0].theo_intensity = 1.0;
        p.real_fragments[0].obs_scale = 0.5; // observed at half of theory
        let data = build(&p);
        let observed_peak = data.fragment_rows[0]
            .intensities
            .iter()
            .cloned()
            .fold(0.0_f32, f32::max);
        // Observed peak reflects obs_scale, not the theoretical ratio.
        assert!((observed_peak - p.height * 0.5).abs() < 1.0);
    }
}
