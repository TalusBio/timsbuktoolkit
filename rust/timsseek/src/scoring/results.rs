use serde::Serialize;

use super::apex_finding::{
    ApexScore,
    PeptideMetadata,
    RelativeIntensities,
};
use super::offsets::MzMobilityOffsets;
use super::pipeline::SecondaryLazyScores;
use crate::errors::DataProcessingError;

use super::{NUM_MS2_IONS, NUM_MS1_IONS};

/// Shared scoring fields produced by Phase 3. Every field is guaranteed populated.
#[derive(Debug, Clone, Serialize)]
pub struct ScoringFields {
    // Identity
    pub sequence: String,
    pub library_id: u32,
    pub decoy_group_id: u32,
    pub precursor_mz: f64,
    pub precursor_charge: u8,
    pub precursor_mobility: f32,
    pub is_target: bool,

    // RT
    pub library_rt: f32,
    pub calibrated_rt_seconds: f32,
    pub obs_rt_seconds: f32,
    pub calibrated_sq_delta_rt: f32,

    // Mobility
    pub obs_mobility: f32,
    pub delta_ms1_ms2_mobility: f32,
    pub sq_delta_ms1_ms2_mobility: f32,

    // Primary scores
    pub main_score: f32,
    pub delta_next: f32,
    pub delta_second_next: f32,

    // Lazyscores
    pub apex_lazyscore: f32,
    pub ms2_lazyscore: f32,
    pub ms2_isotope_lazyscore: f32,
    pub ms2_isotope_lazyscore_ratio: f32,
    pub lazyscore_z: f32,
    pub lazyscore_vs_baseline: f32,

    // Split product (9 components)
    pub split_product_score: f32,
    pub cosine_au: f32,
    pub scribe_au: f32,
    pub cosine_cg: f32,
    pub scribe_cg: f32,
    pub cosine_weighted_coelution: f32,
    pub cosine_gradient_consistency: f32,
    pub scribe_weighted_coelution: f32,
    pub scribe_gradient_consistency: f32,

    // 11 apex features
    pub peak_shape: f32,
    pub ratio_cv: f32,
    pub centered_apex: f32,
    pub precursor_coelution: f32,
    pub fragment_coverage: f32,
    pub precursor_apex_match: f32,
    pub xic_quality: f32,
    pub fragment_apex_agreement: f32,
    pub isotope_correlation: f32,
    pub gaussian_correlation: f32,
    pub per_frag_gaussian_corr: f32,

    // Peak shape
    pub rising_cycles: u8,
    pub falling_cycles: u8,

    // Counts
    pub npeaks: u8,
    pub n_scored_fragments: u8,

    // Intensities
    pub ms2_summed_intensity: f32,
    pub ms1_summed_intensity: f32,

    // Per-ion errors (real arrays, not numbered fields)
    pub ms2_mz_errors: [f32; NUM_MS2_IONS],
    pub ms2_mobility_errors: [f32; NUM_MS2_IONS],
    pub ms1_mz_errors: [f32; NUM_MS1_IONS],
    pub ms1_mobility_errors: [f32; NUM_MS1_IONS],

    // Relative intensities
    pub ms2_intensity_ratios: [f32; NUM_MS2_IONS],
    pub ms1_intensity_ratios: [f32; NUM_MS1_IONS],
}

/// Phase 3 output. All scoring fields guaranteed populated.
#[derive(Debug, Clone, Serialize)]
pub struct ScoredCandidate {
    pub scoring: ScoringFields,
}

/// After target-decoy competition.
#[derive(Debug, Clone, Serialize)]
pub struct CompetedCandidate {
    pub scoring: ScoringFields,
    pub delta_group: f32,
    pub delta_group_ratio: f32,
    /// Scratch field for CrossValidatedScorer (written during rescore)
    pub(crate) discriminant_score: f32,
    /// Scratch field for q-value assignment
    pub(crate) qvalue: f32,
}

/// After rescoring. Written to Parquet.
#[derive(Debug, Clone, Serialize)]
pub struct FinalResult {
    pub scoring: ScoringFields,
    pub delta_group: f32,
    pub delta_group_ratio: f32,
    pub discriminant_score: f32,
    pub qvalue: f32,
}

// ---------------------------------------------------------------------------
// Stage conversions
// ---------------------------------------------------------------------------

impl ScoredCandidate {
    /// Convert into a `CompetedCandidate` with the given delta-group values.
    ///
    /// Items that are alone in their group (no competitor) should pass
    /// `f32::NAN` for both deltas.
    pub fn into_competed(self, delta_group: f32, delta_group_ratio: f32) -> CompetedCandidate {
        CompetedCandidate {
            scoring: self.scoring,
            delta_group,
            delta_group_ratio,
            discriminant_score: f32::NAN,
            qvalue: f32::NAN,
        }
    }
}

impl CompetedCandidate {
    /// Promote to a `FinalResult` (all fields frozen).
    pub fn into_final(self) -> FinalResult {
        FinalResult {
            scoring: self.scoring,
            delta_group: self.delta_group,
            delta_group_ratio: self.delta_group_ratio,
            discriminant_score: self.discriminant_score,
            qvalue: self.qvalue,
        }
    }
}

// ---------------------------------------------------------------------------
// Builder
// ---------------------------------------------------------------------------

/// Tracks whether a field has been set or is still unset.
///
/// Mirrors the `SetField` in `search_results.rs` but lives here so
/// `ScoredCandidateBuilder` is self-contained.
#[derive(Debug, Clone, Copy, Default)]
pub enum SetField<T> {
    Some(T),
    #[default]
    None,
}

impl<T> SetField<T> {
    pub fn is_some(&self) -> bool {
        matches!(self, Self::Some(_))
    }

    pub fn expect_some(self, field_name: &'static str) -> Result<T, DataProcessingError> {
        match self {
            Self::Some(v) => Ok(v),
            Self::None => Err(DataProcessingError::ExpectedSetField {
                field: field_name,
                context: "".into(),
            }),
        }
    }
}

/// Builder that collects all inputs for a `ScoredCandidate` and validates
/// completeness in `finalize()`.
#[derive(Debug, Default)]
pub struct ScoredCandidateBuilder {
    // --- Identity ---
    sequence: SetField<String>,
    library_id: SetField<u32>,
    decoy_group_id: SetField<u32>,
    precursor_mz: SetField<f64>,
    precursor_charge: SetField<u8>,
    precursor_mobility: SetField<f32>,
    is_target: SetField<bool>,

    // --- RT ---
    library_rt: SetField<f32>,
    calibrated_rt_seconds: SetField<f32>,

    // --- Observed RT / mobility ---
    obs_rt_seconds: SetField<f32>,
    obs_mobility: SetField<f32>,
    delta_ms1_ms2_mobility: SetField<f32>,

    // --- Primary scores ---
    main_score: SetField<f32>,
    delta_next: SetField<f32>,
    delta_second_next: SetField<f32>,

    // --- Lazyscores ---
    apex_lazyscore: SetField<f32>,
    ms2_lazyscore: SetField<f32>,
    ms2_isotope_lazyscore: SetField<f32>,
    ms2_isotope_lazyscore_ratio: SetField<f32>,
    lazyscore_z: SetField<f32>,
    lazyscore_vs_baseline: SetField<f32>,

    // --- Split product ---
    split_product_score: SetField<f32>,
    cosine_au: SetField<f32>,
    scribe_au: SetField<f32>,
    cosine_cg: SetField<f32>,
    scribe_cg: SetField<f32>,
    cosine_weighted_coelution: SetField<f32>,
    cosine_gradient_consistency: SetField<f32>,
    scribe_weighted_coelution: SetField<f32>,
    scribe_gradient_consistency: SetField<f32>,

    // --- 11 apex features ---
    peak_shape: SetField<f32>,
    ratio_cv: SetField<f32>,
    centered_apex: SetField<f32>,
    precursor_coelution: SetField<f32>,
    fragment_coverage: SetField<f32>,
    precursor_apex_match: SetField<f32>,
    xic_quality: SetField<f32>,
    fragment_apex_agreement: SetField<f32>,
    isotope_correlation: SetField<f32>,
    gaussian_correlation: SetField<f32>,
    per_frag_gaussian_corr: SetField<f32>,

    // --- Peak shape ---
    rising_cycles: SetField<u8>,
    falling_cycles: SetField<u8>,

    // --- Counts ---
    npeaks: SetField<u8>,
    n_scored_fragments: SetField<u8>,

    // --- Intensities ---
    ms2_summed_intensity: SetField<f32>,
    ms1_summed_intensity: SetField<f32>,

    // --- Per-ion errors ---
    ms2_mz_errors: SetField<[f32; NUM_MS2_IONS]>,
    ms2_mobility_errors: SetField<[f32; NUM_MS2_IONS]>,
    ms1_mz_errors: SetField<[f32; NUM_MS1_IONS]>,
    ms1_mobility_errors: SetField<[f32; NUM_MS1_IONS]>,

    // --- Relative intensities ---
    relative_intensities: SetField<RelativeIntensities>,
}

impl ScoredCandidateBuilder {
    /// Populate identity fields and reference values from peptide metadata.
    pub fn with_metadata(mut self, metadata: &PeptideMetadata) -> Self {
        self.library_id = SetField::Some(metadata.library_id);
        self.sequence = SetField::Some(String::from(metadata.digest.clone()));
        self.is_target = SetField::Some(metadata.digest.decoy.is_target());
        self.decoy_group_id = SetField::Some(metadata.digest.decoy_group);
        self.precursor_charge = SetField::Some(metadata.charge);
        self.precursor_mz = SetField::Some(metadata.ref_precursor_mz);
        self.precursor_mobility = SetField::Some(metadata.ref_mobility_ook0);
        self.library_rt = SetField::Some(metadata.library_rt);
        self.calibrated_rt_seconds = SetField::Some(metadata.calibrated_rt_seconds);
        self
    }

    /// Set the number of scored fragments (ions used during scoring).
    pub fn with_nqueries(mut self, nqueries: u8) -> Self {
        self.n_scored_fragments = SetField::Some(nqueries);
        self
    }

    /// Populate per-ion m/z and mobility error arrays plus observed mobility.
    pub fn with_sorted_offsets(mut self, offsets: &MzMobilityOffsets) -> Self {
        self.ms1_mz_errors = SetField::Some(offsets.ms1_mz_errors());
        self.ms1_mobility_errors = SetField::Some(offsets.ms1_mobility_errors());
        self.ms2_mz_errors = SetField::Some(offsets.ms2_mz_errors());
        self.ms2_mobility_errors = SetField::Some(offsets.ms2_mobility_errors());

        let mob_errors = offsets.avg_delta_mobs();
        let cum_err = mob_errors.0 + mob_errors.1;
        let obs_mob =
            (offsets.ref_mobility + cum_err.mean_mobility().unwrap_or(f64::NAN)) as f32;
        let d_err = match (mob_errors.0.mean_mobility(), mob_errors.1.mean_mobility()) {
            (Ok(mz), Ok(mob)) => mz - mob,
            _ => f64::NAN,
        };
        self.delta_ms1_ms2_mobility = SetField::Some(d_err as f32);
        self.obs_mobility = SetField::Some(obs_mob);
        self
    }

    /// Populate secondary lazyscore fields from the isotope/inner collectors.
    pub fn with_secondary_lazyscores(mut self, lazyscores: SecondaryLazyScores) -> Self {
        self.ms2_lazyscore = SetField::Some(lazyscores.lazyscore);
        self.ms2_isotope_lazyscore = SetField::Some(lazyscores.iso_lazyscore);
        self.ms2_isotope_lazyscore_ratio = SetField::Some(lazyscores.ratio);
        self
    }

    /// Populate relative intensity arrays from the inner collector.
    pub fn with_relative_intensities(mut self, relative_intensities: RelativeIntensities) -> Self {
        self.relative_intensities = SetField::Some(relative_intensities);
        self
    }

    /// Populate all fields derived from the full apex score.
    pub fn with_apex_score(mut self, main_score: &ApexScore) -> Self {
        let ApexScore {
            score,
            retention_time_ms,
            joint_apex_cycle: _,
            split_product,
            features,
            delta_next,
            delta_second_next,
            lazyscore,
            lazyscore_vs_baseline,
            lazyscore_z,
            npeaks,
            ms2_summed_intensity,
            ms1_summed_intensity,
            rising_cycles,
            falling_cycles,
        } = *main_score;

        self.main_score = SetField::Some(score);
        self.delta_next = SetField::Some(delta_next);
        self.delta_second_next = SetField::Some(delta_second_next);
        self.obs_rt_seconds = SetField::Some(retention_time_ms as f32 / 1000.0);

        self.split_product_score = SetField::Some(split_product.base_score);
        self.cosine_au = SetField::Some(split_product.cosine_au);
        self.scribe_au = SetField::Some(split_product.scribe_au);
        self.cosine_cg = SetField::Some(split_product.cosine_cg);
        self.scribe_cg = SetField::Some(split_product.scribe_cg);
        self.cosine_weighted_coelution =
            SetField::Some(split_product.cosine_weighted_coelution);
        self.cosine_gradient_consistency =
            SetField::Some(split_product.cosine_gradient_consistency);
        self.scribe_weighted_coelution =
            SetField::Some(split_product.scribe_weighted_coelution);
        self.scribe_gradient_consistency =
            SetField::Some(split_product.scribe_gradient_consistency);

        self.peak_shape = SetField::Some(features.peak_shape);
        self.ratio_cv = SetField::Some(features.ratio_cv);
        self.centered_apex = SetField::Some(features.centered_apex);
        self.precursor_coelution = SetField::Some(features.precursor_coelution);
        self.fragment_coverage = SetField::Some(features.fragment_coverage);
        self.precursor_apex_match = SetField::Some(features.precursor_apex_match);
        self.xic_quality = SetField::Some(features.xic_quality);
        self.fragment_apex_agreement = SetField::Some(features.fragment_apex_agreement);
        self.isotope_correlation = SetField::Some(features.isotope_correlation);
        self.gaussian_correlation = SetField::Some(features.gaussian_correlation);
        self.per_frag_gaussian_corr = SetField::Some(features.per_frag_gaussian_corr);

        self.apex_lazyscore = SetField::Some(lazyscore);
        self.lazyscore_z = SetField::Some(lazyscore_z);
        self.lazyscore_vs_baseline = SetField::Some(lazyscore_vs_baseline);
        self.npeaks = SetField::Some(npeaks);
        self.ms1_summed_intensity = SetField::Some(ms1_summed_intensity);
        self.ms2_summed_intensity = SetField::Some(ms2_summed_intensity);
        self.rising_cycles = SetField::Some(rising_cycles);
        self.falling_cycles = SetField::Some(falling_cycles);

        self
    }

    /// Validate completeness and construct a `ScoredCandidate`.
    pub fn finalize(self) -> Result<ScoredCandidate, DataProcessingError> {
        macro_rules! expect_some {
            ($field:ident) => {
                self.$field.expect_some(stringify!($field))?
            };
        }

        let obs_rt_seconds = expect_some!(obs_rt_seconds);
        let library_rt = expect_some!(library_rt);
        let calibrated_rt_seconds = expect_some!(calibrated_rt_seconds);
        let calibrated_delta_rt = obs_rt_seconds - calibrated_rt_seconds;
        let calibrated_sq_delta_rt = calibrated_delta_rt * calibrated_delta_rt;

        let delta_ms1_ms2_mobility = expect_some!(delta_ms1_ms2_mobility);
        let sq_delta_ms1_ms2_mobility = delta_ms1_ms2_mobility * delta_ms1_ms2_mobility;

        let relints = expect_some!(relative_intensities);
        let ms1_intensity_ratios = relints.ms1.get_values();
        let ms2_intensity_ratios = relints.ms2.get_values();

        let scoring = ScoringFields {
            // Identity
            sequence: expect_some!(sequence),
            library_id: expect_some!(library_id),
            decoy_group_id: expect_some!(decoy_group_id),
            precursor_mz: expect_some!(precursor_mz),
            precursor_charge: expect_some!(precursor_charge),
            precursor_mobility: expect_some!(precursor_mobility),
            is_target: expect_some!(is_target),

            // RT
            library_rt,
            calibrated_rt_seconds,
            obs_rt_seconds,
            calibrated_sq_delta_rt,

            // Mobility
            obs_mobility: expect_some!(obs_mobility),
            delta_ms1_ms2_mobility,
            sq_delta_ms1_ms2_mobility,

            // Primary scores
            main_score: expect_some!(main_score),
            delta_next: expect_some!(delta_next),
            delta_second_next: expect_some!(delta_second_next),

            // Lazyscores
            apex_lazyscore: expect_some!(apex_lazyscore),
            ms2_lazyscore: expect_some!(ms2_lazyscore),
            ms2_isotope_lazyscore: expect_some!(ms2_isotope_lazyscore),
            ms2_isotope_lazyscore_ratio: expect_some!(ms2_isotope_lazyscore_ratio),
            lazyscore_z: expect_some!(lazyscore_z),
            lazyscore_vs_baseline: expect_some!(lazyscore_vs_baseline),

            // Split product
            split_product_score: expect_some!(split_product_score),
            cosine_au: expect_some!(cosine_au),
            scribe_au: expect_some!(scribe_au),
            cosine_cg: expect_some!(cosine_cg),
            scribe_cg: expect_some!(scribe_cg),
            cosine_weighted_coelution: expect_some!(cosine_weighted_coelution),
            cosine_gradient_consistency: expect_some!(cosine_gradient_consistency),
            scribe_weighted_coelution: expect_some!(scribe_weighted_coelution),
            scribe_gradient_consistency: expect_some!(scribe_gradient_consistency),

            // 11 apex features
            peak_shape: expect_some!(peak_shape),
            ratio_cv: expect_some!(ratio_cv),
            centered_apex: expect_some!(centered_apex),
            precursor_coelution: expect_some!(precursor_coelution),
            fragment_coverage: expect_some!(fragment_coverage),
            precursor_apex_match: expect_some!(precursor_apex_match),
            xic_quality: expect_some!(xic_quality),
            fragment_apex_agreement: expect_some!(fragment_apex_agreement),
            isotope_correlation: expect_some!(isotope_correlation),
            gaussian_correlation: expect_some!(gaussian_correlation),
            per_frag_gaussian_corr: expect_some!(per_frag_gaussian_corr),

            // Peak shape
            rising_cycles: expect_some!(rising_cycles),
            falling_cycles: expect_some!(falling_cycles),

            // Counts
            npeaks: expect_some!(npeaks),
            n_scored_fragments: expect_some!(n_scored_fragments),

            // Intensities
            ms2_summed_intensity: expect_some!(ms2_summed_intensity),
            ms1_summed_intensity: expect_some!(ms1_summed_intensity),

            // Per-ion errors
            ms2_mz_errors: expect_some!(ms2_mz_errors),
            ms2_mobility_errors: expect_some!(ms2_mobility_errors),
            ms1_mz_errors: expect_some!(ms1_mz_errors),
            ms1_mobility_errors: expect_some!(ms1_mobility_errors),

            // Relative intensities
            ms2_intensity_ratios,
            ms1_intensity_ratios,
        };

        Ok(ScoredCandidate { scoring })
    }
}
