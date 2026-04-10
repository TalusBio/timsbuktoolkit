use serde::Serialize;

pub const NUM_MS2_IONS: usize = 7;
pub const NUM_MS1_IONS: usize = 3;

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
    pub query_rt_seconds: f32,
    pub obs_rt_seconds: f32,
    pub delta_rt: f32,
    pub sq_delta_rt: f32,
    pub calibrated_sq_delta_rt: f32,
    pub recalibrated_rt: f32,

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
