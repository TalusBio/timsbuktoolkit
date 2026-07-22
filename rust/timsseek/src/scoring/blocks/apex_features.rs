//! Apex-features family: the 11 apex-local features (METHODS.md §3.4).
//!
//! The struct is reused whole from [`crate::scoring::scores::apex_features`]
//! (it is the return type of `compute_apex_features`, which needs the live
//! chromatogram buffers). Here we give it its [`ScoreBlock`] projection by
//! hand — 11 raw columns, 11 raw features.

use crate::scoring::blocks::{
    ColSink,
    FeatSink,
    ScoreBlock,
};
pub use crate::scoring::scores::apex_features::ApexFeatures;

impl ScoreBlock for ApexFeatures {
    fn columns(&self, o: &mut ColSink) {
        o.f32("peak_shape", self.peak_shape);
        o.f32("ratio_cv", self.ratio_cv);
        o.f32("centered_apex", self.centered_apex);
        o.f32("precursor_coelution", self.precursor_coelution);
        o.f32("fragment_coverage", self.fragment_coverage);
        o.f32("precursor_apex_match", self.precursor_apex_match);
        o.f32("xic_quality", self.xic_quality);
        o.f32("fragment_apex_agreement", self.fragment_apex_agreement);
        o.f32("isotope_correlation", self.isotope_correlation);
        o.f32("gaussian_correlation", self.gaussian_correlation);
        o.f32("per_frag_gaussian_corr", self.per_frag_gaussian_corr);
    }

    fn features(&self, o: &mut FeatSink) {
        o.push("peak_shape", self.peak_shape as f64);
        o.push("ratio_cv", self.ratio_cv as f64);
        o.push("centered_apex", self.centered_apex as f64);
        o.push("precursor_coelution", self.precursor_coelution as f64);
        o.push("fragment_coverage", self.fragment_coverage as f64);
        o.push("precursor_apex_match", self.precursor_apex_match as f64);
        o.push("xic_quality", self.xic_quality as f64);
        o.push(
            "fragment_apex_agreement",
            self.fragment_apex_agreement as f64,
        );
        o.push("isotope_correlation", self.isotope_correlation as f64);
        o.push("gaussian_correlation", self.gaussian_correlation as f64);
        o.push("per_frag_gaussian_corr", self.per_frag_gaussian_corr as f64);
    }
}

impl ApexFeatures {
    pub fn sample() -> Self {
        Self {
            peak_shape: 0.5,
            ratio_cv: 0.5,
            centered_apex: 0.5,
            precursor_coelution: 0.5,
            fragment_coverage: 0.5,
            precursor_apex_match: 0.5,
            xic_quality: 0.5,
            fragment_apex_agreement: 0.5,
            isotope_correlation: 0.5,
            gaussian_correlation: 0.5,
            per_frag_gaussian_corr: 0.5,
        }
    }
}
