//! Relative-intensity arrays (per-ion intensity relative to the MS-level
//! total). Built from the secondary-query collector's sorted values.
//!
//! Same array-family shape as [`super::ion_errors`]: names generated from a
//! prefix + ion count (`f32_array` for columns, `push_indexed` for features),
//! so no name table. Column and feature order both run MS2→MS1.
//!
//! NOTE: adding/renaming a name here also requires updating the parquet schema
//! golden (`parquet_writer`).

use crate::scoring::apex_finding::RelativeIntensityCollector;
use crate::scoring::blocks::{
    ColSink,
    FeatSink,
    NameSink,
    ScoreBlock,
};
use crate::scoring::{
    NUM_MS1_IONS,
    NUM_MS2_IONS,
};

#[derive(Debug, Clone, Copy, serde::Serialize)]
pub struct RelativeIntensities {
    pub ms2_intensity_ratios: [f32; NUM_MS2_IONS],
    pub ms1_intensity_ratios: [f32; NUM_MS1_IONS],
}

impl RelativeIntensities {
    pub fn from_collector(r: &RelativeIntensityCollector) -> Self {
        Self {
            ms2_intensity_ratios: r.ms2.get_values_sorted(),
            ms1_intensity_ratios: r.ms1.get_values_sorted(),
        }
    }

    pub fn sample() -> Self {
        Self {
            ms2_intensity_ratios: [0.5; NUM_MS2_IONS],
            ms1_intensity_ratios: [0.5; NUM_MS1_IONS],
        }
    }
}

impl ScoreBlock for RelativeIntensities {
    fn columns(&self, o: &mut ColSink) {
        o.f32_array("ms2_intensity_ratio", &self.ms2_intensity_ratios);
        o.f32_array("ms1_intensity_ratio", &self.ms1_intensity_ratios);
    }

    fn features(&self, o: &mut FeatSink) {
        o.push_slice(&self.ms2_intensity_ratios);
        o.push_slice(&self.ms1_intensity_ratios);
    }

    fn feature_names(o: &mut NameSink) {
        o.push_indexed("ms2_intensity_ratio", NUM_MS2_IONS);
        o.push_indexed("ms1_intensity_ratio", NUM_MS1_IONS);
    }
}
