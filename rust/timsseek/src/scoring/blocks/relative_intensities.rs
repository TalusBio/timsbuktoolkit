//! Relative-intensity arrays (per-ion intensity relative to the MS-level
//! total). Built from the secondary-query collector's sorted values.
//!
//! Same array-family shape as [`super::ion_errors`]: module-scope name tables +
//! `FeatSink::push_slice`. Column and feature order both run MS2→MS1.
//!
//! NOTE: adding/renaming a name here also requires updating the parquet schema
//! golden (`parquet_writer`) and the feature-name golden (`ml::qvalues`).

use crate::scoring::apex_finding::RelativeIntensityCollector;
use crate::scoring::blocks::{
    ColSink,
    FeatSink,
    ScoreBlock,
};
use crate::scoring::{
    NUM_MS1_IONS,
    NUM_MS2_IONS,
};

const MS2_INTENSITY_RATIO: [&str; NUM_MS2_IONS] = [
    "ms2_intensity_ratio_0",
    "ms2_intensity_ratio_1",
    "ms2_intensity_ratio_2",
    "ms2_intensity_ratio_3",
    "ms2_intensity_ratio_4",
    "ms2_intensity_ratio_5",
    "ms2_intensity_ratio_6",
];
const MS1_INTENSITY_RATIO: [&str; NUM_MS1_IONS] = [
    "ms1_intensity_ratio_0",
    "ms1_intensity_ratio_1",
    "ms1_intensity_ratio_2",
];

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
        o.push_slice(&MS2_INTENSITY_RATIO, &self.ms2_intensity_ratios);
        o.push_slice(&MS1_INTENSITY_RATIO, &self.ms1_intensity_ratios);
    }
}
