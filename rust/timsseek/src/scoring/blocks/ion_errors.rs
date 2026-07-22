//! Per-ion m/z and mobility error arrays. Hand-written array family: the
//! `columns`/`features` bodies walk const name tables. Changing an ion count
//! means editing `NUM_MS*_IONS` *and* the four name tables below (their lengths
//! are checked against `NUM_MS*_IONS` at compile time, so a mismatch is a build
//! error, not silent drift).
//!
//! Column and feature names are identical per quantity (`ms2_mz_error_0`, …) —
//! same table drives both surfaces.
//!
//! NOTE: adding/renaming a name here also requires updating the parquet schema
//! golden (`parquet_writer`) and the feature-name golden (`ml::qvalues`).

use crate::scoring::blocks::{
    ColSink,
    FeatSink,
    ScoreBlock,
};
use crate::scoring::offsets::MzMobilityOffsets;
use crate::scoring::{
    NUM_MS1_IONS,
    NUM_MS2_IONS,
};

const MS2_MZ_ERROR: [&str; NUM_MS2_IONS] = [
    "ms2_mz_error_0",
    "ms2_mz_error_1",
    "ms2_mz_error_2",
    "ms2_mz_error_3",
    "ms2_mz_error_4",
    "ms2_mz_error_5",
    "ms2_mz_error_6",
];
const MS2_MOB_ERROR: [&str; NUM_MS2_IONS] = [
    "ms2_mobility_error_0",
    "ms2_mobility_error_1",
    "ms2_mobility_error_2",
    "ms2_mobility_error_3",
    "ms2_mobility_error_4",
    "ms2_mobility_error_5",
    "ms2_mobility_error_6",
];
const MS1_MZ_ERROR: [&str; NUM_MS1_IONS] = ["ms1_mz_error_0", "ms1_mz_error_1", "ms1_mz_error_2"];
const MS1_MOB_ERROR: [&str; NUM_MS1_IONS] = [
    "ms1_mobility_error_0",
    "ms1_mobility_error_1",
    "ms1_mobility_error_2",
];

#[derive(Debug, Clone, Copy, serde::Serialize)]
pub struct IonErrors {
    pub ms2_mz_errors: [f32; NUM_MS2_IONS],
    pub ms2_mobility_errors: [f32; NUM_MS2_IONS],
    pub ms1_mz_errors: [f32; NUM_MS1_IONS],
    pub ms1_mobility_errors: [f32; NUM_MS1_IONS],
}

impl IonErrors {
    pub fn compute(offsets: &MzMobilityOffsets) -> Self {
        Self {
            ms2_mz_errors: offsets.ms2_mz_errors(),
            ms2_mobility_errors: offsets.ms2_mobility_errors(),
            ms1_mz_errors: offsets.ms1_mz_errors(),
            ms1_mobility_errors: offsets.ms1_mobility_errors(),
        }
    }

    pub fn neutralize(&mut self) {
        self.ms1_mobility_errors.fill(f32::NAN);
        self.ms2_mobility_errors.fill(f32::NAN);
    }

    pub fn sample() -> Self {
        Self {
            ms2_mz_errors: [0.5; NUM_MS2_IONS],
            ms2_mobility_errors: [0.5; NUM_MS2_IONS],
            ms1_mz_errors: [0.5; NUM_MS1_IONS],
            ms1_mobility_errors: [0.5; NUM_MS1_IONS],
        }
    }
}

impl ScoreBlock for IonErrors {
    fn columns(&self, o: &mut ColSink) {
        o.f32_array("ms2_mz_error", &self.ms2_mz_errors);
        o.f32_array("ms2_mobility_error", &self.ms2_mobility_errors);
        o.f32_array("ms1_mz_error", &self.ms1_mz_errors);
        o.f32_array("ms1_mobility_error", &self.ms1_mobility_errors);
    }

    fn features(&self, o: &mut FeatSink) {
        o.push_slice(&MS2_MZ_ERROR, &self.ms2_mz_errors);
        o.push_slice(&MS2_MOB_ERROR, &self.ms2_mobility_errors);
        o.push_slice(&MS1_MZ_ERROR, &self.ms1_mz_errors);
        o.push_slice(&MS1_MOB_ERROR, &self.ms1_mobility_errors);
    }
}
