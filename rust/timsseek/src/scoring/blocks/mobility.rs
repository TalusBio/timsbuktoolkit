//! Mobility family. Finalize-stage: observed 1/K0 and the MS1↔MS2 mobility
//! delta, both from the secondary-query offset collector.
//!
//! `neutralize` NaNs every mobility-derived field, INCLUDING the derived square
//! `sq_delta_ms1_ms2_mobility`, so a future move of the square cannot desync
//! from its source (see [`crate::scoring::results::ScoringFields::neutralize_mobility`]).

use crate::score_block;
use crate::scoring::offsets::MzMobilityOffsets;

score_block! {
    /// Stage: finalize (from `MzMobilityOffsets`).
    pub struct Mobility {
        #[raw] pub obs_mobility: f32,
        #[feat(raw, abs, isna)] pub delta_ms1_ms2_mobility: f32,
        #[raw] pub sq_delta_ms1_ms2_mobility: f32,
    }
}

impl Mobility {
    pub fn compute(offsets: &MzMobilityOffsets) -> Self {
        let (ms1_err, ms2_err) = offsets.avg_delta_mobs();
        let cum_err = ms1_err + ms2_err;
        let obs_mobility =
            (offsets.ref_mobility + cum_err.mean_mobility().unwrap_or(f64::NAN)) as f32;
        let d_err = match (ms1_err.mean_mobility(), ms2_err.mean_mobility()) {
            (Ok(ms1_mob), Ok(ms2_mob)) => ms1_mob - ms2_mob,
            _ => f64::NAN,
        };
        let delta = d_err as f32;
        Self {
            obs_mobility,
            delta_ms1_ms2_mobility: delta,
            sq_delta_ms1_ms2_mobility: delta * delta,
        }
    }

    pub fn neutralize(&mut self) {
        self.obs_mobility = f32::NAN;
        self.delta_ms1_ms2_mobility = f32::NAN;
        self.sq_delta_ms1_ms2_mobility = f32::NAN;
    }
}
