//! Retention-time family. Finalize-stage: reads library/calibrated RT from
//! metadata and the observed apex RT (`inp.apex.retention_time_ms`), so it is
//! cheap arithmetic on an apex scalar — not a straddle.

use crate::score_block;
use crate::scoring::apex_finding::PeptideMetadata;

score_block! {
    /// Stage: finalize (reads apex retention time).
    pub struct Rt {
        #[col_only] pub library_rt: f32,
        #[feat(round, raw)] pub calibrated_rt_seconds: f32,
        #[raw] pub obs_rt_seconds: f32,
        #[raw] pub calibrated_sq_delta_rt: f32,
    }
}

impl Rt {
    pub fn compute(metadata: &PeptideMetadata, obs_rt_seconds: f32) -> Self {
        let calibrated_delta_rt = obs_rt_seconds - metadata.calibrated_rt_seconds;
        Self {
            library_rt: metadata.library_rt,
            calibrated_rt_seconds: metadata.calibrated_rt_seconds,
            obs_rt_seconds,
            calibrated_sq_delta_rt: calibrated_delta_rt * calibrated_delta_rt,
        }
    }
}
