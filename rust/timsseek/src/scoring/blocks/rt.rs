//! Retention-time family. Finalize-stage: reads library/calibrated RT from
//! metadata and the observed apex RT (`inp.apex.retention_time_ms`), so it is
//! cheap arithmetic on an apex scalar — not a straddle.

use timsseek_macros::ScoreBlock;

use crate::scoring::apex_finding::PeptideMetadata;

/// Stage: finalize (reads apex retention time).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct Rt {
    pub library_rt: f32,
    #[feat(round, raw, linear = false)]
    pub calibrated_rt_seconds: f32,
    #[feat(raw, linear = false)]
    pub obs_rt_seconds: f32,
    #[feat(raw, abs)]
    pub calibrated_delta_rt: f32,
}

impl Rt {
    pub fn compute(metadata: &PeptideMetadata, obs_rt_seconds: f32) -> Self {
        Self {
            library_rt: metadata.library_rt,
            calibrated_rt_seconds: metadata.calibrated_rt_seconds,
            obs_rt_seconds,
            calibrated_delta_rt: obs_rt_seconds - metadata.calibrated_rt_seconds,
        }
    }
}
