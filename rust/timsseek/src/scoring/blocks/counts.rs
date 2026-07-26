//! Count families. `ApexCounts` are peak-shape counts at the apex; the single
//! `FinalizeCounts` field is the number of scored fragment ions (known only at
//! finalize).

use timsseek_macros::ScoreBlock;

/// Stage: apex (peak-shape counts).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct ApexCounts {
    #[feat(raw)]
    pub rising_cycles: u8,
    #[feat(raw)]
    pub falling_cycles: u8,
    #[feat(raw)]
    pub npeaks: u8,
}

/// Stage: finalize (number of scored fragment ions).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct FinalizeCounts {
    #[feat(raw)]
    pub n_scored_fragments: u8,
}
