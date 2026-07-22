//! Count families. `ApexCounts` are peak-shape counts at the apex; the single
//! `FinalizeCounts` field is the number of scored fragment ions (known only at
//! finalize).

use crate::score_block;

score_block! {
    /// Stage: apex (peak-shape counts).
    pub struct ApexCounts {
        #[raw] pub rising_cycles: u8,
        #[raw] pub falling_cycles: u8,
        #[raw] pub npeaks: u8,
    }
}

score_block! {
    /// Stage: finalize (number of scored fragment ions).
    pub struct FinalizeCounts {
        #[raw] pub n_scored_fragments: u8,
    }
}
