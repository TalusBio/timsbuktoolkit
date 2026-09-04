mod decoy;
mod digest;
pub mod query_item;
pub mod sequence;

pub use decoy::DecoyMarking;
pub use digest::{
    ProteinSlice,
    deduplicate_digests,
};
pub use query_item::ExpectedIntensities;
/// What a load decides lives in timsquery, since the arena's `seal` is what
/// resolves the decoy half of it. Re-exported so the cli's config and flags name
/// them through `timsseek`.
pub use timsquery::models::capabilities::{
    DecoyPolicy,
    LoadPolicy,
    UnannotatedPeaks,
};
