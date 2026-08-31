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
/// The decoy decision lives in timsquery, since the arena's `seal` is what
/// resolves it. Re-exported so `timsseek::DecoyPolicy` keeps working for the
/// cli's config and flags.
pub use timsquery::models::capabilities::DecoyPolicy;
