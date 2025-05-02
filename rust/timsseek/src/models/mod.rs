mod decoy;
mod digest;

pub use decoy::DecoyMarking;
pub use digest::{
    DigestSlice,
    deduplicate_digests,
};
