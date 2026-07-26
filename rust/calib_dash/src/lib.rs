//! Terminal dashboard for debugging the RT calibration fit.

pub mod bitset;
pub mod frames;
pub mod recording;

pub use bitset::BitSet;
pub use frames::{
    CalibrantPoint,
    FrameIndex,
    FrameStore,
};
pub use recording::{
    DpDecision,
    FitRecording,
};
