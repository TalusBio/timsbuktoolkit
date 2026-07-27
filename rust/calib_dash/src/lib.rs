//! Terminal dashboard for debugging the RT calibration fit.

pub mod app;
pub mod bitset;
pub mod frames;
pub mod metrics;
pub mod recording;

pub use app::{
    App,
    Flow,
    PauseAction,
    Stage,
    Stepper,
    Tab,
};
pub use bitset::BitSet;
pub use frames::{
    CalibrantPoint,
    FrameIndex,
    FrameStore,
};
pub use metrics::BatchMetrics;
pub use recording::{
    DpDecision,
    FitRecording,
};
