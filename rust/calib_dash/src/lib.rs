//! Terminal dashboard for debugging the RT calibration fit.

pub mod app;
pub mod bitset;
pub mod frames;
pub mod metrics;
pub mod recording;
pub mod ui;

pub use app::{
    App,
    CalibDash,
    Flow,
    Layer,
    PauseAction,
    Stage,
    Stepper,
    Tab,
    ToleranceSummary,
};
pub use bitset::BitSet;
/// Re-exported so a crate that records a fit through `FitRecording` can spell
/// `fit_with`'s options without depending on `calibrt` itself.
pub use calibrt::ObserveOpts;
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
