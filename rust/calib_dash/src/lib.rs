//! Terminal dashboard for debugging the RT calibration fit.

pub(crate) mod app;
pub(crate) mod frames;
pub(crate) mod metrics;
pub(crate) mod recording;
pub(crate) mod ui;

pub use app::{
    CalibDash,
    Flow,
    ToleranceSummary,
};
pub use frames::{
    CalibrantPoint,
    DEFAULT_RUN_BUDGET_BYTES,
    REPLAY_BUDGET_BYTES,
};
pub use recording::FitRecording;
