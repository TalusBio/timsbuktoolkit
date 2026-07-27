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
/// Re-exported so a crate that records a fit through `FitRecording` can spell
/// `fit_with`'s options without depending on `calibrt` itself.
pub use calibrt::ObserveOpts;
pub use frames::CalibrantPoint;
pub use recording::FitRecording;
