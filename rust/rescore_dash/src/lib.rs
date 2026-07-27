//! Terminal dashboard for inspecting a rescoring run: score distributions,
//! target/decoy separation, FDR curve and decoy calibration.
//!
//! Two steps, deliberately separate. [`Dashboard::build`] materializes
//! everything on screen from a [`RescoreView`] — see [`precompute`] for what is
//! exact and what is sampled — and [`run`] opens the TUI over the result.
//! Splitting them lets the caller drop the feature matrix, which is gigabytes
//! at a realistic library size, before the TUI blocks for as long as the user
//! leaves it open.
//!
//! ```no_run
//! # fn main() -> std::io::Result<()> {
//! # fn feature_matrix() -> Vec<f64> { unimplemented!() }
//! # let sample = rescore_dash::DEFAULT_SAMPLE;
//! let dash = if rescore_dash::available() {
//!     let matrix = feature_matrix(); // gigabytes
//!     # let view: rescore_dash::RescoreView = unimplemented!();
//!     rescore_dash::Dashboard::build(&view, sample)
//!         .inspect_err(|e| tracing::warn!("rescore dashboard input rejected: {e}"))
//!         .ok()
//! } else {
//!     None
//! };
//! // `matrix` is out of scope here; only the dashboard survives into the TUI.
//! if let Some(dash) = dash {
//!     rescore_dash::run(dash)?;
//! }
//! # Ok(())
//! # }
//! ```

// Private, so `pub` inside them means "visible to the rest of this crate" and
// rustc can prove the rest dead. The whole public surface is the re-exports
// below; nothing outside this crate drives the TUI a frame at a time.
mod app;
mod curves;
mod cycle;
mod precompute;
mod stats;
mod transform;
mod ui;
mod view;

pub use app::{
    available,
    run,
};
pub use curves::ThresholdRow;
pub use precompute::{
    DEFAULT_SAMPLE,
    Dashboard,
};
pub use view::{
    RescoreView,
    ViewError,
};
