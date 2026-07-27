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
//! # let view: rescore_dash::RescoreView = unimplemented!();
//! if rescore_dash::available() {
//!     match rescore_dash::Dashboard::build(&view) {
//!         Ok(dash) => {
//!             drop(view); // the matrix is dead from here on
//!             rescore_dash::run(dash)?;
//!         }
//!         Err(e) => tracing::warn!("rescore dashboard input rejected: {e}"),
//!     }
//! }
//! # Ok(())
//! # }
//! ```

mod cycle;

pub mod app;
pub mod curves;
pub mod precompute;
pub mod stats;
pub mod transform;
pub mod ui;
pub mod view;

/// What a caller embedding the dashboard needs. The panel-level types
/// (`App`, `Tab`, `Flow`, `FeatureColumn`, `Axis`) stay behind their modules —
/// nothing outside this crate drives the TUI a frame at a time.
pub use app::{
    available,
    run,
};
pub use precompute::{
    DEFAULT_SAMPLE,
    Dashboard,
};
pub use view::{
    RescoreView,
    ViewError,
};
