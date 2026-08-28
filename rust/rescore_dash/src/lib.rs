//! Terminal dashboard for inspecting a rescoring run: score distributions,
//! target/decoy separation, FDR curve and decoy calibration.
//!
//! Two steps, deliberately separate. [`Dashboard::build`] materializes
//! everything on screen from a [`RescoreView`] -- see the `precompute` module
//! for what is
//! exact and what is sampled -- and [`run`] opens the TUI over the result.
//! Splitting them lets the caller drop the feature matrix, gigabytes at a
//! realistic library size, before the TUI blocks for as long as the user leaves
//! it open.

// Private, so `pub` inside them means "visible to the rest of this crate" and
// rustc can prove the rest dead. The whole public surface is the re-exports
// below; nothing outside this crate drives the TUI a frame at a time.
mod app;
mod curves;
mod cycle;
mod labels;
mod precompute;
mod stats;
mod transform;
mod ui;
mod view;

pub use app::run;
pub use precompute::{
    DEFAULT_SAMPLE,
    Dashboard,
};
pub use view::{
    RescoreView,
    ThresholdRow,
    ViewError,
};
