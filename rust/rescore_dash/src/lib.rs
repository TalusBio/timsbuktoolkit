//! Terminal dashboard for inspecting a rescoring run: score distributions,
//! target/decoy separation, FDR curve and decoy calibration.

pub mod app;
pub mod curves;
pub mod stats;
pub mod transform;
pub mod ui;
pub mod view;

pub use app::{
    App,
    Flow,
    SortKey,
    Tab,
    run,
};
pub use view::{
    FoldImportance,
    RescoreView,
    ViewError,
};
