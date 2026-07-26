//! Terminal dashboard for inspecting a rescoring run: score distributions,
//! target/decoy separation, FDR curve and decoy calibration.

pub mod view;

pub use view::{
    FoldImportance,
    RescoreView,
    ViewError,
};
