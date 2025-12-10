//! UI module - all user interface components and panels.

pub mod components;
pub mod panel_trait;
pub mod panels;
pub mod tolerance_editor;

pub use panel_trait::{
    CommandSink,
    Panel,
    PanelContext,
};
