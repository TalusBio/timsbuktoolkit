//! apex_sim library surface: synthetic extraction generation + apex scoring,
//! reused by the eframe app (`main.rs`) and the `sensitivity` bench example.

#[cfg(feature = "gui")]
pub mod app;
pub mod bench;
#[cfg(feature = "gui")]
pub mod plots;
pub mod scorer;
pub mod sim;
