// Panel modules for UI organization
// Each panel is responsible for rendering its portion of the UI
// and returning commands for state changes

pub mod left_panel;
pub mod plot_panel;
pub mod spectrum_panel;
pub mod table_panel;

pub use left_panel::LeftPanel;
pub use plot_panel::PlotPanel;
pub use spectrum_panel::SpectrumPanel;
pub use table_panel::TablePanel;
