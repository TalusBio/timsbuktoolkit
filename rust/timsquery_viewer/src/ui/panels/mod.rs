// Panel modules for UI organization
// Each panel is responsible for rendering its portion of the UI
// and returning commands for state changes

pub mod left_settings_panel;
pub mod precursor_table_panel;
pub mod spectrum_display_panel;

pub use left_settings_panel::ConfigPanel;
pub use precursor_table_panel::TablePanel;
pub use spectrum_display_panel::SpectrumPanel;
