use eframe::egui;
use egui::Color32;
use egui_dock::{
    DockArea,
    DockState,
    Style,
    TabViewer,
};
use std::path::PathBuf;
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::tolerance::Tolerance;

use crate::chromatogram_processor::SmoothingMethod;
use crate::cli::Cli;
use crate::computed_state::ComputedState;
use crate::file_loader::{
    ElutionGroupData,
    FileLoader,
};
use crate::plot_renderer::AutoZoomMode;
use crate::ui::panels::{
    ConfigPanel,
    SpectrumPanel,
    TablePanel,
};

/// Pane types for the tile layout
#[derive(Debug, Clone, Copy, serde::Deserialize, serde::Serialize)]
enum Pane {
    ConfigPanel,
    TablePanel,
    MS2Spectrum,
    PrecursorPlot,
    FragmentPlot,
    ScoresPlot,
}

/// State to be persisted across restarts
#[derive(serde::Deserialize, serde::Serialize)]
struct PersistentState {
    file_loader: FileLoader,
    ui_state: UiState,
    tolerance: Tolerance,
    smoothing: SmoothingMethod,
    dock_state: DockState<Pane>,
}

/// State of indexed raw data loading
#[derive(Debug, Default)]
pub enum IndexedDataState {
    /// No data loaded
    #[default]
    None,
    /// Currently loading data from path
    Loading(PathBuf),
    /// Loading failed with error message
    Failed(PathBuf, String),
    /// Data successfully loaded
    Loaded { index: Arc<IndexedTimstofPeaks> },
}

/// Domain/data state - represents loaded data and analysis parameters
#[derive(Debug, Default)]
pub struct DataState {
    /// Loaded elution groups
    pub elution_groups: Option<ElutionGroupData>,
    /// Path from which the elution groups were loaded
    pub elution_groups_source: Option<PathBuf>,
    /// Indexed timsTOF data with loading state
    pub indexed_data: IndexedDataState,
    /// Tolerance settings
    pub tolerance: Tolerance,
    /// Smoothing method configuration
    pub smoothing: SmoothingMethod,
    /// Auto-zoom mode for plots
    pub auto_zoom_mode: AutoZoomMode,
}

/// UI-specific state - transient UI state that doesn't affect data
#[derive(Debug, serde::Deserialize, serde::Serialize, Default)]
pub struct UiState {
    /// Filter text for precursor table
    pub table_filter: String,
    /// Currently selected elution group index
    pub selected_index: Option<usize>,
    /// Vim mode: search mode active
    pub search_mode: bool,
    /// Vim mode: search input buffer
    pub search_input: String,
}

/// Main application state
pub struct ViewerApp {
    /// File loader for handling file dialogs and loading
    file_loader: FileLoader,

    /// Domain/data state
    data: DataState,

    /// UI-specific state
    ui: UiState,

    /// Computed/cached state
    computed: ComputedState,

    /// Dock state for layout management
    dock_state: DockState<Pane>,

    /// UI Panels
    config_panel: ConfigPanel,
    table_panel: TablePanel,
    spectrum_panel: SpectrumPanel,
}

impl ViewerApp {
    /// Create a new test instance without eframe context
    #[cfg(test)]
    pub fn new_test() -> Self {
        let tabs = vec![
            Pane::ConfigPanel,
            Pane::TablePanel,
            Pane::PrecursorPlot,
            Pane::FragmentPlot,
            Pane::MS2Spectrum,
        ];

        let dock_state = DockState::new(tabs);

        Self {
            file_loader: FileLoader::new(),
            data: DataState::default(),
            ui: UiState::default(),
            computed: ComputedState::default(),
            dock_state,
            config_panel: ConfigPanel::new(),
            table_panel: TablePanel::new(),
            spectrum_panel: SpectrumPanel::new(),
        }
    }

    pub fn new(cc: &eframe::CreationContext<'_>, args: &Cli) -> Self {
        // Try to load previous state
        if let Some(storage) = cc.storage {
            if let Some(state_string) = storage.get_string(eframe::APP_KEY) {
                // Try RON first (new format), then fallback to JSON (legacy)
                let state = match ron::from_str::<PersistentState>(&state_string) {
                    Ok(state) => {
                        tracing::info!("Loaded persistent state successfully (RON).");
                        Some(state)
                    }
                    Err(ron_err) => {
                        // Fallback to JSON for backward compatibility
                        match serde_json::from_str::<PersistentState>(&state_string) {
                            Ok(state) => {
                                tracing::info!("Loaded persistent state successfully (JSON).");
                                Some(state)
                            }
                            Err(json_err) => {
                                tracing::warn!(
                                    "Failed to deserialize persistent state. RON: {:?}, JSON: {:?}",
                                    ron_err,
                                    json_err
                                );
                                None
                            }
                        }
                    }
                };

                if let Some(state) = state {
                    return Self {
                        file_loader: state
                            .file_loader
                            .with_initial_paths(&args.raw_data_path, &args.elution_groups_path),
                        data: DataState {
                            tolerance: state.tolerance,
                            smoothing: state.smoothing,
                            ..DataState::default()
                        },
                        ui: state.ui_state,
                        computed: ComputedState::default(),
                        dock_state: state.dock_state,
                        config_panel: ConfigPanel::new(),
                        table_panel: TablePanel::new(),
                        spectrum_panel: SpectrumPanel::new(),
                    };
                }
            } else {
                tracing::info!("No persistent state found.");
            }
        }

        // Create initial tabs: Settings, Table, Precursors, Fragments, MS2
        let tabs = vec![
            Pane::ConfigPanel,
            Pane::TablePanel,
            Pane::PrecursorPlot,
            Pane::FragmentPlot,
            Pane::MS2Spectrum,
            Pane::ScoresPlot,
        ];

        let dock_state = DockState::new(tabs);

        Self {
            file_loader: FileLoader::new()
                .with_initial_paths(&args.raw_data_path, &args.elution_groups_path),
            data: DataState::default(),
            ui: UiState::default(),
            computed: ComputedState::default(),
            dock_state,
            config_panel: ConfigPanel::new(),
            table_panel: TablePanel::new(),
            spectrum_panel: SpectrumPanel::new(),
        }
    }

    fn handle_vim_keys(&mut self, ctx: &egui::Context) {
        if ctx.memory(|mem| mem.focused().is_some()) && !self.ui.search_mode {
            return;
        }

        ctx.input(|i| {
            if self.ui.search_mode {
                if i.key_pressed(egui::Key::Escape) {
                    self.ui.search_mode = false;
                }
                if i.key_pressed(egui::Key::Enter) {
                    self.ui.search_mode = false;
                }
                return;
            }

            if i.key_pressed(egui::Key::Slash) && !i.modifiers.any() {
                self.ui.search_mode = true;
                self.ui.search_input = self.ui.table_filter.clone();
                return;
            }

            if self.data.elution_groups.is_none() {
                return;
            }

            // Cache filtered indices - computed once and reused
            let filtered_indices = &self.table_panel.filtered_indices();
            let cursor = &mut self.ui.selected_index;

            if i.key_pressed(egui::Key::J) && !i.modifiers.any() {
                Self::move_selection_down(cursor, filtered_indices);
            }

            if i.key_pressed(egui::Key::K) && !i.modifiers.any() {
                Self::move_selection_up(cursor, filtered_indices);
            }

            if i.key_pressed(egui::Key::G) && !i.modifiers.any() {
                Self::move_selection_to_first(cursor, filtered_indices);
            }

            if i.key_pressed(egui::Key::G) && i.modifiers.shift_only() {
                Self::move_selection_to_last(cursor, filtered_indices);
            }
        });
    }

    fn generate_chromatogram(&mut self) {
        let IndexedDataState::Loaded { index } = &self.data.indexed_data else {
            return;
        };

        let selected_idx = match self.ui.selected_index {
            Some(idx) => idx,
            None => return,
        };

        if let Some(elution_groups) = &self.data.elution_groups {
            self.computed.update(
                elution_groups,
                selected_idx,
                index,
                &self.data.tolerance,
                &self.data.smoothing,
            )
        };
    }

    fn render_elution_groups_section_static(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
    ) {
        ui.label("Elution Groups:");
        if ui.button("Load Elution Groups...").clicked() {
            file_loader.open_elution_groups_dialog();
        }

        if let Some(path) = &file_loader.elution_groups_path {
            Self::display_filename(ui, path);
        }

        Self::load_elution_groups_if_needed(ui, file_loader, data);

        if let Some(egs) = &data.elution_groups {
            ui.label(format!("✓ Loaded: {} elution groups", egs.len()));
        }
    }

    fn render_raw_data_section_static(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
        _ui_state: &mut UiState,
    ) {
        ui.label("Raw Data File (.d):");
        if ui.button("Load Raw Data...").clicked() {
            file_loader.open_raw_data_dialog();
        }

        if let Some(path) = &file_loader.raw_data_path {
            Self::display_filename(ui, path);
        }

        Self::load_raw_data_if_needed(ui, file_loader, data);

        if matches!(data.indexed_data, IndexedDataState::Loaded { .. }) {
            ui.label("✓ Raw data indexed");
        }
    }

    fn render_tolerance_loading_section_static(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
    ) {
        ui.label("Tolerance Settings:");
        if ui.button("Load Tolerances...").clicked() {
            file_loader.open_tolerance_dialog();
        }

        if let Some(path) = &file_loader.tolerance_path {
            Self::display_filename(ui, path);
        }

        Self::load_tolerance_if_needed(file_loader, data);

        ui.label("Tolerance settings loaded");
    }

    fn display_filename(ui: &mut egui::Ui, path: &std::path::Path) {
        let filename = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("Unknown");
        ui.label(egui::RichText::new(filename).small().italics())
            .on_hover_text(path.display().to_string());
    }

    fn load_elution_groups_if_needed(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
    ) {
        if let Some(path) = &file_loader.elution_groups_path {
            let should_load = match &data.elution_groups_source {
                Some(current_path) => current_path != path,
                None => true,
            };

            if should_load {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Loading elution groups...");
                });
                match file_loader.load_elution_groups(path) {
                    Ok(egs) => {
                        tracing::info!("Loaded {} elution groups", egs.len());
                        data.elution_groups = Some(egs);
                        data.elution_groups_source = Some(path.clone());
                    }
                    Err(e) => {
                        tracing::error!("Failed to load elution groups: {:?}", e);
                    }
                }
            }
        }
    }

    fn load_raw_data_if_needed(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
    ) {
        if let Some(path) = &file_loader.raw_data_path {
            // Transition from None to Loading
            if matches!(data.indexed_data, IndexedDataState::None) {
                data.indexed_data = IndexedDataState::Loading(path.clone());
                ui.ctx().request_repaint();
            }
        }

        match &data.indexed_data {
            IndexedDataState::Loading(path) => {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Indexing raw data... (this may take 10-30 seconds)");
                });

                match file_loader.load_raw_data(path) {
                    Ok(index) => {
                        data.indexed_data = IndexedDataState::Loaded { index };
                        file_loader.raw_data_path = None;
                        tracing::info!("Raw data indexing completed");
                    }
                    Err(e) => {
                        let error_msg = format!("{:?}", e);
                        tracing::error!("Failed to load raw data: {}", error_msg);
                        data.indexed_data = IndexedDataState::Failed(path.clone(), error_msg);
                        file_loader.raw_data_path = None;
                    }
                }
            }
            IndexedDataState::Failed(path, error) => {
                ui.label(
                    egui::RichText::new(format!(
                        "Failed to load raw data from {}:",
                        path.display()
                    ))
                    .color(egui::Color32::RED),
                );
                ui.label(egui::RichText::new(error).color(egui::Color32::RED).small());
                if ui.button("Clear Error").clicked() {
                    data.indexed_data = IndexedDataState::None;
                }
            }
            _ => {}
        }
    }

    fn load_tolerance_if_needed(file_loader: &mut FileLoader, data: &mut DataState) {
        if let Some(path) = &file_loader.tolerance_path {
            match file_loader.load_tolerance(path) {
                Ok(tol) => {
                    data.tolerance = tol;
                }
                Err(e) => {
                    tracing::error!("Failed to load tolerance: {:?}", e);
                }
            }
        }
    }

    fn move_selection_down(cursor: &mut Option<usize>, filtered_indices: &[usize]) {
        if filtered_indices.is_empty() {
            return;
        }

        match cursor {
            None => {
                Self::select_elution_group(cursor, filtered_indices[0]);
            }
            Some(current) => {
                if let Some(pos) = filtered_indices.iter().position(|&idx| idx == *current)
                    && pos + 1 < filtered_indices.len()
                {
                    Self::select_elution_group(cursor, filtered_indices[pos + 1]);
                }
            }
        }
    }

    fn move_selection_up(cursor: &mut Option<usize>, filtered_indices: &[usize]) {
        if filtered_indices.is_empty() {
            return;
        }

        match cursor {
            None => {
                Self::select_elution_group(cursor, filtered_indices[0]);
            }
            Some(current) => {
                if let Some(pos) = filtered_indices.iter().position(|&idx| idx == *current)
                    && pos > 0
                {
                    Self::select_elution_group(cursor, filtered_indices[pos - 1]);
                }
            }
        }
    }

    fn move_selection_to_first(cursor: &mut Option<usize>, filtered_indices: &[usize]) {
        if !filtered_indices.is_empty() {
            Self::select_elution_group(cursor, filtered_indices[0]);
        }
    }

    fn move_selection_to_last(cursor: &mut Option<usize>, filtered_indices: &[usize]) {
        if !filtered_indices.is_empty() {
            Self::select_elution_group(cursor, filtered_indices[filtered_indices.len() - 1]);
        }
    }

    fn select_elution_group(cursor: &mut Option<usize>, idx: usize) {
        *cursor = Some(idx);
    }
}

impl eframe::App for ViewerApp {
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        tracing::info!("Saving application state...");
        let state = PersistentState {
            file_loader: FileLoader {
                elution_groups_path: self.file_loader.elution_groups_path.clone(),
                raw_data_path: self.file_loader.raw_data_path.clone(),
                tolerance_path: self.file_loader.tolerance_path.clone(),
            },
            ui_state: UiState {
                table_filter: self.ui.table_filter.clone(),
                selected_index: self.ui.selected_index,
                search_mode: self.ui.search_mode,
                search_input: self.ui.search_input.clone(),
            },
            tolerance: self.data.tolerance.clone(),
            smoothing: self.data.smoothing,
            dock_state: self.dock_state.clone(),
        };

        if let Ok(value) = ron::to_string(&state) {
            storage.set_string(eframe::APP_KEY, value);
        } else {
            tracing::error!("Failed to serialize state to RON");
        }
    }

    fn auto_save_interval(&self) -> std::time::Duration {
        std::time::Duration::from_secs(5)
    }

    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.handle_vim_keys(ctx);

        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.heading("TimsQuery Viewer");
            ui.separator();
        });
        self.generate_chromatogram();

        let mut tab_viewer = AppTabViewer {
            file_loader: &mut self.file_loader,
            data: &mut self.data,
            ui: &mut self.ui,
            computed: &mut self.computed,
            left_panel: &mut self.config_panel,
            table_panel: &mut self.table_panel,
            spectrum_panel: &mut self.spectrum_panel,
        };

        egui::CentralPanel::default().show(ctx, |ui| {
            DockArea::new(&mut self.dock_state)
                .style(Style::from_egui(ui.style().as_ref()))
                .show_inside(ui, &mut tab_viewer);
        });
    }
}

struct AppTabViewer<'a> {
    file_loader: &'a mut FileLoader,
    data: &'a mut DataState,
    ui: &'a mut UiState,
    computed: &'a mut ComputedState,
    left_panel: &'a mut ConfigPanel,
    table_panel: &'a mut TablePanel,
    spectrum_panel: &'a mut SpectrumPanel,
}

impl<'a> AppTabViewer<'a> {
    fn render_left_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Data Loading");
        ui.separator();

        ViewerApp::render_elution_groups_section_static(ui, self.file_loader, self.data);
        ui.add_space(10.0);

        ViewerApp::render_raw_data_section_static(ui, self.file_loader, self.data, self.ui);
        ui.add_space(10.0);

        ViewerApp::render_tolerance_loading_section_static(ui, self.file_loader, self.data);

        ui.add_space(20.0);
        ui.separator();

        self.left_panel.render(
            ui,
            &mut self.data.tolerance,
            &mut self.data.smoothing,
            &mut self.data.auto_zoom_mode,
        );
    }
}

impl<'a> TabViewer for AppTabViewer<'a> {
    type Tab = Pane;

    fn title(&mut self, tab: &mut Self::Tab) -> egui::WidgetText {
        match tab {
            Pane::ConfigPanel => self.left_panel.title().into(),
            Pane::TablePanel => self.table_panel.title().into(),
            Pane::MS2Spectrum => self.spectrum_panel.title().into(),
            Pane::PrecursorPlot => "Precursors".into(),
            Pane::FragmentPlot => "Fragments".into(),
            Pane::ScoresPlot => "Scores".into(),
        }
    }

    fn ui(&mut self, ui: &mut egui::Ui, tab: &mut Self::Tab) {
        let mode = self.data.auto_zoom_mode;

        // TODO: figure out how to prevent this allocation per frame...
        let ref_lines: Vec<(String, f64, Color32)> = self
            .computed
            .reference_lines()
            .iter()
            .map(|(k, v)| (k.clone(), v.0, v.1))
            .collect();
        match tab {
            Pane::ConfigPanel => {
                // Wrap settings in a scroll area
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        self.render_left_panel(ui);
                    });
            }
            Pane::TablePanel => {
                self.table_panel.render(
                    ui,
                    &self.data.elution_groups,
                    self.ui.search_mode,
                    &mut self.ui.search_input,
                    &mut self.ui.selected_index,
                );
            }
            Pane::MS2Spectrum => {
                self.spectrum_panel.render(
                    ui,
                    &self.computed.ms2_spectrum,
                    &self.computed.expected_intensities,
                );
            }
            Pane::PrecursorPlot => {
                if let Some(chromatogram) = &self.computed.chromatogram_lines {
                    // Use shared link_id for synchronized X-axis with Fragments
                    let click_response = crate::plot_renderer::render_chromatogram_plot(
                        ui,
                        chromatogram,
                        crate::plot_renderer::PlotMode::PrecursorsOnly,
                        Some("precursor_fragment_x_axis"),
                        true,
                        &mut self.computed.auto_zoom_frame_counter,
                        &mode,
                        &ref_lines,
                    );
                    if let Some(clicked_rt) = click_response {
                        self.computed.clicked_rt = Some(clicked_rt);
                    }
                } else if self.ui.selected_index.is_some() {
                    ui.label("Generating chromatogram...");
                } else {
                    ui.label("Select a precursor from the table to view chromatogram");
                }
            }
            Pane::FragmentPlot => {
                if let Some(chromatogram) = &self.computed.chromatogram_lines {
                    // Use shared link_id for synchronized X-axis with Precursors
                    let response = crate::plot_renderer::render_chromatogram_plot(
                        ui,
                        chromatogram,
                        crate::plot_renderer::PlotMode::FragmentsOnly,
                        Some("precursor_fragment_x_axis"),
                        false,
                        &mut self.computed.auto_zoom_frame_counter,
                        &mode,
                        &ref_lines,
                    );
                    if let Some(clicked_rt) = response {
                        self.computed.clicked_rt = Some(clicked_rt);
                    }
                } else if self.ui.selected_index.is_some() {
                    ui.label("Generating chromatogram...");
                } else {
                    ui.label("Select a precursor to view fragment traces");
                }
            }
            Pane::ScoresPlot => {
                if let Some(score_lines) = &self.computed.score_lines {
                    let response = score_lines.render(
                        ui,
                        Some("precursor_fragment_x_axis"),
                        &mut self.computed.auto_zoom_frame_counter,
                        &mode,
                        &ref_lines,
                    );
                    if let Some(clicked_rt) = response {
                        self.computed.clicked_rt = Some(clicked_rt);
                    }
                } else if self.ui.selected_index.is_some() {
                    ui.label("Generating score plot...");
                } else {
                    ui.label("Select a precursor to view score traces");
                }
            }
        }
    }

    fn is_closeable(&self, _tab: &Self::Tab) -> bool {
        false // Tabs cannot be closed
    }
}
