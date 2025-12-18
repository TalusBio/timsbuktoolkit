use eframe::egui;
use egui_dock::{
    DockArea,
    DockState,
    Style,
    TabViewer,
};
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::tolerance::Tolerance;

use crate::chromatogram_processor::{
    self,
    ChromatogramOutput,
    SmoothingMethod,
};
use crate::domain::ChromatogramService;
use crate::file_loader::{
    ElutionGroupData,
    FileLoader,
};
use crate::plot_renderer::{
    ChromatogramLines,
    MS2Spectrum,
};
use crate::ui::panels::{
    LeftPanel,
    SpectrumPanel,
    TablePanel,
};
use crate::ui::{
    Panel,
    PanelContext,
};

/// Commands that trigger state changes in the application
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub enum AppCommand {
    /// Regenerate the chromatogram for the current selection
    RegenerateChromatogram,
    /// Select a specific elution group by index
    SelectElutionGroup(usize),
    /// Update the tolerance settings
    UpdateTolerance,
    /// Update the smoothing settings
    UpdateSmoothing,
    /// Query MS2 spectrum at given RT (seconds)
    QueryMS2Spectrum(f64),
}

/// Pane types for the tile layout
#[derive(Debug, Clone, Copy, serde::Deserialize, serde::Serialize)]
enum Pane {
    LeftPanel,
    TablePanel,
    MS2Spectrum,
    PrecursorPlot,
    FragmentPlot,
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
    Loaded {
        index: Arc<IndexedTimstofPeaks>,
        ms1_rts: Arc<[u32]>,
    },
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
}

/// UI-specific state - transient UI state that doesn't affect data
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct UiState {
    /// Filter text for precursor table
    pub table_filter: String,
    /// Currently selected elution group index
    pub selected_index: Option<usize>,
    /// Vim mode: search mode active
    pub search_mode: bool,
    /// Vim mode: search input buffer
    pub search_input: String,
    /// Master toggle for MS2 spectrum feature
    pub show_ms2_spectrum: bool,
}

impl Default for UiState {
    fn default() -> Self {
        Self {
            table_filter: String::new(),
            selected_index: None,
            search_mode: false,
            search_input: String::new(),
            show_ms2_spectrum: true, // Enable MS2 by default
        }
    }
}

/// Computed/cached state - derived from data and UI state
#[derive(Debug, Default)]
pub struct ComputedState {
    /// Computed chromatogram for the selected elution group (plot data)
    pub chromatogram: Option<ChromatogramLines>,
    /// X-axis bounds to apply on next plot render (min_rt, max_rt)
    pub chromatogram_x_bounds: Option<(f64, f64)>,
    /// Raw chromatogram output data (for MS2 extraction)
    pub chromatogram_output: Option<ChromatogramOutput>,
    /// Computed MS2 spectrum at selected RT
    pub ms2_spectrum: Option<MS2Spectrum>,
    /// Frame counter for auto-zooming the chromatogram
    pub auto_zoom_frame_counter: u8,
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

    /// Pending commands to be executed
    pending_commands: Vec<AppCommand>,

    /// Dock state for layout management
    dock_state: DockState<Pane>,

    /// UI Panels
    left_panel: LeftPanel,
    table_panel: TablePanel,
    spectrum_panel: SpectrumPanel,

    /// Optional session log writer
    session_log: Option<Box<dyn Write + Send + Sync>>,
}

impl ViewerApp {
    /// Create a new test instance without eframe context
    #[cfg(test)]
    pub fn new_test() -> Self {
        let tabs = vec![
            Pane::LeftPanel,
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
            pending_commands: Vec::new(),
            dock_state,
            left_panel: LeftPanel::new(),
            table_panel: TablePanel::new(),
            spectrum_panel: SpectrumPanel::new(),
            session_log: None,
        }
    }

    pub fn new(
        cc: &eframe::CreationContext<'_>,
        session_log: Option<Box<dyn Write + Send + Sync>>,
    ) -> Self {
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
                        file_loader: state.file_loader,
                        data: DataState {
                            tolerance: state.tolerance,
                            smoothing: state.smoothing,
                            ..DataState::default()
                        },
                        ui: state.ui_state,
                        computed: ComputedState::default(),
                        pending_commands: Vec::new(),
                        dock_state: state.dock_state,
                        left_panel: LeftPanel::new(),
                        table_panel: TablePanel::new(),
                        spectrum_panel: SpectrumPanel::new(),
                        session_log,
                    };
                }
            } else {
                tracing::info!("No persistent state found.");
            }
        }

        // Create initial tabs: Settings, Table, Precursors, Fragments, MS2
        let tabs = vec![
            Pane::LeftPanel,
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
            pending_commands: Vec::new(),
            dock_state,
            left_panel: LeftPanel::new(),
            table_panel: TablePanel::new(),
            spectrum_panel: SpectrumPanel::new(),
            session_log,
        }
    }

    pub(crate) fn handle_commands(&mut self) {
        let commands = std::mem::take(&mut self.pending_commands);

        for cmd in commands {
            tracing::debug!("Handling command: {:?}", cmd);

            if let Some(logger) = &mut self.session_log {
                if let Ok(json) = serde_json::to_string(&cmd) {
                    let _ = writeln!(logger, "{}", json);
                }
            }

            match cmd {
                AppCommand::RegenerateChromatogram => {
                    self.generate_chromatogram();
                }
                AppCommand::SelectElutionGroup(idx) => {
                    self.ui.selected_index = Some(idx);
                    self.pending_commands
                        .push(AppCommand::RegenerateChromatogram);
                }
                AppCommand::UpdateTolerance => {
                    if self.ui.selected_index.is_some() {
                        self.pending_commands
                            .push(AppCommand::RegenerateChromatogram);
                    }
                }
                AppCommand::UpdateSmoothing => {
                    if self.ui.selected_index.is_some() {
                        self.pending_commands
                            .push(AppCommand::RegenerateChromatogram);
                    }
                }
                AppCommand::QueryMS2Spectrum(rt_seconds) => {
                    if let Some(chrom_output) = &self.computed.chromatogram_output {
                        match chromatogram_processor::extract_ms2_spectrum_from_chromatogram(
                            chrom_output,
                            rt_seconds,
                        ) {
                            Ok(spectrum) => {
                                let num_peaks = spectrum.mz_values.len();
                                self.computed.ms2_spectrum = Some(spectrum);
                                tracing::info!(
                                    "Extracted MS2 spectrum at RT {:.2}s with {} peaks",
                                    rt_seconds,
                                    num_peaks
                                );
                            }
                            Err(e) => {
                                tracing::error!("Failed to extract MS2 spectrum: {:?}", e);
                                self.computed.ms2_spectrum = None;
                            }
                        }
                    } else {
                        tracing::warn!("No chromatogram data available for MS2 extraction");
                    }
                }
            }
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
                    self.ui.search_input.clear();
                }
                if i.key_pressed(egui::Key::Enter) {
                    self.ui.table_filter = self.ui.search_input.clone();
                    self.ui.search_mode = false;
                    self.ui.search_input.clear();
                }
                return;
            }

            if i.key_pressed(egui::Key::Slash) && !i.modifiers.any() {
                self.ui.search_mode = true;
                self.ui.search_input = self.ui.table_filter.clone();
                return;
            }

            let Some(egs) = &self.data.elution_groups else {
                return;
            };

            // Cache filtered indices - computed once and reused
            let filtered_indices = egs.matching_indices_for_id_filter(&self.ui.table_filter);

            if filtered_indices.is_empty() {
                return;
            }

            if i.key_pressed(egui::Key::J) && !i.modifiers.any() {
                self.move_selection_down(&filtered_indices);
            }

            if i.key_pressed(egui::Key::K) && !i.modifiers.any() {
                self.move_selection_up(&filtered_indices);
            }

            if i.key_pressed(egui::Key::G) && !i.modifiers.any() {
                self.move_selection_to_first(&filtered_indices);
            }

            if i.key_pressed(egui::Key::G) && i.modifiers.shift_only() {
                self.move_selection_to_last(&filtered_indices);
            }
        });
    }

    fn generate_chromatogram(&mut self) {
        let IndexedDataState::Loaded { index, ms1_rts } = &self.data.indexed_data else {
            return;
        };

        let selected_idx = match self.ui.selected_index {
            Some(idx) => idx,
            None => return,
        };

        if let Some(elution_groups) = &self.data.elution_groups {
            macro_rules! process_chromatogram {
                ($egs:expr) => {
                    match ChromatogramService::generate(
                        &$egs[selected_idx],
                        index,
                        Arc::clone(ms1_rts),
                        &self.data.tolerance,
                        &self.data.smoothing,
                    ) {
                        Ok(mut chrom) => {
                            // If we have library extras sidecar, attach library fragment intensities
                            if let Some(extras_map) = &elution_groups.extras {
                                if let Some(extras) = extras_map.get(&chrom.id) {
                                    let lib_ints: Vec<f32> = chrom
                                        .fragment_labels
                                        .iter()
                                        .map(|lbl| {
                                            extras.fragment_intensities
                                                .iter()
                                                .find(|(l, _)| l == lbl)
                                                .map(|(_, v)| *v)
                                                .unwrap_or(0.0)
                                        })
                                        .collect();
                                    chrom.library_fragment_intensities = Some(lib_ints);
                                }
                            }

                            let chrom_lines = ChromatogramLines::from_chromatogram(&chrom);
                            self.computed.chromatogram_x_bounds =
                                Some(chrom_lines.rt_seconds_range);
                            self.computed.chromatogram = Some(chrom_lines);
                            self.computed.chromatogram_output = Some(chrom);
                            // Reset auto zoom frame counter to 5 to force bounds update for a few frames
                            self.computed.auto_zoom_frame_counter = 5;
                        }
                        Err(e) => {
                            tracing::error!("Failed to generate chromatogram: {:?}", e);
                            self.computed.chromatogram = None;
                            self.computed.chromatogram_output = None;
                            self.computed.chromatogram_x_bounds = None;
                        }
                    }
                };
            }

            crate::with_elution_collection!(elution_groups, process_chromatogram);
        };

        if self.ui.show_ms2_spectrum
            && let Some(chrom_output) = &self.computed.chromatogram_output
        {
            let reference_rt = chrom_output.rt_seconds as f64;
            self.pending_commands
                .push(AppCommand::QueryMS2Spectrum(reference_rt));
        }
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
                    Ok((index, ms1_rts)) => {
                        data.indexed_data = IndexedDataState::Loaded { index, ms1_rts };
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

    fn move_selection_down(&mut self, filtered_indices: &[usize]) {
        if filtered_indices.is_empty() {
            return;
        }

        match self.ui.selected_index {
            None => {
                self.select_elution_group(filtered_indices[0]);
            }
            Some(current) => {
                if let Some(pos) = filtered_indices.iter().position(|&idx| idx == current)
                    && pos + 1 < filtered_indices.len()
                {
                    self.select_elution_group(filtered_indices[pos + 1]);
                }
            }
        }
    }

    fn move_selection_up(&mut self, filtered_indices: &[usize]) {
        if filtered_indices.is_empty() {
            return;
        }

        match self.ui.selected_index {
            None => {
                self.select_elution_group(filtered_indices[0]);
            }
            Some(current) => {
                if let Some(pos) = filtered_indices.iter().position(|&idx| idx == current)
                    && pos > 0
                {
                    self.select_elution_group(filtered_indices[pos - 1]);
                }
            }
        }
    }

    fn move_selection_to_first(&mut self, filtered_indices: &[usize]) {
        if !filtered_indices.is_empty() {
            self.select_elution_group(filtered_indices[0]);
        }
    }

    fn move_selection_to_last(&mut self, filtered_indices: &[usize]) {
        if !filtered_indices.is_empty() {
            self.select_elution_group(filtered_indices[filtered_indices.len() - 1]);
        }
    }

    fn select_elution_group(&mut self, idx: usize) {
        self.pending_commands
            .push(AppCommand::SelectElutionGroup(idx));
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
                show_ms2_spectrum: self.ui.show_ms2_spectrum,
            },
            tolerance: self.data.tolerance.clone(),
            smoothing: self.data.smoothing.clone(),
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
        self.handle_commands();

        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.heading("TimsQuery Viewer");
            ui.separator();
        });

        let mut tab_viewer = AppTabViewer {
            file_loader: &mut self.file_loader,
            data: &mut self.data,
            ui: &mut self.ui,
            computed: &mut self.computed,
            pending_commands: &mut self.pending_commands,
            left_panel: &mut self.left_panel,
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
    pending_commands: &'a mut Vec<AppCommand>,
    left_panel: &'a mut LeftPanel,
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

        let mut ctx = PanelContext::new(
            self.data,
            self.ui,
            self.computed,
            self.file_loader,
            self.pending_commands,
        );
        self.left_panel.render(ui, &mut ctx);
    }
}

impl<'a> TabViewer for AppTabViewer<'a> {
    type Tab = Pane;

    fn title(&mut self, tab: &mut Self::Tab) -> egui::WidgetText {
        match tab {
            Pane::LeftPanel => self.left_panel.title().into(),
            Pane::TablePanel => self.table_panel.title().into(),
            Pane::MS2Spectrum => self.spectrum_panel.title().into(),
            Pane::PrecursorPlot => "Precursors".into(),
            Pane::FragmentPlot => "Fragments".into(),
        }
    }

    fn ui(&mut self, ui: &mut egui::Ui, tab: &mut Self::Tab) {
        match tab {
            Pane::LeftPanel => {
                // Wrap settings in a scroll area
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        self.render_left_panel(ui);
                    });
            }
            Pane::TablePanel => {
                let mut ctx = PanelContext::new(
                    self.data,
                    self.ui,
                    self.computed,
                    self.file_loader,
                    self.pending_commands,
                );
                self.table_panel.render(ui, &mut ctx);
            }
            Pane::MS2Spectrum => {
                let mut ctx = PanelContext::new(
                    self.data,
                    self.ui,
                    self.computed,
                    self.file_loader,
                    self.pending_commands,
                );
                self.spectrum_panel.render(ui, &mut ctx);
            }
            Pane::PrecursorPlot => {
                if let Some(chromatogram) = &self.computed.chromatogram {
                    // Use shared link_id for synchronized X-axis with Fragments
                    let click_response = crate::plot_renderer::render_chromatogram_plot(
                        ui,
                        chromatogram,
                        crate::plot_renderer::PlotMode::PrecursorsOnly,
                        Some("precursor_fragment_x_axis"),
                        true,
                        &mut self.computed.auto_zoom_frame_counter,
                    );
                    if let Some(clicked_rt) = click_response {
                        self.pending_commands
                            .push(AppCommand::QueryMS2Spectrum(clicked_rt));
                    }
                } else if self.ui.selected_index.is_some() {
                    ui.label("Generating chromatogram...");
                } else {
                    ui.label("Select a precursor from the table to view chromatogram");
                }
            }
            Pane::FragmentPlot => {
                if let Some(chromatogram) = &self.computed.chromatogram {
                    // Use shared link_id for synchronized X-axis with Precursors
                    let response = crate::plot_renderer::render_chromatogram_plot(
                        ui,
                        chromatogram,
                        crate::plot_renderer::PlotMode::FragmentsOnly,
                        Some("precursor_fragment_x_axis"),
                        false,
                        &mut self.computed.auto_zoom_frame_counter,
                    );
                    if let Some(clicked_rt) = response {
                        self.pending_commands
                            .push(AppCommand::QueryMS2Spectrum(clicked_rt));
                    }
                } else if self.ui.selected_index.is_some() {
                    ui.label("Generating chromatogram...");
                } else {
                    ui.label("Select a precursor to view fragment traces");
                }
            }
        }
    }

    fn is_closeable(&self, _tab: &Self::Tab) -> bool {
        false // Tabs cannot be closed
    }
}

#[cfg(test)]
mod tests;
