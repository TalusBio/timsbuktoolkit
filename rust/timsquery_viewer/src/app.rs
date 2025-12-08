use eframe::egui;
use egui_dock::{
    DockArea,
    DockState,
    Style,
    TabViewer,
};
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::KeyLike;
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::Tolerance;

use crate::chromatogram_processor::{
    self,
    ChromatogramOutput,
    SmoothingMethod,
};
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

/// Pane types for the tile layout
#[derive(Debug)]
enum Pane {
    LeftPanel,
    TablePanel,
    MS2Spectrum,
    PrecursorPlot,
    FragmentPlot,
}

/// TabViewer implementation for the dock layout
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

/// Commands that trigger state changes in the application
#[derive(Debug, Clone)]
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

/// Domain/data state - represents loaded data and analysis parameters
#[derive(Debug, Default)]
pub struct DataState {
    /// Loaded elution groups
    pub elution_groups: Option<ElutionGroupData>,
    /// Loaded and indexed timsTOF data
    pub indexed_data: Option<Arc<IndexedTimstofPeaks>>,
    /// MS1 retention times (in milliseconds)
    pub ms1_rts: Option<Arc<[u32]>>,
    /// Tolerance settings
    pub tolerance: Tolerance,
    /// Smoothing method configuration
    pub smoothing: SmoothingMethod,
}

/// UI-specific state - transient UI state that doesn't affect data
#[derive(Debug)]
pub struct UiState {
    /// Filter text for precursor table
    pub table_filter: String,
    /// Currently selected elution group index
    pub selected_index: Option<usize>,
    /// Vim mode: search mode active
    pub search_mode: bool,
    /// Vim mode: search input buffer
    pub search_input: String,
    /// Track if we're currently loading/indexing data
    pub is_indexing: bool,
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
            is_indexing: false,
            show_ms2_spectrum: true, // Enable MS2 by default
        }
    }
}

/// Computed/cached state - derived from data and UI state
#[derive(Debug, Default)]
pub struct ComputedState {
    /// Computed chromatogram for the selected elution group (plot data)
    pub chromatogram: Option<ChromatogramLines>,
    /// Whether auto-zoom has been applied to the current chromatogram
    pub chromatogram_auto_zoom_applied: bool,
    /// Raw chromatogram output data (for MS2 extraction)
    pub chromatogram_output: Option<ChromatogramOutput>,
    /// Computed MS2 spectrum at selected RT
    pub ms2_spectrum: Option<MS2Spectrum>,
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
}

impl ViewerApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
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
        }
    }
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

        let commands = self.left_panel.render_tolerance_editor(ui, self.data);
        self.pending_commands.extend(commands);

        ui.add_space(20.0);
        ui.separator();

        let commands = self.left_panel.render_smoothing(ui, self.data);
        self.pending_commands.extend(commands);

        ui.add_space(20.0);
        ui.separator();

        // Plot visibility controls
        ui.heading("Plot Settings");
        ui.checkbox(&mut self.ui.show_ms2_spectrum, "Show MS2 Spectrum");
    }
}

impl<'a> TabViewer for AppTabViewer<'a> {
    type Tab = Pane;

    fn title(&mut self, tab: &mut Self::Tab) -> egui::WidgetText {
        match tab {
            Pane::LeftPanel => "Settings".into(),
            Pane::TablePanel => "Table".into(),
            Pane::MS2Spectrum => "MS2".into(),
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
                if let Some(egs) = &self.data.elution_groups {
                    if let Some(cmd) = self.table_panel.render(ui, egs, self.ui) {
                        self.pending_commands.push(cmd);
                    }
                } else {
                    ui.heading("Precursor Table");
                    ui.separator();
                    ui.label("Load elution groups to see the table");
                }
            }
            Pane::MS2Spectrum => {
                self.spectrum_panel
                    .render(ui, self.computed.ms2_spectrum.as_ref());
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
                        !self.computed.chromatogram_auto_zoom_applied,
                    );
                    self.computed.chromatogram_auto_zoom_applied = true;
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

impl eframe::App for ViewerApp {
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

impl ViewerApp {
    /// Process all pending commands
    fn handle_commands(&mut self) {
        let commands = std::mem::take(&mut self.pending_commands);

        for cmd in commands {
            tracing::debug!("Handling command: {:?}", cmd);

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
        ui_state: &mut UiState,
    ) {
        ui.label("Raw Data File (.d):");
        if ui.button("Load Raw Data...").clicked() {
            file_loader.open_raw_data_dialog();
        }

        if let Some(path) = &file_loader.raw_data_path {
            Self::display_filename(ui, path);
        }

        Self::load_raw_data_if_needed(ui, file_loader, data, ui_state);

        if data.indexed_data.is_some() && !ui_state.is_indexing {
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
        if let Some(path) = &file_loader.elution_groups_path
            && data.elution_groups.is_none()
        {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label("Loading elution groups...");
            });
            match file_loader.load_elution_groups(path) {
                Ok(egs) => {
                    tracing::info!("Loaded {} elution groups", egs.len());
                    data.elution_groups = Some(egs);
                }
                Err(e) => {
                    tracing::error!("Failed to load elution groups: {:?}", e);
                }
            }
        }
    }

    fn load_raw_data_if_needed(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
        ui_state: &mut UiState,
    ) {
        if let Some(path) = &file_loader.raw_data_path {
            if data.indexed_data.is_none() && !ui_state.is_indexing {
                ui_state.is_indexing = true;
                ui.ctx().request_repaint();
            }

            if ui_state.is_indexing {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Indexing raw data... (this may take 10-30 seconds)");
                });

                match file_loader.load_raw_data(path) {
                    Ok((index, rts)) => {
                        data.indexed_data = Some(index);
                        data.ms1_rts = Some(rts);
                        ui_state.is_indexing = false;
                        tracing::info!("Raw data indexing completed");
                    }
                    Err(e) => {
                        tracing::error!("Failed to load raw data: {:?}", e);
                        ui_state.is_indexing = false;
                    }
                }
            }
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

    fn generate_chromatogram(&mut self) {
        let Some(index) = &self.data.indexed_data else {
            return;
        };
        let Some(ms1_rts) = &self.data.ms1_rts else {
            return;
        };

        let selected_idx = match self.ui.selected_index {
            Some(idx) => idx,
            None => return,
        };

        macro_rules! process_chromatogram {
            ($egs:expr) => {
                Self::process_chromatogram(
                    &mut self.computed.chromatogram,
                    &mut self.computed.chromatogram_output,
                    &$egs[selected_idx],
                    index,
                    ms1_rts,
                    &self.data.tolerance,
                    &self.data.smoothing,
                );
            };
        }

        if let Some(x) = &self.data.elution_groups {
            // This feels pretty dirty ... I am pretty sure I can implement a getter
            // that just casts to the lowest common denominator trait object or something
            match &x.inner {
                timsquery::serde::ElutionGroupCollection::StringLabels(egs) => {
                    process_chromatogram!(egs);
                }
                timsquery::serde::ElutionGroupCollection::MzpafLabels(egs) => {
                    process_chromatogram!(egs);
                }
                timsquery::serde::ElutionGroupCollection::TinyIntLabels(egs) => {
                    process_chromatogram!(egs);
                }
                timsquery::serde::ElutionGroupCollection::IntLabels(egs) => {
                    process_chromatogram!(egs);
                }
            }
        };
        self.computed.chromatogram_auto_zoom_applied = false;

        if self.ui.show_ms2_spectrum
            && let Some(chrom_output) = &self.computed.chromatogram_output
        {
            let reference_rt = chrom_output.rt_seconds as f64;
            self.pending_commands
                .push(AppCommand::QueryMS2Spectrum(reference_rt));
        }
    }

    fn process_chromatogram<T: KeyLike + std::fmt::Display>(
        chromatogram_lines: &mut Option<ChromatogramLines>,
        chromatogram_output: &mut Option<ChromatogramOutput>,
        eg: &TimsElutionGroup<T>,
        index: &Arc<IndexedTimstofPeaks>,
        ms1_rts: &Arc<[u32]>,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) {
        match chromatogram_processor::generate_chromatogram(
            eg,
            index,
            Arc::clone(ms1_rts),
            tolerance,
            smoothing,
        ) {
            Ok(chrom) => {
                tracing::info!(
                    "Generated chromatogram for elution group {} with {} precursors and {} fragments",
                    chrom.id,
                    chrom.precursor_mzs.len(),
                    chrom.fragment_mzs.len()
                );
                *chromatogram_lines = Some(ChromatogramLines::from_chromatogram(&chrom));
                *chromatogram_output = Some(chrom);
            }
            Err(e) => {
                tracing::error!("Failed to generate chromatogram: {:?}", e);
                *chromatogram_lines = None;
                *chromatogram_output = None;
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
