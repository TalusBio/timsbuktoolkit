use eframe::egui;
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::KeyLike;
use timsquery::models::elution_group::ElutionGroup;
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
use crate::panels::{
    LeftPanel,
    PlotPanel,
    SpectrumPanel,
    TablePanel,
};
use crate::plot_renderer::{ChromatogramLines, MS2Spectrum};

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
#[derive(Debug)]
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

impl Default for DataState {
    fn default() -> Self {
        Self {
            elution_groups: None,
            indexed_data: None,
            ms1_rts: None,
            tolerance: Tolerance::default(),
            smoothing: SmoothingMethod::default(),
        }
    }
}

/// UI-specific state - transient UI state that doesn't affect data
#[derive(Debug, Default)]
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
    /// Toggle for precursor/fragment split view
    pub show_split_xic: bool,
}

/// Computed/cached state - derived from data and UI state
#[derive(Debug, Default)]
pub struct ComputedState {
    /// Computed chromatogram for the selected elution group (plot data)
    pub chromatogram: Option<ChromatogramLines>,
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

    /// UI Panels
    left_panel: LeftPanel,
    table_panel: TablePanel,
    plot_panel: PlotPanel,
    spectrum_panel: SpectrumPanel,
}

impl ViewerApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            file_loader: FileLoader::new(),
            data: DataState::default(),
            ui: UiState::default(),
            computed: ComputedState::default(),
            pending_commands: Vec::new(),
            left_panel: LeftPanel::new(),
            table_panel: TablePanel::new(),
            plot_panel: PlotPanel::new(),
            spectrum_panel: SpectrumPanel::new(),
        }
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

        egui::SidePanel::left("left_panel")
            .default_width(300.0)
            .max_width(350.0)
            .resizable(true)
            .show(ctx, |ui| {
                self.render_left_panel(ui);
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            if let Some(egs) = &self.data.elution_groups {
                if let Some(cmd) = self.table_panel.render(ui, egs, &mut self.ui) {
                    self.pending_commands.push(cmd);
                }
            } else {
                ui.heading("Precursor Table");
                ui.separator();
                ui.label("Load elution groups to see the table");
            }
        });

        egui::SidePanel::right("right_panel")
            .default_width(600.0)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.checkbox(&mut self.ui.show_ms2_spectrum, "Show MS2 Spectrum");
                    ui.checkbox(&mut self.ui.show_split_xic, "Split XIC View");
                });
                ui.separator();

                if self.ui.show_ms2_spectrum {
                    let available_height = ui.available_height();
                    let ms2_height = available_height * 0.3;
                    let xic_height = available_height * 0.7;

                    ui.allocate_ui_with_layout(
                        egui::vec2(ui.available_width(), ms2_height),
                        egui::Layout::top_down(egui::Align::LEFT),
                        |ui| {
                            self.spectrum_panel
                                .render(ui, self.computed.ms2_spectrum.as_ref());
                        },
                    );

                    ui.separator();

                    ui.allocate_ui_with_layout(
                        egui::vec2(ui.available_width(), xic_height),
                        egui::Layout::top_down(egui::Align::LEFT),
                        |ui| {
                            if let Some(clicked_rt) = self.plot_panel.render(
                                ui,
                                &self.computed,
                                self.ui.selected_index,
                                self.ui.show_split_xic,
                            ) {
                                self.pending_commands
                                    .push(AppCommand::QueryMS2Spectrum(clicked_rt));
                            }
                        },
                    );
                } else {
                    if let Some(clicked_rt) = self.plot_panel.render(
                        ui,
                        &self.computed,
                        self.ui.selected_index,
                        self.ui.show_split_xic,
                    ) {
                        self.pending_commands
                            .push(AppCommand::QueryMS2Spectrum(clicked_rt));
                    }
                }
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
                                self.computed.ms2_spectrum = Some(spectrum);
                                tracing::info!(
                                    "Extracted MS2 spectrum at RT {:.2}s with {} peaks",
                                    rt_seconds,
                                    self.computed.ms2_spectrum.as_ref().unwrap().mz_values.len()
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

            let filtered_indices = egs.filter_by_id(&self.ui.table_filter);

            if filtered_indices.is_empty() {
                return;
            }

            if i.key_pressed(egui::Key::J) && !i.modifiers.any() {
                match self.ui.selected_index {
                    None => {
                        self.pending_commands
                            .push(AppCommand::SelectElutionGroup(filtered_indices[0]));
                    }
                    Some(current) => {
                        if let Some(pos) = filtered_indices.iter().position(|&idx| idx == current) {
                            if pos + 1 < filtered_indices.len() {
                                self.pending_commands.push(AppCommand::SelectElutionGroup(
                                    filtered_indices[pos + 1],
                                ));
                            }
                        }
                    }
                }
            }

            if i.key_pressed(egui::Key::K) && !i.modifiers.any() {
                match self.ui.selected_index {
                    None => {
                        self.pending_commands
                            .push(AppCommand::SelectElutionGroup(filtered_indices[0]));
                    }
                    Some(current) => {
                        if let Some(pos) = filtered_indices.iter().position(|&idx| idx == current) {
                            if pos > 0 {
                                self.pending_commands.push(AppCommand::SelectElutionGroup(
                                    filtered_indices[pos - 1],
                                ));
                            }
                        }
                    }
                }
            }

            if i.key_pressed(egui::Key::G) && !i.modifiers.any() {
                self.pending_commands
                    .push(AppCommand::SelectElutionGroup(filtered_indices[0]));
            }

            if i.key_pressed(egui::Key::G) && i.modifiers.shift_only() {
                self.pending_commands.push(AppCommand::SelectElutionGroup(
                    filtered_indices[filtered_indices.len() - 1],
                ));
            }
        });
    }

    fn render_left_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Data Loading");
        ui.separator();

        self.render_elution_groups_section(ui);
        ui.add_space(10.0);

        self.render_raw_data_section(ui);
        ui.add_space(10.0);

        self.render_tolerance_loading_section(ui);

        ui.add_space(20.0);
        ui.separator();

        let commands = self.left_panel.render_tolerance_editor(ui, &mut self.data);
        self.pending_commands.extend(commands);

        ui.add_space(20.0);
        ui.separator();

        let commands = self.left_panel.render_smoothing(ui, &mut self.data);
        self.pending_commands.extend(commands);
    }

    fn render_elution_groups_section(&mut self, ui: &mut egui::Ui) {
        ui.label("Elution Groups:");
        if ui.button("Load Elution Groups...").clicked() {
            self.file_loader.open_elution_groups_dialog();
        }

        if let Some(path) = &self.file_loader.elution_groups_path {
            Self::display_filename(ui, path);
        }

        Self::load_elution_groups_if_needed(ui, &mut self.file_loader, &mut self.data);

        if let Some(egs) = &self.data.elution_groups {
            ui.label(format!("✓ Loaded: {} elution groups", egs.len()));
        }
    }

    fn render_raw_data_section(&mut self, ui: &mut egui::Ui) {
        ui.label("Raw Data File (.d):");
        if ui.button("Load Raw Data...").clicked() {
            self.file_loader.open_raw_data_dialog();
        }

        if let Some(path) = &self.file_loader.raw_data_path {
            Self::display_filename(ui, path);
        }

        Self::load_raw_data_if_needed(ui, &mut self.file_loader, &mut self.data, &mut self.ui);

        if self.data.indexed_data.is_some() && !self.ui.is_indexing {
            ui.label("✓ Raw data indexed");
        }
    }

    fn render_tolerance_loading_section(&mut self, ui: &mut egui::Ui) {
        ui.label("Tolerance Settings:");
        if ui.button("Load Tolerances...").clicked() {
            self.file_loader.open_tolerance_dialog();
        }

        if let Some(path) = &self.file_loader.tolerance_path {
            Self::display_filename(ui, path);
        }

        Self::load_tolerance_if_needed(&mut self.file_loader, &mut self.data);

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
        let Some(index) = self.data.indexed_data.clone() else {
            return;
        };
        let Some(ms1_rts) = self.data.ms1_rts.clone() else {
            return;
        };

        let selected_idx = match self.ui.selected_index {
            Some(idx) => idx,
            None => return,
        };

        match &self.data.elution_groups {
            Some(ElutionGroupData::StringLabels(egs)) => Self::process_chromatogram(
                &mut self.computed.chromatogram,
                &mut self.computed.chromatogram_output,
                &egs[selected_idx],
                &index,
                ms1_rts,
                &self.data.tolerance,
                &self.data.smoothing,
            ),
            Some(ElutionGroupData::MzpafLabels(egs)) => Self::process_chromatogram(
                &mut self.computed.chromatogram,
                &mut self.computed.chromatogram_output,
                &egs[selected_idx],
                &index,
                ms1_rts,
                &self.data.tolerance,
                &self.data.smoothing,
            ),
            None => (),
        };

        if self.ui.show_ms2_spectrum {
            if let Some(chrom_output) = &self.computed.chromatogram_output {
                let reference_rt = chrom_output.rt_seconds as f64;
                self.pending_commands
                    .push(AppCommand::QueryMS2Spectrum(reference_rt));
            }
        }
    }

    fn process_chromatogram<T: KeyLike + std::fmt::Display>(
        chromatogram_lines: &mut Option<ChromatogramLines>,
        chromatogram_output: &mut Option<ChromatogramOutput>,
        eg: &ElutionGroup<T>,
        index: &Arc<IndexedTimstofPeaks>,
        ms1_rts: Arc<[u32]>,
        tolerance: &Tolerance,
        smoothing: &SmoothingMethod,
    ) {
        match chromatogram_processor::generate_chromatogram(
            eg, index, ms1_rts, tolerance, smoothing,
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
}
