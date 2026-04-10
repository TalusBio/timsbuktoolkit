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
use std::time::Instant;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::IndexedPeaksHandle;

use crate::calibration::ViewerCalibrationState;
use crate::chromatogram_processor::SmoothingMethod;
use crate::cli::Cli;
use crate::computed_state::{
    ChromatogramComputationResult,
    ComputedState,
};
use crate::file_loader::{
    ElutionGroupData,
    FileLoader,
};
use crate::plot_renderer::AutoZoomMode;
use crate::ui::panels::{
    ConfigPanel,
    MobilityPanel,
    ScreenshotAction,
    ScreenshotState,
    SpectrumPanel,
    TablePanel,
};
use std::sync::atomic::{
    AtomicBool,
    Ordering,
};
use std::sync::mpsc::{
    Receiver,
    channel,
};

/// Result from background chromatogram computation
type ChromatogramComputeResult = Result<
    (
        crate::chromatogram_processor::ChromatogramOutput,
        crate::chromatogram_processor::ChromatogramCollector<String, f32>,
        timsseek::ExpectedIntensities<String>,
        u64, // selected_idx as cache key
        timsquery::models::elution_group::TimsElutionGroup<String>,
    ),
    String,
>;

/// Handle to cancel a running background computation
#[derive(Clone)]
struct CancellationToken {
    cancelled: Arc<AtomicBool>,
}

impl CancellationToken {
    fn new() -> Self {
        Self {
            cancelled: Arc::new(AtomicBool::new(false)),
        }
    }

    fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
    }

    fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }
}

// UI spacing and layout constants
const SECTION_MARGIN: i8 = 10;
const SECTION_SPACING: f32 = 12.0;
const INTERNAL_SPACING: f32 = 8.0;
const SMALL_SPACING: f32 = 4.0;
const SMALL_SPACING_VERTICAL: f32 = 6.0;
const SEPARATOR_SPACING: f32 = 16.0;
const LOCATION_DISPLAY_MAX_LEN: usize = 40;
const CORNER_RADIUS: u8 = 4;
const LOCATION_FRAME_PADDING_H: i8 = 8;
const LOCATION_FRAME_PADDING_V: i8 = 4;

/// Pane types for the tile layout
#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize, serde::Serialize)]
enum Pane {
    ConfigPanel,
    TablePanel,
    MS2Spectrum,
    PrecursorPlot,
    FragmentPlot,
    ScoresPlot,
    Mobility,
    Calibration,
}

/// All pane variants that should exist in the dock.
/// New variants added here will be auto-injected into saved layouts.
const ALL_PANES: &[Pane] = &[
    Pane::ConfigPanel,
    Pane::TablePanel,
    Pane::MS2Spectrum,
    Pane::PrecursorPlot,
    Pane::FragmentPlot,
    Pane::ScoresPlot,
    Pane::Mobility,
    Pane::Calibration,
];

/// State to be persisted across restarts
#[derive(serde::Deserialize, serde::Serialize)]
struct PersistentState {
    file_loader: FileLoader,
    ui_state: UiState,
    tolerance: Tolerance,
    smoothing: SmoothingMethod,
    dock_state: DockState<Pane>,
}

fn default_true() -> bool {
    true
}

/// State of indexed raw data loading
#[derive(Debug, Default)]
pub enum IndexedDataState {
    /// No data loaded
    #[default]
    None,
    /// Load requested but not yet started (triggers load on next frame)
    LoadRequested(String),
    /// Load in progress - shows spinner, does NOT re-trigger load
    Loading(String),
    /// Loading failed with error message
    Failed(String, String),
    /// Data successfully loaded
    Loaded {
        index: Arc<IndexedPeaksHandle>,
        source: String,
    },
}

/// State of elution group library loading
#[derive(Debug, Default)]
pub enum ElutionGroupState {
    /// No library loaded
    #[default]
    None,
    /// Load requested but not yet started (triggers load on next frame)
    LoadRequested(PathBuf),
    /// Load in progress - shows spinner, does NOT re-trigger load
    Loading(PathBuf),
    /// Loading failed with error message (path, error)
    Failed(PathBuf, String),
    /// Library successfully loaded
    Loaded {
        data: Arc<ElutionGroupData>,
        source: PathBuf,
    },
}

impl ElutionGroupState {
    pub fn as_ref(&self) -> Option<&ElutionGroupData> {
        match self {
            ElutionGroupState::Loaded { data, .. } => Some(data),
            _ => None,
        }
    }

    pub fn is_none(&self) -> bool {
        matches!(self, ElutionGroupState::None)
    }
}

/// State of tolerance file loading
#[derive(Debug, Default)]
pub enum ToleranceState {
    /// No tolerance file loaded (using defaults or user-edited values)
    #[default]
    None,
    /// Load requested but not yet started (triggers load on next frame)
    LoadRequested(PathBuf),
    /// Tolerance file successfully loaded (user can freely edit via UI)
    Loaded { source: PathBuf },
    /// Loading failed with error message (path, error)
    Failed(PathBuf, String),
}

/// Domain/data state - represents loaded data and analysis parameters
#[derive(Debug, Default)]
pub struct DataState {
    /// Elution group library with loading state
    pub elution_groups: ElutionGroupState,
    /// Indexed timsTOF data with loading state
    pub indexed_data: IndexedDataState,
    /// Tolerance file loading state
    pub tolerance_state: ToleranceState,
    /// Tolerance settings (live-editable values)
    pub tolerance: Tolerance,
    /// Smoothing method configuration
    pub smoothing: SmoothingMethod,
    /// Auto-zoom mode for plots
    pub auto_zoom_mode: AutoZoomMode,
}

/// Input mode for raw data loading
#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize, serde::Serialize, Default)]
pub enum RawDataInputMode {
    #[default]
    Local,
    Cloud,
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
    /// Raw data input mode (Local or Cloud)
    pub raw_data_input_mode: RawDataInputMode,
    /// Dark mode preference (true = dark, false = light)
    #[serde(default = "default_true")]
    pub dark_mode: bool,
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
    mobility_panel: MobilityPanel,

    /// Receiver for background chromatogram computation
    chromatogram_receiver: Option<Receiver<ChromatogramComputeResult>>,
    /// Token to cancel current background computation
    cancellation_token: Option<CancellationToken>,

    /// Screenshot capture state machine
    screenshot_state: ScreenshotState,
    /// Countdown duration in seconds before capture
    screenshot_delay_secs: f32,

    /// Calibration state machine (background scoring + calibrt)
    calibration: ViewerCalibrationState,
}

fn apply_theme(ctx: &egui::Context, dark_mode: bool) {
    ctx.set_visuals(if dark_mode {
        egui::Visuals::dark()
    } else {
        egui::Visuals::light()
    });
}

impl ViewerApp {
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
                    let mut dock_state = state.dock_state;
                    // Inject any pane variants added since the state was saved
                    for &pane in ALL_PANES {
                        if dock_state.find_tab(&pane).is_none() {
                            tracing::info!("Adding missing pane {:?} to saved layout", pane);
                            dock_state.main_surface_mut().push_to_first_leaf(pane);
                        }
                    }

                    apply_theme(&cc.egui_ctx, state.ui_state.dark_mode);

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
                        dock_state,
                        config_panel: ConfigPanel::new(),
                        table_panel: TablePanel::new(),
                        spectrum_panel: SpectrumPanel::new(),
                        mobility_panel: MobilityPanel::new(),
                        chromatogram_receiver: None,
                        cancellation_token: None,
                        screenshot_state: ScreenshotState::default(),
                        screenshot_delay_secs: 3.0,
                        calibration: ViewerCalibrationState::default(),
                    };
                }
            } else {
                tracing::info!("No persistent state found.");
            }
        }

        // Create initial tabs: Settings, Table, Precursors, Fragments, MS2, Scores, Mobilogram
        let tabs = vec![
            Pane::ConfigPanel,
            Pane::TablePanel,
            Pane::PrecursorPlot,
            Pane::FragmentPlot,
            Pane::MS2Spectrum,
            Pane::ScoresPlot,
            Pane::Mobility,
            Pane::Calibration,
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
            mobility_panel: MobilityPanel::new(),
            chromatogram_receiver: None,
            cancellation_token: None,
            screenshot_state: ScreenshotState::default(),
            screenshot_delay_secs: 3.0,
            calibration: ViewerCalibrationState::default(),
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

    /// Respond to a user RT click: generate MS2 spectrum and mobility data.
    fn handle_rt_click(&mut self) {
        let Some(requested_rt) = self.computed.clicked_rt() else {
            return;
        };

        if self.computed.generate_spectrum_at_rt(requested_rt) {
            self.computed
                .insert_reference_line("Clicked RT".into(), requested_rt, Color32::GREEN);
        }

        if let IndexedDataState::Loaded { index, .. } = &self.data.indexed_data {
            self.computed
                .generate_mobility_at_rt(requested_rt, index, &self.data.tolerance);
        }
    }

    fn generate_chromatogram(&mut self, ctx: &egui::Context) {
        let IndexedDataState::Loaded { index, .. } = &self.data.indexed_data else {
            return;
        };

        let selected_idx = match self.ui.selected_index {
            Some(idx) => idx,
            None => return,
        };

        let Some(elution_groups) = self.data.elution_groups.as_ref() else {
            return;
        };

        if !self
            .computed
            .is_cache_valid(selected_idx, &self.data.tolerance, &self.data.smoothing)
        {
            // Cache miss - chromatogram will be regenerated
            tracing::debug!(
                "Cache MISS for index {} - triggering chromatogram regeneration",
                selected_idx
            );

            let is_new_request = self
                .computed
                .computing_index()
                .map(|computing_idx| computing_idx != selected_idx as u64)
                .unwrap_or(true);

            if self.computed.is_computing() && !is_new_request {
                tracing::trace!(
                    "Chromatogram computation already in progress for index {}, skipping",
                    selected_idx
                );
                return;
            }

            let Ok((elution_group, expected_intensities)) = elution_groups.get_elem(selected_idx)
            else {
                tracing::error!("Invalid elution group index: {}", selected_idx);
                return;
            };

            if self.computed.is_computing() && is_new_request {
                tracing::debug!(
                    "Cancelling previous chromatogram computation (switching to index {})",
                    selected_idx
                );
                if let Some(token) = &self.cancellation_token {
                    token.cancel();
                }
                self.chromatogram_receiver = None;
                self.cancellation_token = None;
            }

            tracing::debug!(
                "Starting new chromatogram computation for index {}",
                selected_idx
            );
            self.computed.start_computing(selected_idx as u64);

            let cancel_token = CancellationToken::new();
            self.cancellation_token = Some(cancel_token.clone());

            let (tx, rx) = channel();
            self.chromatogram_receiver = Some(rx);

            let index_owned = index.clone();
            let elution_group_owned = elution_group.clone();
            let expected_intensities_owned = expected_intensities.clone();
            let tolerance_owned = self.data.tolerance.clone();
            let smoothing_owned = self.data.smoothing;
            let ctx_clone = ctx.clone();
            match std::thread::Builder::new()
                .name(format!("chrom-{}", selected_idx))
                .spawn(move || {
                    tracing::debug!(
                        "Starting chromatogram computation for elution group {}",
                        selected_idx
                    );
                    let result = Self::compute_chromatogram_background(
                        elution_group_owned,
                        expected_intensities_owned,
                        selected_idx,
                        index_owned,
                        tolerance_owned,
                        smoothing_owned,
                        cancel_token,
                    );
                    if let Err(e) = &result {
                        tracing::error!("Chromatogram computation failed: {}", e);
                    } else {
                        tracing::debug!(
                            "Chromatogram computation completed for elution group {}",
                            selected_idx
                        );
                    }
                    let _ = tx.send(result);
                    ctx_clone.request_repaint();
                }) {
                Ok(_) => {
                    tracing::debug!("Successfully spawned chromatogram computation thread");
                }
                Err(e) => {
                    tracing::error!("Failed to spawn chromatogram thread: {}", e);
                    let index = self.computed.computing_index().unwrap_or(0);
                    self.computed.fail_computing(
                        index,
                        &self.data.tolerance,
                        &self.data.smoothing,
                        format!("Failed to spawn computation thread: {}", e),
                    );
                    self.chromatogram_receiver = None;
                    self.cancellation_token = None;
                }
            }
        }
    }

    /// Compute chromatogram in background thread
    fn compute_chromatogram_background(
        elution_group: timsquery::models::elution_group::TimsElutionGroup<String>,
        expected_intensities: timsseek::ExpectedIntensities<String>,
        selected_idx: usize,
        index: Arc<IndexedPeaksHandle>,
        tolerance: Tolerance,
        smoothing: SmoothingMethod,
        cancel_token: CancellationToken,
    ) -> ChromatogramComputeResult {
        // Check if cancelled before starting
        if cancel_token.is_cancelled() {
            tracing::debug!("Chromatogram computation cancelled before starting");
            return Err("Computation cancelled".to_string());
        }

        // Build collector
        let mut collector = ComputedState::build_collector(&index, elution_group.clone())
            .map_err(|e| format!("Failed to build collector: {:?}", e))?;

        // Check if cancelled after building collector
        if cancel_token.is_cancelled() {
            tracing::debug!("Chromatogram computation cancelled after building collector");
            return Err("Computation cancelled".to_string());
        }

        // Generate chromatogram
        let output = ComputedState::generate_chromatogram(
            &mut collector,
            &elution_group,
            &index,
            &tolerance,
            &smoothing,
        )
        .map_err(|e| format!("Failed to generate chromatogram: {:?}", e))?;

        // Check if cancelled after generating chromatogram
        if cancel_token.is_cancelled() {
            tracing::debug!("Chromatogram computation cancelled after generation");
            return Err("Computation cancelled".to_string());
        }

        Ok((
            output,
            collector,
            expected_intensities,
            selected_idx as u64,
            elution_group,
        ))
    }

    /// Check if background chromatogram computation completed
    fn check_chromatogram_completion(&mut self) {
        let Some(rx) = &self.chromatogram_receiver else {
            return;
        };

        // Try to receive result (non-blocking)
        if let Ok(result) = rx.try_recv() {
            self.chromatogram_receiver = None;
            self.cancellation_token = None;

            match result {
                Ok((output, collector, expected_intensities, selected_idx, elution_group)) => {
                    tracing::debug!(
                        "Chromatogram computation result received for index {}",
                        selected_idx
                    );

                    if let IndexedDataState::Loaded { index, .. } = &self.data.indexed_data {
                        let result = ChromatogramComputationResult {
                            selected_idx,
                            output,
                            collector,
                            expected_intensities,
                            elution_group,
                        };

                        self.computed.complete_chromatogram_computation(
                            result,
                            index,
                            &self.data.tolerance,
                            self.data.smoothing,
                        );
                    }

                    // Add library RT reference line
                    if let Some(elution_groups) = self.data.elution_groups.as_ref()
                        && let Ok((elution_group, _)) =
                            elution_groups.get_elem(selected_idx as usize)
                    {
                        self.computed.insert_reference_line(
                            "Library RT".into(),
                            elution_group.rt_seconds() as f64,
                            Color32::BLUE,
                        );
                    }
                }
                Err(e) => {
                    tracing::error!("Chromatogram computation failed: {}", e);
                    // Capture the computing index before clearing state
                    let index = self.computed.computing_index().unwrap_or(0);
                    self.computed.fail_computing(
                        index,
                        &self.data.tolerance,
                        &self.data.smoothing,
                        e,
                    );
                }
            }
        }
    }

    fn render_elution_groups_section_static(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
        ui_state: &mut UiState,
    ) {
        egui::Frame::group(ui.style())
            .inner_margin(egui::Margin::same(SECTION_MARGIN))
            .show(ui, |ui| {
                ui.heading("Elution Groups");
                ui.add_space(INTERNAL_SPACING);

                if ui.button("Load Elution Groups...").clicked() {
                    file_loader.open_elution_groups_dialog();
                }

                if let Some(path) = &file_loader.elution_groups_path {
                    ui.add_space(SMALL_SPACING);
                    Self::display_filename(ui, path);
                }

                Self::load_elution_groups_if_needed(ui, file_loader, data, ui_state);

                if let Some(egs) = data.elution_groups.as_ref() {
                    ui.add_space(SMALL_SPACING);
                    ui.label(
                        egui::RichText::new(format!("✓ {} groups loaded", egs.len()))
                            .color(egui::Color32::DARK_GREEN),
                    );
                }
            });
    }

    fn render_raw_data_section_static(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
        ui_state: &mut UiState,
        computed: &mut ComputedState,
    ) {
        egui::Frame::group(ui.style())
            .inner_margin(egui::Margin::same(SECTION_MARGIN))
            .show(ui, |ui| {
                ui.heading("Raw Data");
                ui.add_space(INTERNAL_SPACING);

                // Clearer tab selector with better visual style
                ui.horizontal(|ui| {
                    ui.label("Source:");
                    ui.radio_value(
                        &mut ui_state.raw_data_input_mode,
                        RawDataInputMode::Local,
                        "Local File",
                    );
                    ui.radio_value(
                        &mut ui_state.raw_data_input_mode,
                        RawDataInputMode::Cloud,
                        "Cloud URL",
                    );
                });

                ui.add_space(INTERNAL_SPACING);

                match ui_state.raw_data_input_mode {
                    RawDataInputMode::Local => {
                        // Existing folder picker button
                        if ui.button("Browse for Raw Data...").clicked() {
                            file_loader.open_raw_data_dialog();
                        }
                        ui.label(
                            egui::RichText::new("Select a .d folder or .idx cache")
                                .small()
                                .weak(),
                        );
                    }
                    RawDataInputMode::Cloud => {
                        // URL input field
                        ui.vertical(|ui| {
                            ui.label("Enter cloud storage URL:");
                            let mut url_buffer =
                                file_loader.raw_data_url.clone().unwrap_or_default();
                            let response = ui.add(
                                egui::TextEdit::singleline(&mut url_buffer)
                                    .hint_text("s3://bucket/data.idx or gs://bucket/data.idx")
                                    .desired_width(f32::INFINITY),
                            );

                            if response.changed() {
                                if !url_buffer.is_empty() {
                                    file_loader.set_raw_data_url(url_buffer);
                                } else {
                                    file_loader.clear_raw_data();
                                }
                            }

                            ui.add_space(SMALL_SPACING);
                            // Show examples in a collapsible section
                            ui.collapsing("Examples", |ui| {
                                ui.label(egui::RichText::new("S3:").strong().small());
                                ui.label(
                                    egui::RichText::new("  s3://bucket/experiment.d.idx")
                                        .code()
                                        .small(),
                                );
                                ui.add_space(2.0);
                                ui.label(egui::RichText::new("Google Cloud:").strong().small());
                                ui.label(
                                    egui::RichText::new("  gs://bucket/experiment.d.idx")
                                        .code()
                                        .small(),
                                );
                            });
                        });
                    }
                }

                // Show current location with clear button
                if let Some(location) = file_loader.get_raw_data_location() {
                    ui.add_space(SMALL_SPACING_VERTICAL);
                    egui::Frame::new()
                        .fill(ui.visuals().faint_bg_color)
                        .corner_radius(egui::CornerRadius::same(CORNER_RADIUS))
                        .inner_margin(egui::Margin::symmetric(
                            LOCATION_FRAME_PADDING_H,
                            LOCATION_FRAME_PADDING_V,
                        ))
                        .show(ui, |ui| {
                            ui.horizontal(|ui| {
                                let display_text =
                                    Self::truncate_middle(&location, LOCATION_DISPLAY_MAX_LEN);
                                ui.label(egui::RichText::new(display_text).small())
                                    .on_hover_text(&location);

                                ui.with_layout(
                                    egui::Layout::right_to_left(egui::Align::Center),
                                    |ui| {
                                        if ui.small_button("Clear").clicked() {
                                            file_loader.clear_raw_data();
                                        }
                                    },
                                );
                            });
                        });
                }

                Self::load_raw_data_if_needed(ui, file_loader, data, computed);

                if let IndexedDataState::Loaded { source, .. } = &data.indexed_data {
                    ui.add_space(SMALL_SPACING);
                    ui.label(
                        egui::RichText::new("✓ Raw data indexed").color(egui::Color32::DARK_GREEN),
                    );

                    // Show data source with middle truncation for long paths
                    let display_source = Self::truncate_middle(source, 60);

                    ui.label(
                        egui::RichText::new(format!("Source: {}", display_source))
                            .small()
                            .color(egui::Color32::GRAY),
                    );
                }
            });
    }

    /// Helper to truncate long paths/URLs for display
    fn truncate_middle(s: &str, max_len: usize) -> String {
        if s.len() <= max_len {
            s.to_string()
        } else {
            let start = &s[..max_len / 2 - 2];
            let end = &s[s.len() - max_len / 2 + 2..];
            format!("{}...{}", start, end)
        }
    }

    fn render_tolerance_loading_section_static(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
    ) {
        egui::Frame::group(ui.style())
            .inner_margin(egui::Margin::same(SECTION_MARGIN))
            .show(ui, |ui| {
                ui.heading("Tolerance Settings");
                ui.add_space(INTERNAL_SPACING);

                if ui.button("Load Tolerances...").clicked() {
                    file_loader.open_tolerance_dialog();
                }

                if let Some(path) = &file_loader.tolerance_path {
                    ui.add_space(SMALL_SPACING);
                    Self::display_filename(ui, path);
                }

                Self::load_tolerance_if_needed(ui, file_loader, data);

                match &data.tolerance_state {
                    ToleranceState::Loaded { .. } => {
                        ui.add_space(SMALL_SPACING);
                        ui.horizontal(|ui| {
                            ui.label(
                                egui::RichText::new("✓ Tolerance settings loaded")
                                    .color(egui::Color32::DARK_GREEN),
                            );
                            if ui.button("Reset to file").clicked()
                                && let Some(path) = &file_loader.tolerance_path
                            {
                                data.tolerance_state = ToleranceState::LoadRequested(path.clone());
                            }
                        });
                    }
                    ToleranceState::None => {
                        ui.add_space(SMALL_SPACING);
                        ui.label(
                            egui::RichText::new("Using current tolerance values")
                                .small()
                                .weak(),
                        );
                    }
                    // Failed renders error UI inside load_tolerance_if_needed above.
                    // LoadRequested always transitions out within the same frame, so it never
                    // reaches this match.
                    ToleranceState::Failed(_, _) | ToleranceState::LoadRequested(_) => {}
                }
            });
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
        ui_state: &mut UiState,
    ) {
        if let Some(path) = &file_loader.elution_groups_path {
            // Check if we need to request a new load
            let should_request_load = match &data.elution_groups {
                ElutionGroupState::None => true,
                ElutionGroupState::LoadRequested(current_path) => current_path != path,
                ElutionGroupState::Loading(current_path) => current_path != path,
                ElutionGroupState::Loaded { source, .. } => source != path,
                ElutionGroupState::Failed(_, _) => true,
            };

            // Transition to LoadRequested if new path
            if should_request_load {
                tracing::info!("Requesting load of elution groups from: {}", path.display());
                data.elution_groups = ElutionGroupState::LoadRequested(path.clone());
                // Reset selected index when loading a new library to avoid out-of-bounds errors
                ui_state.selected_index = None;
                ui.ctx().request_repaint();
            }
        }

        // Handle state machine transitions
        match &data.elution_groups {
            ElutionGroupState::LoadRequested(path) => {
                // Transition to Loading and perform the actual load
                let path = path.clone();
                tracing::info!("Starting to load elution groups from: {}", path.display());
                data.elution_groups = ElutionGroupState::Loading(path.clone());

                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Loading elution groups...");
                });

                match file_loader.load_elution_groups(&path) {
                    Ok(egs) => {
                        let count = egs.len();
                        tracing::info!("Loaded {} elution groups", count);
                        data.elution_groups = ElutionGroupState::Loaded {
                            data: Arc::new(egs),
                            source: path,
                        };
                        // Validate selected_index is within bounds of new library
                        if let Some(selected) = ui_state.selected_index
                            && selected >= count
                        {
                            tracing::warn!(
                                "Selected index {} is out of bounds for library with {} groups, resetting to None",
                                selected,
                                count
                            );
                            ui_state.selected_index = None;
                        }
                    }
                    Err(e) => {
                        let error_msg = format!("{:?}", e);
                        tracing::error!("Failed to load elution groups: {}", error_msg);
                        data.elution_groups = ElutionGroupState::Failed(path, error_msg);
                    }
                }
            }
            ElutionGroupState::Loading(_) => {
                // Load already in progress - just show spinner, don't re-trigger
                // (This state is transient - we transition out of it in the same frame)
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Loading elution groups...");
                });
            }
            ElutionGroupState::Failed(path, error) => {
                ui.label(
                    egui::RichText::new(format!("Failed to load library from {}:", path.display()))
                        .color(egui::Color32::from_rgb(255, 100, 100)),
                );
                ui.label(egui::RichText::new(error).color(egui::Color32::RED).small());
                if ui.button("Clear Error").clicked() {
                    data.elution_groups = ElutionGroupState::None;
                    file_loader.elution_groups_path = None;
                }
            }
            _ => {}
        }
    }

    fn load_raw_data_if_needed(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
        computed: &mut ComputedState,
    ) {
        if let Some(location) = file_loader.get_raw_data_location() {
            // Check if we need to request a new load
            let should_request_load = match &data.indexed_data {
                IndexedDataState::None => true,
                IndexedDataState::LoadRequested(current_location) => current_location != &location,
                IndexedDataState::Loading(current_location) => current_location != &location,
                IndexedDataState::Loaded { source, .. } => source != &location,
                IndexedDataState::Failed(_, _) => true,
            };

            // Transition to LoadRequested if new location
            if should_request_load {
                tracing::info!("Requesting load of raw data from: {}", location);
                // Clear computed state to avoid showing stale chromatograms from old index
                computed.clear();
                data.indexed_data = IndexedDataState::LoadRequested(location.clone());
                ui.ctx().request_repaint();
            }
        }

        // Handle state machine transitions
        match &data.indexed_data {
            IndexedDataState::LoadRequested(location) => {
                // Transition to Loading and perform the actual load
                let location = location.clone();
                tracing::info!("Starting to load raw data from: {}", location);
                data.indexed_data = IndexedDataState::Loading(location.clone());

                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Indexing raw data... (this may take 10-30 seconds)");
                });

                match file_loader.load_raw_data_from_location(&location) {
                    Ok(index) => {
                        data.indexed_data = IndexedDataState::Loaded {
                            index,
                            source: location,
                        };
                        file_loader.clear_raw_data();
                        tracing::info!("Raw data indexing completed");
                    }
                    Err(e) => {
                        let error_msg = format!("{:?}", e);
                        tracing::error!("Failed to load raw data: {}", error_msg);
                        data.indexed_data = IndexedDataState::Failed(location, error_msg);
                        file_loader.clear_raw_data();
                    }
                }
            }
            IndexedDataState::Loading(_) => {
                // Load already in progress - just show spinner, don't re-trigger
                // (This state is transient - we transition out of it in the same frame)
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Indexing raw data... (this may take 10-30 seconds)");
                });
            }
            IndexedDataState::Failed(location, error) => {
                ui.label(
                    egui::RichText::new(format!("Failed to load raw data from {}:", location))
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

    fn load_tolerance_if_needed(
        ui: &mut egui::Ui,
        file_loader: &mut FileLoader,
        data: &mut DataState,
    ) {
        if let Some(path) = &file_loader.tolerance_path {
            // Check if we need to request a new load
            let should_request_load = match &data.tolerance_state {
                ToleranceState::None => true,
                ToleranceState::LoadRequested(current_path) => current_path != path,
                ToleranceState::Loaded { source } => source != path,
                // Only retry if the path changed; same-path failures require the user to clear the error
                ToleranceState::Failed(failed_path, _) => failed_path != path,
            };

            if should_request_load {
                tracing::info!("Requesting load of tolerance from: {}", path.display());
                data.tolerance_state = ToleranceState::LoadRequested(path.clone());
                ui.ctx().request_repaint();
            }
        }

        // Handle state machine transitions
        match &data.tolerance_state {
            ToleranceState::LoadRequested(path) => {
                let path = path.clone();
                tracing::info!("Loading tolerance from: {}", path.display());

                match file_loader.load_tolerance(&path) {
                    Ok(tol) => {
                        tracing::info!("Tolerance loaded successfully");
                        data.tolerance = tol;
                        data.tolerance_state = ToleranceState::Loaded { source: path };
                    }
                    Err(e) => {
                        let error_msg = format!("{:?}", e);
                        tracing::error!("Failed to load tolerance: {}", error_msg);
                        data.tolerance_state = ToleranceState::Failed(path, error_msg);
                    }
                }
            }
            ToleranceState::Failed(path, error) => {
                ui.label(
                    egui::RichText::new(format!(
                        "Failed to load tolerance from {}:",
                        path.display()
                    ))
                    .color(egui::Color32::from_rgb(255, 100, 100)),
                );
                ui.label(egui::RichText::new(error).color(egui::Color32::RED).small());
                if ui.button("Clear Error").clicked() {
                    data.tolerance_state = ToleranceState::None;
                    file_loader.tolerance_path = None;
                }
            }
            _ => {}
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

    /// Handle screenshot state machine transitions and capture
    fn handle_screenshot(&mut self, ctx: &egui::Context) {
        // Check for Escape to cancel countdown
        if matches!(self.screenshot_state, ScreenshotState::Countdown { .. })
            && ctx.input(|i| i.key_pressed(egui::Key::Escape))
        {
            tracing::info!("Screenshot countdown cancelled by user");
            self.screenshot_state = ScreenshotState::Idle;
            return;
        }

        match &self.screenshot_state {
            ScreenshotState::Idle => {}
            ScreenshotState::Countdown { deadline } => {
                let now = Instant::now();
                if now >= *deadline {
                    // Deadline reached — capture at native resolution
                    ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot(
                        egui::UserData::default(),
                    ));
                    self.screenshot_state = ScreenshotState::Capturing;
                    tracing::info!("Screenshot capturing at native resolution");
                } else {
                    // Keep repainting during countdown
                    ctx.request_repaint();
                }
            }
            ScreenshotState::Capturing => {
                // Waiting for the screenshot event
            }
            ScreenshotState::Saving(image) => {
                let image = Arc::clone(image);

                // Open save dialog and write PNG
                if let Some(path) = rfd::FileDialog::new()
                    .set_title("Save Screenshot")
                    .add_filter("PNG Image", &["png"])
                    .set_file_name("screenshot.png")
                    .save_file()
                {
                    if let Err(e) = save_color_image_as_png(&image, &path) {
                        tracing::error!("Failed to save screenshot: {}", e);
                    } else {
                        tracing::info!("Screenshot saved to {}", path.display());
                    }
                }

                self.screenshot_state = ScreenshotState::Idle;
            }
        }

        // Check for screenshot events from egui
        let screenshot_image = ctx.input(|i| {
            for event in &i.raw.events {
                if let egui::Event::Screenshot { image, .. } = event {
                    return Some(Arc::clone(image));
                }
            }
            None
        });

        if let Some(image) = screenshot_image {
            tracing::info!(
                "Screenshot received: {}x{} pixels",
                image.width(),
                image.height()
            );
            self.screenshot_state = ScreenshotState::Saving(image);
            // Request repaint so the Saving state is processed next frame
            ctx.request_repaint();
        }
    }

    /// Draw countdown overlay when screenshot timer is running
    fn draw_screenshot_countdown_overlay(&self, ctx: &egui::Context) {
        let ScreenshotState::Countdown { deadline } = &self.screenshot_state else {
            return;
        };

        let remaining = deadline.saturating_duration_since(Instant::now());
        let secs = remaining.as_secs_f32().ceil() as u32;

        egui::Area::new(egui::Id::new("screenshot_countdown"))
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .order(egui::Order::Foreground)
            .show(ctx, |ui| {
                egui::Frame::new()
                    .fill(egui::Color32::from_black_alpha(180))
                    .corner_radius(egui::CornerRadius::same(12))
                    .inner_margin(egui::Margin::symmetric(40, 24))
                    .show(ui, |ui| {
                        ui.vertical_centered(|ui| {
                            ui.label(
                                egui::RichText::new(secs.to_string())
                                    .size(72.0)
                                    .color(egui::Color32::WHITE)
                                    .strong(),
                            );
                            ui.label(
                                egui::RichText::new("Press Escape to cancel")
                                    .size(14.0)
                                    .color(egui::Color32::from_white_alpha(180)),
                            );
                        });
                    });
            });
    }
}

/// Encode an egui ColorImage as PNG and write to disk
fn save_color_image_as_png(image: &egui::ColorImage, path: &std::path::Path) -> Result<(), String> {
    let width = image.width() as u32;
    let height = image.height() as u32;
    let pixels: Vec<u8> = image
        .pixels
        .iter()
        .flat_map(|c| [c.r(), c.g(), c.b(), c.a()])
        .collect();

    let img_buf: image::ImageBuffer<image::Rgba<u8>, Vec<u8>> =
        image::ImageBuffer::from_raw(width, height, pixels)
            .ok_or_else(|| "Failed to create image buffer from pixel data".to_string())?;

    img_buf
        .save(path)
        .map_err(|e| format!("Failed to write PNG: {}", e))
}

impl eframe::App for ViewerApp {
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        tracing::info!("Saving application state...");
        let state = PersistentState {
            file_loader: FileLoader {
                elution_groups_path: self.file_loader.elution_groups_path.clone(),
                raw_data_path: self.file_loader.raw_data_path.clone(),
                raw_data_url: self.file_loader.raw_data_url.clone(),
                tolerance_path: self.file_loader.tolerance_path.clone(),
            },
            ui_state: UiState {
                table_filter: self.ui.table_filter.clone(),
                selected_index: self.ui.selected_index,
                search_mode: self.ui.search_mode,
                search_input: self.ui.search_input.clone(),
                raw_data_input_mode: self.ui.raw_data_input_mode,
                dark_mode: self.ui.dark_mode,
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

        // Handle screenshot state machine
        self.handle_screenshot(ctx);

        // Check if background computation completed
        self.check_chromatogram_completion();

        // Poll calibration background thread
        if self.calibration.poll() {
            ctx.request_repaint();
        }

        // Generate MS2 spectrum and mobility data if RT was clicked
        self.handle_rt_click();

        // Generate chromatogram if needed
        self.generate_chromatogram(ctx);

        let mut tab_viewer = AppTabViewer {
            file_loader: &mut self.file_loader,
            data: &mut self.data,
            ui: &mut self.ui,
            computed: &mut self.computed,
            left_panel: &mut self.config_panel,
            table_panel: &mut self.table_panel,
            spectrum_panel: &mut self.spectrum_panel,
            mobility_panel: &mut self.mobility_panel,
            screenshot_delay_secs: &mut self.screenshot_delay_secs,
            screenshot_state: &mut self.screenshot_state,
            calibration: &mut self.calibration,
        };

        egui::CentralPanel::default().show(ctx, |ui| {
            DockArea::new(&mut self.dock_state)
                .style(Style::from_egui(ui.style().as_ref()))
                .show_inside(ui, &mut tab_viewer);
        });

        // Draw countdown overlay on top of everything
        self.draw_screenshot_countdown_overlay(ctx);
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
    mobility_panel: &'a mut MobilityPanel,
    screenshot_delay_secs: &'a mut f32,
    screenshot_state: &'a mut ScreenshotState,
    calibration: &'a mut ViewerCalibrationState,
}

impl<'a> AppTabViewer<'a> {
    fn render_left_panel(&mut self, ui: &mut egui::Ui) {
        // Data Loading Section
        ui.label(egui::RichText::new("DATA LOADING").strong().size(13.0));
        ui.add_space(INTERNAL_SPACING);

        ViewerApp::render_elution_groups_section_static(ui, self.file_loader, self.data, self.ui);
        ui.add_space(SECTION_SPACING);

        ViewerApp::render_raw_data_section_static(
            ui,
            self.file_loader,
            self.data,
            self.ui,
            self.computed,
        );
        ui.add_space(SECTION_SPACING);

        ViewerApp::render_tolerance_loading_section_static(ui, self.file_loader, self.data);

        ui.add_space(SEPARATOR_SPACING);
        ui.separator();
        ui.add_space(INTERNAL_SPACING);

        self.left_panel.render(
            ui,
            &mut self.data.tolerance,
            &mut self.data.smoothing,
            &mut self.data.auto_zoom_mode,
        );

        ui.add_space(SEPARATOR_SPACING);
        ui.separator();
        ui.add_space(INTERNAL_SPACING);

        // Export section
        let action = ConfigPanel::render_export_section(
            ui,
            self.screenshot_delay_secs,
            self.screenshot_state,
        );
        match action {
            ScreenshotAction::None => {}
            ScreenshotAction::Start => {
                let deadline = Instant::now()
                    + std::time::Duration::from_secs_f32(*self.screenshot_delay_secs);
                *self.screenshot_state = ScreenshotState::Countdown { deadline };
            }
            ScreenshotAction::Cancel => {
                *self.screenshot_state = ScreenshotState::Idle;
            }
        }

        ui.add_space(SEPARATOR_SPACING);
        ui.separator();
        ui.add_space(INTERNAL_SPACING);

        // Theme toggle
        ui.horizontal(|ui| {
            ui.label("Theme:");
            let label = if self.ui.dark_mode {
                "Switch to Light"
            } else {
                "Switch to Dark"
            };
            if ui.button(label).clicked() {
                self.ui.dark_mode = !self.ui.dark_mode;
                apply_theme(ui.ctx(), self.ui.dark_mode);
            }
        });
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
            Pane::Mobility => self.mobility_panel.title().into(),
            Pane::Calibration => "Calibration".into(),
        }
    }

    fn ui(&mut self, ui: &mut egui::Ui, tab: &mut Self::Tab) {
        tracing::trace!("rendering pane: {:?}", tab);
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
                    self.computed.ms2_spectrum(),
                    self.computed.expected_intensities(),
                );
            }
            Pane::PrecursorPlot => {
                if self.computed.is_computing() {
                    ui.centered_and_justified(|ui| {
                        ui.spinner();
                        ui.label("Computing chromatogram...");
                    });
                } else if let Some(error) = self.computed.failure_error() {
                    ui.label(
                        egui::RichText::new(format!("Chromatogram computation failed: {}", error))
                            .color(egui::Color32::RED),
                    );
                } else {
                    let click = if let Some((chromatogram, apex, counter)) =
                        self.computed.chromatogram_render_context()
                    {
                        crate::plot_renderer::render_chromatogram_plot(
                            ui,
                            chromatogram,
                            crate::plot_renderer::PlotMode::PrecursorsOnly,
                            Some("precursor_fragment_x_axis"),
                            true,
                            counter,
                            &mode,
                            &ref_lines,
                            apex,
                        )
                    } else if self.ui.selected_index.is_some() {
                        ui.label("Generating chromatogram...");
                        None
                    } else {
                        ui.label("Select a precursor from the table to view chromatogram");
                        None
                    };
                    if let Some(rt) = click {
                        self.computed.set_clicked_rt(rt);
                    }
                }
            }
            Pane::FragmentPlot => {
                if self.computed.is_computing() {
                    ui.centered_and_justified(|ui| {
                        ui.spinner();
                        ui.label("Computing chromatogram...");
                    });
                } else if let Some(error) = self.computed.failure_error() {
                    ui.label(
                        egui::RichText::new(format!("Chromatogram computation failed: {}", error))
                            .color(egui::Color32::RED),
                    );
                } else {
                    let click = if let Some((chromatogram, apex, counter)) =
                        self.computed.chromatogram_render_context()
                    {
                        crate::plot_renderer::render_chromatogram_plot(
                            ui,
                            chromatogram,
                            crate::plot_renderer::PlotMode::FragmentsOnly,
                            Some("precursor_fragment_x_axis"),
                            false,
                            counter,
                            &mode,
                            &ref_lines,
                            apex,
                        )
                    } else if self.ui.selected_index.is_some() {
                        ui.label("Generating chromatogram...");
                        None
                    } else {
                        ui.label("Select a precursor to view fragment traces");
                        None
                    };
                    if let Some(rt) = click {
                        self.computed.set_clicked_rt(rt);
                    }
                }
            }
            Pane::ScoresPlot => {
                if self.computed.is_computing() {
                    ui.centered_and_justified(|ui| {
                        ui.spinner();
                        ui.label("Computing scores...");
                    });
                } else if let Some(error) = self.computed.failure_error() {
                    ui.label(
                        egui::RichText::new(format!("Score computation failed: {}", error))
                            .color(egui::Color32::RED),
                    );
                } else {
                    let click = if let Some((score_lines, counter)) =
                        self.computed.score_render_context()
                    {
                        score_lines.render(
                            ui,
                            Some("precursor_fragment_x_axis"),
                            counter,
                            &mode,
                            &ref_lines,
                        )
                    } else if self.ui.selected_index.is_some() {
                        ui.label("Generating score plot...");
                        None
                    } else {
                        ui.label("Select a precursor to view score traces");
                        None
                    };
                    if let Some(rt) = click {
                        self.computed.set_clicked_rt(rt);
                    }
                }
            }
            Pane::Mobility => {
                self.mobility_panel
                    .render(ui, self.computed.mobility_data());
            }
            Pane::Calibration => {
                self.calibration.render_panel(
                    ui,
                    &self.data.indexed_data,
                    &self.data.elution_groups,
                    &mut self.data.tolerance,
                );
            }
        }
    }

    fn is_closeable(&self, _tab: &Self::Tab) -> bool {
        false // Tabs cannot be closed
    }
}
