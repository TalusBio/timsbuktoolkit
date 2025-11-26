use eframe::egui;
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::elution_group::ElutionGroup;
use timsquery::models::tolerance::Tolerance;

use crate::file_loader::FileLoader;
use crate::plot_renderer::ChromatogramLines;
use crate::{
    chromatogram_processor,
    plot_renderer,
    precursor_table,
    tolerance_editor,
};

/// Main application state
pub struct ViewerApp {
    /// File loader for handling file dialogs and loading
    file_loader: FileLoader,

    /// Loaded elution groups
    elution_groups: Option<Vec<ElutionGroup<usize>>>,

    /// Loaded and indexed timsTOF data
    indexed_data: Option<Arc<IndexedTimstofPeaks>>,

    /// MS1 retention times (in milliseconds)
    ms1_rts: Option<Arc<[u32]>>,

    /// Tolerance settings
    tolerance: Option<Tolerance>,

    /// Currently selected elution group index
    selected_index: Option<usize>,

    /// Computed chromatogram for the selected elution group
    chromatogram: Option<ChromatogramLines>,

    /// Track if we need to regenerate the chromatogram
    needs_regeneration: bool,

    /// Filter text for precursor table
    table_filter: String,

    /// Track if we're currently loading/indexing data
    is_indexing: bool,
}

impl ViewerApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            file_loader: FileLoader::new(),
            elution_groups: None,
            indexed_data: None,
            ms1_rts: None,
            tolerance: Some(Tolerance::default()),
            selected_index: None,
            chromatogram: None,
            needs_regeneration: false,
            table_filter: String::new(),
            is_indexing: false,
        }
    }
}

impl eframe::App for ViewerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Generate chromatogram if needed
        if self.needs_regeneration {
            self.generate_chromatogram();
            self.needs_regeneration = false;
        }

        // Top panel with title
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.heading("TimsQuery Viewer");
            ui.separator();
        });

        // Left sidebar for file loading and settings
        egui::SidePanel::left("left_panel")
            .default_width(300.0)
            .max_width(350.0)
            .resizable(true)
            .show(ctx, |ui| {
                self.render_left_panel(ui);
            });

        // Center panel for precursor table
        egui::CentralPanel::default().show(ctx, |ui| {
            self.render_table_panel(ui);
        });

        // Right panel for chromatogram plot
        egui::SidePanel::right("right_panel")
            .default_width(600.0)
            .show(ctx, |ui| {
                self.render_plot_panel(ui);
            });

    }
}

impl ViewerApp {
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

        self.render_tolerance_editor_section(ui);
    }

    fn render_elution_groups_section(&mut self, ui: &mut egui::Ui) {
        ui.label("Elution Groups:");
        if ui.button("Load Elution Groups...").clicked() {
            self.file_loader.open_elution_groups_dialog();
        }

        if let Some(path) = &self.file_loader.elution_groups_path {
            Self::display_filename(ui, path);
        }

        Self::load_elution_groups_if_needed(ui, &mut self.file_loader, &mut self.elution_groups);

        if let Some(egs) = &self.elution_groups {
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

        Self::load_raw_data_if_needed(
            ui,
            &mut self.file_loader,
            &mut self.indexed_data,
            &mut self.ms1_rts,
            &mut self.is_indexing,
        );

        if self.indexed_data.is_some() && !self.is_indexing {
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

        Self::load_tolerance_if_needed(&mut self.file_loader, &mut self.tolerance);

        if self.tolerance.is_some() {
            ui.label("Tolerance settings loaded");
        } else {
            ui.horizontal(|ui| {
                if ui.button("Use Default Tolerance").clicked() {
                    self.tolerance = Some(Tolerance::default());
                }
            });
        }
    }

    fn render_tolerance_editor_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Tolerance Settings");

        if let Some(tolerance) = &mut self.tolerance {
            let changed = tolerance_editor::render_tolerance_editor(ui, tolerance);
            if changed && self.selected_index.is_some() {
                self.needs_regeneration = true;
            }
        } else {
            ui.label("Load or create tolerance settings to edit");
        }
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
        elution_groups: &mut Option<Vec<ElutionGroup<usize>>>,
    ) {
        if let Some(path) = &file_loader.elution_groups_path
            && elution_groups.is_none()
        {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label("Loading elution groups...");
            });
            match file_loader.load_elution_groups(path) {
                Ok(egs) => {
                    tracing::info!("Loaded {} elution groups", egs.len());
                    *elution_groups = Some(egs);
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
        indexed_data: &mut Option<Arc<IndexedTimstofPeaks>>,
        ms1_rts: &mut Option<Arc<[u32]>>,
        is_indexing: &mut bool,
    ) {
        if let Some(path) = &file_loader.raw_data_path {
            if indexed_data.is_none() && !*is_indexing {
                *is_indexing = true;
                ui.ctx().request_repaint();
            }

            if *is_indexing {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Indexing raw data... (this may take 10-30 seconds)");
                });

                match file_loader.load_raw_data(path) {
                    Ok((index, rts)) => {
                        *indexed_data = Some(index);
                        *ms1_rts = Some(rts);
                        *is_indexing = false;
                        tracing::info!("Raw data indexing completed");
                    }
                    Err(e) => {
                        tracing::error!("Failed to load raw data: {:?}", e);
                        *is_indexing = false;
                    }
                }
            }
        }
    }

    fn load_tolerance_if_needed(
        file_loader: &mut FileLoader,
        tolerance: &mut Option<Tolerance>,
    ) {
        if let Some(path) = &file_loader.tolerance_path
            && tolerance.is_none()
        {
            match file_loader.load_tolerance(path) {
                Ok(tol) => {
                    *tolerance = Some(tol);
                }
                Err(e) => {
                    tracing::error!("Failed to load tolerance: {:?}", e);
                }
            }
        }
    }

    fn generate_chromatogram(&mut self) {
        let selected_idx = match self.selected_index {
            Some(idx) => idx,
            None => return,
        };

        // Clone necessary data to avoid borrowing conflicts
        let eg = match &self.elution_groups {
            Some(egs) => match egs.get(selected_idx) {
                Some(eg) => eg.clone(),
                None => return,
            },
            None => return,
        };

        let Some(index) = self.indexed_data.clone() else {
            return;
        };
        let Some(ms1_rts) = self.ms1_rts.clone() else {
            return;
        };
        let Some(tolerance) = self.tolerance.clone() else {
            return;
        };

        Self::process_chromatogram(&mut self.chromatogram, &eg, &index, ms1_rts, &tolerance);
    }

    fn process_chromatogram(
        chromatogram: &mut Option<ChromatogramLines>,
        eg: &ElutionGroup<usize>,
        index: &Arc<IndexedTimstofPeaks>,
        ms1_rts: Arc<[u32]>,
        tolerance: &Tolerance,
    ) {
        match chromatogram_processor::generate_chromatogram(eg, index, ms1_rts, tolerance) {
            Ok(chrom) => {
                tracing::info!(
                    "Generated chromatogram for elution group {} with {} precursors and {} fragments",
                    chrom.id,
                    chrom.precursor_mzs.len(),
                    chrom.fragment_mzs.len()
                );
                *chromatogram = Some(ChromatogramLines::from_chromatogram(&chrom));
            }
            Err(e) => {
                tracing::error!("Failed to generate chromatogram: {:?}", e);
                *chromatogram = None;
            }
        }
    }

    fn render_table_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Precursor Table");
        ui.separator();

        // Clone elution groups to avoid borrowing conflicts
        let Some(egs) = self.elution_groups.clone() else {
            ui.label("Load elution groups to see the table");
            return;
        };

        Self::render_table_filter_ui(ui, &mut self.table_filter);
        ui.separator();

        let filtered_egs = Self::apply_table_filter(&self.table_filter, &egs);
        let total_count = egs.len();

        ui.label(format!(
            "Showing {} of {} precursors",
            filtered_egs.len(),
            total_count
        ));

        Self::render_precursor_table_with_selection(
            ui,
            &filtered_egs,
            &mut self.selected_index,
            &mut self.needs_regeneration,
        );
    }

    fn render_table_filter_ui(ui: &mut egui::Ui, table_filter: &mut String) {
        ui.horizontal(|ui| {
            ui.label("Filter by ID:");
            ui.text_edit_singleline(table_filter);
            if ui.button("Clear").clicked() {
                table_filter.clear();
            }
        });
    }

    fn apply_table_filter<'a>(
        table_filter: &str,
        egs: &'a [ElutionGroup<usize>],
    ) -> Vec<(usize, &'a ElutionGroup<usize>)> {
        egs.iter()
            .enumerate()
            .filter(|(_, eg)| {
                table_filter.is_empty() || eg.id.to_string().contains(table_filter)
            })
            .collect()
    }

    fn render_precursor_table_with_selection(
        ui: &mut egui::Ui,
        filtered_egs: &[(usize, &ElutionGroup<usize>)],
        selected_index: &mut Option<usize>,
        needs_regeneration: &mut bool,
    ) {
        let old_selection = *selected_index;

        egui::ScrollArea::vertical()
            .auto_shrink([false; 2])
            .show(ui, |ui| {
                precursor_table::render_precursor_table_filtered(
                    ui,
                    filtered_egs,
                    selected_index,
                );
            });

        if old_selection != *selected_index && selected_index.is_some() {
            *needs_regeneration = true;
        }
    }

    fn render_plot_panel(&mut self, ui: &mut egui::Ui) {
        Self::render_plot_header(ui);
        ui.separator();
        self.render_plot_content(ui);
    }

    fn render_plot_header(ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.heading("Chromatogram");
        });
    }

    fn render_plot_content(&self, ui: &mut egui::Ui) {
        if let Some(chromatogram) = &self.chromatogram {
            plot_renderer::render_chromatogram_plot(ui, chromatogram);
        } else if self.selected_index.is_some() {
            ui.label("Generating chromatogram...");
        } else {
            ui.label("Select a precursor from the table to view chromatogram");
        }
    }
}
