use serde::Deserialize;
use std::path::PathBuf;

use crate::cli::Cli;

// ── Default value helpers ────────────────────────────────────────────────────

fn default_enzyme() -> String {
    "trypsin".to_string()
}
fn default_min_length() -> usize {
    7
}
fn default_max_length() -> usize {
    25
}
fn default_missed_cleavages() -> usize {
    1
}
fn default_max_variable_mods() -> usize {
    2
}
fn default_min_charge() -> u8 {
    2
}
fn default_max_charge() -> u8 {
    4
}
fn default_decoy_strategy() -> String {
    "none".to_string()
}
fn default_fragment_model() -> String {
    "Prosit_2020_intensity_HCD".to_string()
}
fn default_rt_model() -> String {
    "Prosit_2019_irt".to_string()
}
fn default_koina_url() -> String {
    "https://koina.wilhelmlab.org/v2/models".to_string()
}
fn default_batch_size() -> usize {
    1000
}
fn default_nce() -> f32 {
    0.3
}
fn default_max_ions() -> usize {
    10
}
fn default_min_mz() -> f32 {
    400.0
}
fn default_max_mz() -> f32 {
    2000.0
}
fn default_min_ion_mz() -> f32 {
    250.0
}
fn default_max_ion_mz() -> f32 {
    2000.0
}
fn default_min_ions() -> usize {
    3
}
fn default_output() -> PathBuf {
    PathBuf::from("library.msgpack.zst")
}

// ── Sub-structs ──────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct DigestionConfig {
    /// Enzyme name (currently only "trypsin" is recognised).
    pub enzyme: String,
    pub min_length: usize,
    pub max_length: usize,
    pub missed_cleavages: usize,
}

impl Default for DigestionConfig {
    fn default() -> Self {
        Self {
            enzyme: default_enzyme(),
            min_length: default_min_length(),
            max_length: default_max_length(),
            missed_cleavages: default_missed_cleavages(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct ModificationsConfig {
    /// Fixed modifications applied to every matching residue, e.g. ["Carbamidomethyl@C"].
    pub fixed: Vec<String>,
    /// Variable modifications considered during peptide generation, e.g. ["Oxidation@M"].
    pub variable: Vec<String>,
    /// Maximum number of variable modifications per peptide.
    pub max_variable: usize,
}

impl Default for ModificationsConfig {
    fn default() -> Self {
        Self {
            fixed: Vec::new(),
            variable: Vec::new(),
            max_variable: default_max_variable_mods(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct ChargesConfig {
    pub min: u8,
    pub max: u8,
}

impl Default for ChargesConfig {
    fn default() -> Self {
        Self {
            min: default_min_charge(),
            max: default_max_charge(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct DecoysConfig {
    /// Strategy for generating decoy entries: "none", "reverse", or "edge_mutate".
    pub strategy: String,
}

impl Default for DecoysConfig {
    fn default() -> Self {
        Self {
            strategy: default_decoy_strategy(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct PredictionConfig {
    pub fragment_model: String,
    pub rt_model: String,
    pub koina_url: String,
    pub batch_size: usize,
    /// Normalised collision energy sent to the Koina model (0.0–1.0).
    pub nce: f32,
}

impl Default for PredictionConfig {
    fn default() -> Self {
        Self {
            fragment_model: default_fragment_model(),
            rt_model: default_rt_model(),
            koina_url: default_koina_url(),
            batch_size: default_batch_size(),
            nce: default_nce(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct FiltersConfig {
    pub max_ions: usize,
    pub min_mz: f32,
    pub max_mz: f32,
    pub min_ion_mz: f32,
    pub max_ion_mz: f32,
    /// Discard precursors with fewer than this many fragment ions surviving filters.
    pub min_ions: usize,
}

impl Default for FiltersConfig {
    fn default() -> Self {
        Self {
            max_ions: default_max_ions(),
            min_mz: default_min_mz(),
            max_mz: default_max_mz(),
            min_ion_mz: default_min_ion_mz(),
            max_ion_mz: default_max_ion_mz(),
            min_ions: default_min_ions(),
        }
    }
}

// ── Top-level config ─────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct SpeclibBuildConfig {
    // Inputs — not directly deserialised from TOML but set after merging CLI args.
    #[serde(skip)]
    pub fasta: Option<PathBuf>,
    #[serde(skip)]
    pub peptide_list: Option<PathBuf>,

    pub output: PathBuf,
    pub digestion: DigestionConfig,
    pub modifications: ModificationsConfig,
    pub charges: ChargesConfig,
    pub decoys: DecoysConfig,
    pub prediction: PredictionConfig,
    pub filters: FiltersConfig,
}

impl Default for SpeclibBuildConfig {
    fn default() -> Self {
        Self {
            fasta: None,
            peptide_list: None,
            output: default_output(),
            digestion: DigestionConfig::default(),
            modifications: ModificationsConfig::default(),
            charges: ChargesConfig::default(),
            decoys: DecoysConfig::default(),
            prediction: PredictionConfig::default(),
            filters: FiltersConfig::default(),
        }
    }
}

impl SpeclibBuildConfig {
    /// Build a config by loading an optional TOML file and then overlaying CLI args.
    ///
    /// Precedence: CLI flags > TOML file > compiled-in defaults.
    pub fn from_cli(cli: &Cli) -> Result<Self, Box<dyn std::error::Error>> {
        // Start from TOML file if provided, otherwise use defaults.
        let mut cfg: SpeclibBuildConfig = if let Some(path) = &cli.config {
            let raw = std::fs::read_to_string(path)
                .map_err(|e| format!("cannot read config file {}: {e}", path.display()))?;
            toml::from_str(&raw)
                .map_err(|e| format!("invalid TOML in {}: {e}", path.display()))?
        } else {
            SpeclibBuildConfig::default()
        };

        // ── Inputs ──
        cfg.fasta = cli.fasta.clone();
        cfg.peptide_list = cli.peptide_list.clone();

        // ── Output ──
        if let Some(v) = &cli.output {
            cfg.output = v.clone();
        }

        // ── Digestion ──
        if let Some(v) = cli.min_length {
            cfg.digestion.min_length = v;
        }
        if let Some(v) = cli.max_length {
            cfg.digestion.max_length = v;
        }
        if let Some(v) = cli.missed_cleavages {
            cfg.digestion.missed_cleavages = v;
        }

        // ── Modifications ──
        if !cli.fixed_mods.is_empty() {
            cfg.modifications.fixed = cli.fixed_mods.clone();
        }
        if !cli.var_mods.is_empty() {
            cfg.modifications.variable = cli.var_mods.clone();
        }
        if let Some(v) = cli.max_var_mods {
            cfg.modifications.max_variable = v;
        }

        // ── Charges ──
        if let Some(v) = cli.min_charge {
            cfg.charges.min = v;
        }
        if let Some(v) = cli.max_charge {
            cfg.charges.max = v;
        }

        // ── Decoys ──
        if let Some(ref v) = cli.decoy_strategy {
            cfg.decoys.strategy = v.clone();
        }

        // ── Prediction ──
        if let Some(ref v) = cli.fragment_model {
            cfg.prediction.fragment_model = v.clone();
        }
        if let Some(ref v) = cli.rt_model {
            cfg.prediction.rt_model = v.clone();
        }
        if let Some(ref v) = cli.koina_url {
            cfg.prediction.koina_url = v.clone();
        }
        if let Some(v) = cli.batch_size {
            cfg.prediction.batch_size = v;
        }
        if let Some(v) = cli.nce {
            cfg.prediction.nce = v;
        }

        // ── Filters ──
        if let Some(v) = cli.max_ions {
            cfg.filters.max_ions = v;
        }
        if let Some(v) = cli.min_mz {
            cfg.filters.min_mz = v;
        }
        if let Some(v) = cli.max_mz {
            cfg.filters.max_mz = v;
        }
        if let Some(v) = cli.min_ion_mz {
            cfg.filters.min_ion_mz = v;
        }
        if let Some(v) = cli.max_ion_mz {
            cfg.filters.max_ion_mz = v;
        }
        if let Some(v) = cli.min_ions {
            cfg.filters.min_ions = v;
        }

        Ok(cfg)
    }

    /// Validate logical constraints across the merged config.
    pub fn validate(&self) -> Result<(), String> {
        // Exactly one input source must be provided.
        match (&self.fasta, &self.peptide_list) {
            (None, None) => {
                return Err(
                    "provide exactly one input: --fasta or --peptide-list".to_string()
                )
            }
            (Some(_), Some(_)) => {
                return Err(
                    "--fasta and --peptide-list are mutually exclusive".to_string()
                )
            }
            _ => {}
        }

        if self.digestion.min_length == 0 {
            return Err("min_length must be >= 1".to_string());
        }
        if self.digestion.max_length < self.digestion.min_length {
            return Err(format!(
                "max_length ({}) must be >= min_length ({})",
                self.digestion.max_length, self.digestion.min_length
            ));
        }
        if self.charges.max < self.charges.min {
            return Err(format!(
                "max_charge ({}) must be >= min_charge ({})",
                self.charges.max, self.charges.min
            ));
        }
        let valid_strategies = ["none", "reverse", "edge_mutate"];
        if !valid_strategies.contains(&self.decoys.strategy.as_str()) {
            return Err(format!(
                "unknown decoy strategy {:?}; valid values: {}",
                self.decoys.strategy,
                valid_strategies.join(", ")
            ));
        }
        if !(0.0..=1.0).contains(&self.prediction.nce) {
            return Err(format!(
                "nce ({}) must be in [0.0, 1.0]",
                self.prediction.nce
            ));
        }
        if self.filters.min_ions == 0 {
            return Err("min_ions must be >= 1".to_string());
        }
        if self.filters.max_ions < self.filters.min_ions {
            return Err(format!(
                "max_ions ({}) must be >= min_ions ({})",
                self.filters.max_ions, self.filters.min_ions
            ));
        }
        if self.filters.min_mz >= self.filters.max_mz {
            return Err(format!(
                "min_mz ({}) must be < max_mz ({})",
                self.filters.min_mz, self.filters.max_mz
            ));
        }
        if self.filters.min_ion_mz >= self.filters.max_ion_mz {
            return Err(format!(
                "min_ion_mz ({}) must be < max_ion_mz ({})",
                self.filters.min_ion_mz, self.filters.max_ion_mz
            ));
        }

        Ok(())
    }
}
