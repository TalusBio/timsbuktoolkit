use clap::Parser;
use std::path::PathBuf;

/// Build a spectral library from FASTA or peptide list using Koina predictions.
#[derive(Debug, Parser)]
#[command(name = "speclib_build", version, about)]
pub struct Cli {
    // ── Input ──────────────────────────────────────────────────────────────
    /// FASTA file to digest into peptides.
    #[arg(long)]
    pub fasta: Option<PathBuf>,

    /// Pre-digested peptide list (one bare sequence per line).
    #[arg(long)]
    pub peptide_list: Option<PathBuf>,

    // ── Config ─────────────────────────────────────────────────────────────
    /// TOML config file; CLI flags override values from this file.
    #[arg(long)]
    pub config: Option<PathBuf>,

    // ── Output ─────────────────────────────────────────────────────────────
    /// Output path for the spectral library (default: library.msgpack.zst).
    #[arg(long, short = 'o')]
    pub output: Option<PathBuf>,

    // ── Modifications ──────────────────────────────────────────────────────
    /// Fixed modification, e.g. `Carbamidomethyl@C`. Repeatable.
    #[arg(long = "fixed-mod")]
    pub fixed_mods: Vec<String>,

    /// Variable modification, e.g. `Oxidation@M`. Repeatable.
    #[arg(long = "var-mod")]
    pub var_mods: Vec<String>,

    /// Maximum number of variable modifications per peptide.
    #[arg(long)]
    pub max_var_mods: Option<usize>,

    // ── Digestion ──────────────────────────────────────────────────────────
    /// Minimum peptide length (amino acids).
    #[arg(long)]
    pub min_length: Option<usize>,

    /// Maximum peptide length (amino acids).
    #[arg(long)]
    pub max_length: Option<usize>,

    /// Maximum missed cleavages allowed.
    #[arg(long)]
    pub missed_cleavages: Option<usize>,

    // ── Charges ────────────────────────────────────────────────────────────
    /// Minimum precursor charge state.
    #[arg(long)]
    pub min_charge: Option<u8>,

    /// Maximum precursor charge state.
    #[arg(long)]
    pub max_charge: Option<u8>,

    // ── Decoys ─────────────────────────────────────────────────────────────
    /// Decoy generation strategy: none | reverse | edge_mutate.
    #[arg(long)]
    pub decoy_strategy: Option<String>,

    // ── Prediction ─────────────────────────────────────────────────────────
    /// Koina fragment intensity model name.
    #[arg(long)]
    pub fragment_model: Option<String>,

    /// Koina retention-time model name.
    #[arg(long)]
    pub rt_model: Option<String>,

    /// Base URL for the Koina prediction service.
    #[arg(long)]
    pub koina_url: Option<String>,

    /// Number of peptides per Koina request batch.
    #[arg(long)]
    pub batch_size: Option<usize>,

    /// Normalised collision energy (0.0–1.0).
    #[arg(long)]
    pub nce: Option<f32>,

    // ── Fragment filters ───────────────────────────────────────────────────
    /// Maximum number of fragment ions retained per precursor.
    #[arg(long)]
    pub max_ions: Option<usize>,

    /// Minimum precursor m/z to include.
    #[arg(long)]
    pub min_mz: Option<f32>,

    /// Maximum precursor m/z to include.
    #[arg(long)]
    pub max_mz: Option<f32>,

    /// Minimum fragment ion m/z to retain.
    #[arg(long)]
    pub min_ion_mz: Option<f32>,

    /// Maximum fragment ion m/z to retain.
    #[arg(long)]
    pub max_ion_mz: Option<f32>,

    /// Minimum number of fragment ions required to keep a precursor.
    #[arg(long)]
    pub min_ions: Option<usize>,
}
