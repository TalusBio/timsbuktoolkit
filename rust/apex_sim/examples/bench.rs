//! Apex-finder + score benches.
//!
//! Runs three suites:
//!   * canonical (historical 245-cyc scenarios; apex-recovery only),
//!   * broad (production Phase-1 window ~1695 cyc; apex-recovery), and
//!   * narrow (production Phase-3 window ~150 cyc; apex-recovery + score AUC).
//! "recovery" = fraction of runs whose pass-2 apex lands within `tol` cycles of
//! the (jittered) truth. "AUC" = ROC-AUC separating the Pass-2 score of a real
//! peptide from a pure-noise twin (signal-vs-noise discrimination).
//!
//! Run:
//!   cargo run -p apex_sim --release --example bench
//!   cargo run -p apex_sim --release --example bench -- 2000 2 500
//!     (positional: <n_runs> <tolerance_cycles> <broad_n_runs>)

use apex_sim::bench::{
    SensitivityReport,
    broad_suite,
    canonical_suite,
    narrow_suite,
    run_discrimination,
    run_sensitivity,
};
use apex_sim::sim::SimParams;

/// Run a suite and print a comparison summary table under `=== {title} ===`.
/// The canonical suite passes `title = "summary"` so `bench_out/score.sh`
/// (which keys on `=== summary`) keeps working.
fn run_and_print_recovery(
    suite: &[(&'static str, SimParams)],
    n_runs: usize,
    tol: i64,
    title: &str,
) {
    let rows: Vec<(&str, SensitivityReport)> = suite
        .iter()
        .map(|(name, cfg)| (*name, run_sensitivity(cfg, n_runs, tol)))
        .collect();

    println!("=== {title} (n={n_runs}, tol=±{tol} cycles) ===");
    println!(
        "  {:<28} {:>7} {:>7} {:>8} {:>9}",
        "scenario", "pass2%", "pass1%", "medErr", "us/run"
    );
    for (name, r) in &rows {
        println!(
            "  {:<28} {:>7.1} {:>7.1} {:>8} {:>9.2}",
            name,
            r.pass2_pct(),
            r.pass1_pct(),
            r.median_err(),
            r.score_us_per_run(),
        );
    }
    println!();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n_runs: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(2000);
    let tol: i64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2);
    // Broad sims are ~42x heavier per run; default to fewer.
    let broad_runs: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(500);

    // Canonical (historical) apex-recovery suite. Header stays `=== summary`.
    run_and_print_recovery(&canonical_suite(), n_runs, tol, "summary");

    // Broad apex-finding (production Phase-1 window ~1695 cyc).
    run_and_print_recovery(&broad_suite(), broad_runs, tol, "broad apex-finding");

    // Narrow scoring (production Phase-3 window ~150 cyc): apex recovery.
    run_and_print_recovery(&narrow_suite(), n_runs, tol, "narrow recovery");

    // Narrow score discrimination: signal-present vs pure-noise (ROC-AUC).
    println!("=== narrow score discrimination (AUC, n_pairs={n_runs}) ===");
    println!(
        "  {:<28} {:>7} {:>12} {:>12}",
        "scenario", "AUC", "med+signal", "med-noise"
    );
    for (name, cfg) in &narrow_suite() {
        let d = run_discrimination(cfg, n_runs);
        println!(
            "  {:<28} {:>7.3} {:>12.3e} {:>12.3e}",
            name, d.auc, d.median_present, d.median_absent
        );
    }
}
