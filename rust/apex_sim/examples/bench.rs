//! Apex-finder bench: runs the whole catalog of canonical scenarios
//! (`apex_sim::bench::canonical_suite`) at once and prints per-scenario detail
//! plus a comparison summary table.
//!
//! Each scenario runs `n_runs` seeds of one fixed config; "sensitivity" is the
//! fraction of runs whose pass-2 apex lands within `tol` cycles of truth.
//!
//! Run:
//!   cargo run -p apex_sim --release --example bench
//!   cargo run -p apex_sim --release --example bench -- 500 2
//!     (positional: <n_runs> <tolerance_cycles>)

use apex_sim::bench::{
    canonical_suite,
    run_sensitivity,
};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n_runs: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(2000);
    let tol: i64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2);

    let suite = canonical_suite();
    let mut rows = Vec::with_capacity(suite.len());
    for (name, cfg) in &suite {
        let report = run_sensitivity(cfg, n_runs, tol);
        report.print(name);
        println!();
        rows.push((*name, report));
    }

    println!("=== summary (n={n_runs}, tol=±{tol} cycles) ===");
    println!(
        "  {:<26} {:>7} {:>7} {:>8} {:>9}",
        "scenario", "pass2%", "pass1%", "medErr", "us/run"
    );
    for (name, r) in &rows {
        println!(
            "  {:<26} {:>7.1} {:>7.1} {:>8} {:>9.2}",
            name,
            r.pass2_pct(),
            r.pass1_pct(),
            r.median_err(),
            r.score_us_per_run(),
        );
    }
}
