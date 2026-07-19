//! Shared sensitivity-bench harness used by the `examples/` benches.
//!
//! Runs the SAME simulation config over many seeds (each a fresh noise +
//! random-peak realization), scoring every one with the real `timsseek` apex
//! finder, and measures how reliably the apex is recovered plus how fast the
//! scoring runs. Simulation is built up front so its cost stays OUT of the
//! scoring timer; the scorer is reused across runs (no per-run alloc/clone),
//! mirroring production's `ScoringWorker`.

use std::time::Instant;

use timsseek::scoring::apex_finding::TraceScorer;

use crate::scorer;
use crate::sim::{
    self,
    SimData,
    SimParams,
};

/// Aggregate results of a sensitivity run.
pub struct SensitivityReport {
    pub n: usize,
    pub tol: i64,
    pub true_apex: i64,
    pub build_ms: f64,
    pub score_ms: f64,
    pub errors: usize,
    pub pass1_hits: usize,
    pub pass2_hits: usize,
    /// Sorted pass-2 |apex - truth| per successful run.
    pub abs_err: Vec<i64>,
    // Echoed config for the printout.
    pub n_cycles: usize,
    pub width_sigma: f32,
    pub noise_floor: f32,
    pub random_peaks: usize,
}

fn pct(num: usize, den: usize) -> f64 {
    if den == 0 {
        0.0
    } else {
        100.0 * num as f64 / den as f64
    }
}

/// The canonical benchmark scenarios, run together by `examples/bench.rs`.
/// All share the same peptide geometry (245 cycles, apex@30, sigma 1, seed
/// 680); each isolates one stressor so the summary table reads as an ablation.
pub fn canonical_suite() -> Vec<(&'static str, SimParams)> {
    let base = || {
        let mut p = SimParams::default();
        p.n_cycles = 245;
        p.apex_cycle = 30.0;
        p.width_sigma = 1.0;
        p.seed = 680;
        p
    };

    // Baseline: low noise, no interference -> should be near-perfect.
    let mut clean = base();
    clean.noise_floor = 0.02;
    clean.random_peaks.enabled = false;

    // Moderate noise + light interference.
    let mut moderate = base();
    moderate.noise_floor = 0.2;
    moderate.random_peaks.count = 50;

    // High noise + heavy interference (the main stress condition).
    let mut stress = base();
    stress.noise_floor = 0.5;
    stress.random_peaks.count = 100;

    // Even more injected peaks at the same noise.
    let mut heavy = base();
    heavy.noise_floor = 0.5;
    heavy.random_peaks.count = 300;

    // Mismatched library at stress level: flat theoretical ratios, observed
    // keeps the descending ladder (cross-instrument speclib case).
    let mut mismatched = base();
    mismatched.noise_floor = 0.5;
    mismatched.random_peaks.count = 100;
    for f in &mut mismatched.real_fragments {
        f.obs_scale = f.theo_intensity;
        f.theo_intensity = 1.0;
    }

    // Absent top fragment at stress level (expected but observed as noise).
    let mut absent = base();
    absent.noise_floor = 0.5;
    absent.random_peaks.count = 100;
    absent.real_fragments[0].obs_scale = 0.0;

    // Absent precursor at stress level: MS1 precursor undetected (intensity 0),
    // fragments intact. Guards against over-reliance on precursor signal.
    let mut absent_prec = base();
    absent_prec.noise_floor = 0.5;
    absent_prec.random_peaks.count = 100;
    absent_prec.precursor_intensity = 0.0;

    vec![
        ("clean", clean),
        ("moderate_noise", moderate),
        ("high_noise+interference", stress),
        ("heavy_interference", heavy),
        ("mismatched_library", mismatched),
        ("absent_top_fragment", absent),
        ("absent_precursor", absent_prec),
    ]
}

/// Production-realistic BROAD apex-finding suite (Phase-1 window ~1695 cycles,
/// RT-unrestricted). Interferents are DENSITY-scaled (fixed per-cycle rate) and
/// the apex is jittered per seed with a sub-cycle offset, so widening the window
/// does not make the task artificially easy and the finder faces a moving,
/// off-grid target. Each scenario is its own line so it can be commented out.
pub fn broad_suite() -> Vec<(&'static str, SimParams)> {
    let base = || {
        let mut p = SimParams::default();
        p.n_cycles = 1695;
        p.apex_cycle = 847.0;
        p.apex_jitter = Some(800.0);
        p.width_sigma = 1.0;
        p.seed = 680;
        p.random_peaks.count = 0;
        p.random_peaks.density_per_cycle = 2.5; // ~ current 100/40 density
        p
    };

    let mut clean = base();
    clean.noise_floor = 0.02;
    clean.random_peaks.enabled = false;

    let mut moderate = base();
    moderate.noise_floor = 0.2;

    let mut stress = base();
    stress.noise_floor = 0.5;

    let mut hard = base();
    hard.noise_floor = 0.5;
    hard.random_peaks.hardness = 3.0; // 3x interferent density

    let mut mismatched = base();
    mismatched.noise_floor = 0.5;
    for f in &mut mismatched.real_fragments {
        f.obs_scale = f.theo_intensity;
        f.theo_intensity = 1.0;
    }

    vec![
        ("broad_clean", clean),
        ("broad_moderate_noise", moderate),
        ("broad_high_noise+interf", stress),
        ("broad_hard_3x_density", hard),
        ("broad_mismatched_library", mismatched),
    ]
}

/// Production-realistic NARROW scoring suite (Phase-3 calibrated window ~150
/// cycles). Same density-scaled interferents + jittered sub-cycle apex; used
/// for BOTH apex recovery and score discrimination (`run_discrimination`).
pub fn narrow_suite() -> Vec<(&'static str, SimParams)> {
    let base = || {
        let mut p = SimParams::default();
        p.n_cycles = 150;
        p.apex_cycle = 75.0;
        p.apex_jitter = Some(60.0);
        p.width_sigma = 1.0;
        p.seed = 680;
        p.random_peaks.count = 0;
        p.random_peaks.density_per_cycle = 2.5;
        p
    };

    let mut clean = base();
    clean.noise_floor = 0.02;
    clean.random_peaks.enabled = false;

    let mut moderate = base();
    moderate.noise_floor = 0.2;

    let mut stress = base();
    stress.noise_floor = 0.5;

    let mut hard = base();
    hard.noise_floor = 0.5;
    hard.random_peaks.hardness = 3.0;

    let mut mismatched = base();
    mismatched.noise_floor = 0.5;
    for f in &mut mismatched.real_fragments {
        f.obs_scale = f.theo_intensity;
        f.theo_intensity = 1.0;
    }

    vec![
        ("narrow_clean", clean),
        ("narrow_moderate_noise", moderate),
        ("narrow_high_noise+interf", stress),
        ("narrow_hard_3x_density", hard),
        ("narrow_mismatched_library", mismatched),
    ]
}

/// Run `n_runs` scored simulations of `base` (varying only the seed) and
/// collect apex-recovery + timing stats. A hit = pass-2 apex within `tol`
/// cycles of the configured `apex_cycle`.
pub fn run_sensitivity(base: &SimParams, n_runs: usize, tol: i64) -> SensitivityReport {
    // Echoed in the report header; each run is scored against its OWN realized
    // apex (jitter varies it per seed), not this nominal value.
    let nominal_apex = base.apex_cycle.round() as i64;
    // rt mapping is seed-independent, so one mapper serves every run.
    let map = base.rt_mapper();

    // Phase 1: build ALL extractions up front (kept out of the timed loop).
    let t_build = Instant::now();
    let sims: Vec<SimData> = (0..n_runs)
        .map(|i| {
            let mut p = base.clone();
            p.seed = base.seed.wrapping_add(i as u64); // independent realization
            sim::build(&p)
        })
        .collect();
    let build_ms = t_build.elapsed().as_secs_f64() * 1e3;

    // Phase 2: time ONLY the scoring. One reusable scorer, no per-run alloc.
    let mut errors = 0;
    let mut pass1_hits = 0;
    let mut pass2_hits = 0;
    let mut abs_err = Vec::with_capacity(n_runs);
    let mut scorer = TraceScorer::new(base.n_cycles, base.real_fragments.len().max(1));
    let t_score = Instant::now();
    for data in &sims {
        let true_apex = data.realized_apex_cycle.round() as i64;
        match scorer::run_with(&mut scorer, &data.extraction, &map) {
            Ok((pass1, pass2)) => {
                if (pass1.apex_cycle as i64 - true_apex).abs() <= tol {
                    pass1_hits += 1;
                }
                let e2 = (pass2.joint_apex_cycle as i64 - true_apex).abs();
                if e2 <= tol {
                    pass2_hits += 1;
                }
                abs_err.push(e2);
            }
            Err(_) => errors += 1,
        }
    }
    let score_ms = t_score.elapsed().as_secs_f64() * 1e3;

    abs_err.sort_unstable();
    SensitivityReport {
        n: n_runs,
        tol,
        true_apex: nominal_apex,
        build_ms,
        score_ms,
        errors,
        pass1_hits,
        pass2_hits,
        abs_err,
        n_cycles: base.n_cycles,
        width_sigma: base.width_sigma,
        noise_floor: base.noise_floor,
        random_peaks: base.random_peaks.count,
    }
}

impl SensitivityReport {
    pub fn pass1_pct(&self) -> f64 {
        pct(self.pass1_hits, self.n)
    }

    pub fn pass2_pct(&self) -> f64 {
        pct(self.pass2_hits, self.n)
    }

    pub fn median_err(&self) -> i64 {
        self.abs_err
            .get(self.abs_err.len() / 2)
            .copied()
            .unwrap_or(-1)
    }

    pub fn score_us_per_run(&self) -> f64 {
        self.score_ms * 1e3 / self.n.max(1) as f64
    }

    /// Print a human-readable summary under `title`.
    pub fn print(&self, title: &str) {
        let median = self
            .abs_err
            .get(self.abs_err.len() / 2)
            .copied()
            .unwrap_or(-1);
        let mean = if self.abs_err.is_empty() {
            f64::NAN
        } else {
            self.abs_err.iter().sum::<i64>() as f64 / self.abs_err.len() as f64
        };
        let max = self.abs_err.last().copied().unwrap_or(-1);
        let score_us = self.score_ms * 1e3 / self.n.max(1) as f64;
        let build_us = self.build_ms * 1e3 / self.n.max(1) as f64;
        let throughput = self.n as f64 / (self.score_ms / 1e3).max(1e-12);

        println!("=== {title} ===");
        println!(
            "  config: n_cycles={} sigma={} noise_floor={} random_peaks={} true_apex={}",
            self.n_cycles, self.width_sigma, self.noise_floor, self.random_peaks, self.true_apex,
        );
        println!("  runs={} tol=±{} cycles", self.n, self.tol);
        println!(
            "  build : {:>8.2} ms total  ({:.2} us/run)",
            self.build_ms, build_us
        );
        println!(
            "  score : {:>8.2} ms total  ({:.2} us/run, {:.0} runs/s)",
            self.score_ms, score_us, throughput
        );
        println!("  scorer errors:    {}", self.errors);
        println!(
            "  pass1 sensitivity: {:>6.2}%  ({}/{})",
            pct(self.pass1_hits, self.n),
            self.pass1_hits,
            self.n
        );
        println!(
            "  pass2 sensitivity: {:>6.2}%  ({}/{})",
            pct(self.pass2_hits, self.n),
            self.pass2_hits,
            self.n
        );
        println!(
            "  pass2 |apex-truth| cycles: mean={:.2} median={} max={}",
            mean, median, max
        );

        let mut buckets = [0usize; 8];
        for &e in &self.abs_err {
            buckets[(e.max(0) as usize).min(7)] += 1;
        }
        print!("  pass2 err histogram: ");
        for (i, c) in buckets.iter().enumerate() {
            let label = if i == 7 { "7+".into() } else { i.to_string() };
            print!("{label}:{c} ");
        }
        println!();
    }
}

/// Present-vs-absent score-discrimination result for one scenario.
pub struct DiscriminationReport {
    pub n_pairs: usize,
    /// ROC-AUC = P(score_present > score_absent), 0.5 tie credit. 1.0 = perfect
    /// signal/noise separation, 0.5 = useless.
    pub auc: f64,
    pub median_present: f32,
    pub median_absent: f32,
}

/// ROC-AUC of `present` vs `absent` score populations via the average-rank
/// Mann-Whitney statistic: `P(present > absent)` with 0.5 credit for ties.
pub fn roc_auc(present: &[f32], absent: &[f32]) -> f64 {
    let (np, na) = (present.len(), absent.len());
    if np == 0 || na == 0 {
        return f64::NAN;
    }
    // Combined values tagged present(true)/absent(false), sorted ascending.
    let mut all: Vec<(f32, bool)> = present
        .iter()
        .map(|&v| (v, true))
        .chain(absent.iter().map(|&v| (v, false)))
        .collect();
    all.sort_by(|a, b| a.0.total_cmp(&b.0));
    // Sum of average ranks (1-based) assigned to present values.
    let mut rank_sum_present = 0.0f64;
    let mut i = 0usize;
    while i < all.len() {
        let mut j = i + 1;
        while j < all.len() && all[j].0 == all[i].0 {
            j += 1;
        }
        // Mean rank of the tie block covering positions (i+1..=j).
        let avg_rank = ((i + 1 + j) as f64) / 2.0;
        for entry in &all[i..j] {
            if entry.1 {
                rank_sum_present += avg_rank;
            }
        }
        i = j;
    }
    let u = rank_sum_present - (np * (np + 1)) as f64 / 2.0;
    u / (np as f64 * na as f64)
}

/// Score `n_pairs` matched present/absent realizations and report how well the
/// production Pass-2 score separates real signal from pure noise. Each pair
/// shares a seed, so the "absent" twin (all real fragments `obs_scale = 0`) has
/// IDENTICAL noise + interferents — only the real peak differs.
pub fn run_discrimination(base: &SimParams, n_pairs: usize) -> DiscriminationReport {
    let map = base.rt_mapper();
    let mut scorer = TraceScorer::new(base.n_cycles, base.real_fragments.len().max(1));
    let mut present: Vec<f32> = Vec::with_capacity(n_pairs);
    let mut absent: Vec<f32> = Vec::with_capacity(n_pairs);
    for i in 0..n_pairs {
        let mut pp = base.clone();
        pp.seed = base.seed.wrapping_add(i as u64);
        // Pure-noise twin: NO real signal at all (fragments AND precursor
        // zeroed), identical seed => identical noise + interferents. Zeroing
        // only fragments would leave the precursor isotope peaks in place and
        // the score would still see real signal.
        let mut ap = pp.clone();
        for f in &mut ap.real_fragments {
            f.obs_scale = 0.0;
        }
        ap.precursor_intensity = 0.0;
        if let Ok((_, s)) = scorer::run_with(&mut scorer, &sim::build(&pp).extraction, &map) {
            present.push(s.score);
        }
        if let Ok((_, s)) = scorer::run_with(&mut scorer, &sim::build(&ap).extraction, &map) {
            absent.push(s.score);
        }
    }
    let median = |mut v: Vec<f32>| {
        v.sort_by(f32::total_cmp);
        v.get(v.len() / 2).copied().unwrap_or(f32::NAN)
    };
    let auc = roc_auc(&present, &absent);
    DiscriminationReport {
        n_pairs,
        auc,
        median_present: median(present),
        median_absent: median(absent),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roc_auc_basic() {
        assert!((roc_auc(&[3.0, 4.0, 5.0], &[0.0, 1.0, 2.0]) - 1.0).abs() < 1e-9);
        assert!((roc_auc(&[0.0, 1.0], &[2.0, 3.0]) - 0.0).abs() < 1e-9);
        assert!((roc_auc(&[1.0, 1.0], &[1.0, 1.0]) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn discrimination_high_when_signal_present() {
        // Genuinely clean: low noise, no injected interferents. A real peptide
        // vs pure noise must separate near-perfectly.
        let mut base = SimParams::default();
        base.n_cycles = 150;
        base.noise_floor = 0.05;
        base.random_peaks.enabled = false;
        let rep = run_discrimination(&base, 200);
        assert!(
            rep.auc > 0.9,
            "clean signal should separate: auc={}",
            rep.auc
        );
        assert!(rep.median_present > rep.median_absent);
    }
}
