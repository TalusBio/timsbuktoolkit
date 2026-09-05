//! Search-free stress test for the RT path finder.
//!
//! Emits CSV to stdout. A known monotone ridge is mixed with equal-weight
//! uniform background; no target labels or search scores are involved.

use calibrt::{
    CalibrationState,
    LibraryRT,
    ObservedRTSeconds,
    RidgeSummary,
};
use rand::rngs::StdRng;
use rand::{
    Rng,
    SeedableRng,
};

const SPAN: f64 = 1_000.0;
const N_SAMPLES: usize = 200;
const N_REPLICATES: u64 = 10;
const SIGNAL_COUNTS: &[usize] = &[250, 1_000, 2_000];
const NOISE_RATIOS: &[usize] = &[0, 1, 2, 4, 8, 16, 32, 64, 128];
const GEOMETRIES: &[(usize, usize)] = &[(50, 15), (100, 30), (200, 60)];

fn truth(x: f64) -> f64 {
    0.8 * x + 100.0
}

fn main() {
    println!(
        "bins,lookback,signal,noise,replicate,seed,path_nodes,coverage,mean_abs_error,max_abs_error,ridge_half_width,in_ridge_ratio"
    );

    for &(bins, lookback) in GEOMETRIES {
        let signal_jitter = SPAN / bins as f64 / 4.0;
        for &signal_count in SIGNAL_COUNTS {
            for &noise_ratio in NOISE_RATIOS {
                let noise_count = signal_count * noise_ratio;
                for replicate in 0..N_REPLICATES {
                    let seed = 0x1300_0000
                        ^ (bins as u64).rotate_left(7)
                        ^ (signal_count as u64).rotate_left(19)
                        ^ (noise_ratio as u64).rotate_left(31)
                        ^ replicate;
                    let mut rng = StdRng::seed_from_u64(seed);
                    let mut points = Vec::with_capacity(signal_count + noise_count);

                    for i in 0..signal_count {
                        let x = (i as f64 + 0.5) * SPAN / signal_count as f64;
                        let y = truth(x) + rng.random_range(-signal_jitter..signal_jitter);
                        points.push((LibraryRT(x), ObservedRTSeconds(y), 1.0));
                    }
                    for _ in 0..noise_count {
                        points.push((
                            LibraryRT(rng.random_range(0.0..SPAN)),
                            ObservedRTSeconds(rng.random_range(0.0..SPAN)),
                            1.0,
                        ));
                    }

                    let mut state =
                        CalibrationState::new(bins, (0.0, SPAN), (0.0, SPAN), lookback).unwrap();
                    state.update(points.into_iter()).unwrap();
                    state.fit();

                    let mut n = 0usize;
                    let mut sum_abs = 0.0;
                    let mut max_abs = 0.0f64;
                    if let Some(curve) = state.curve() {
                        for i in 0..N_SAMPLES {
                            let x = (i as f64 + 0.5) * SPAN / N_SAMPLES as f64;
                            if let Ok(y) = curve.predict(LibraryRT(x)) {
                                let error = (y.0 - truth(x)).abs();
                                n += 1;
                                sum_abs += error;
                                max_abs = max_abs.max(error);
                            }
                        }
                    }
                    let coverage = n as f64 / N_SAMPLES as f64;
                    let mean_abs = if n == 0 { f64::NAN } else { sum_abs / n as f64 };
                    let ridge = RidgeSummary::of(state.ridge_widths());

                    println!(
                        "{bins},{lookback},{signal_count},{noise_count},{replicate},{seed},{},{coverage},{mean_abs},{max_abs},{},{}",
                        state.path_indices().len(),
                        ridge.map_or(f64::NAN, |r| r.weighted_half_width),
                        ridge.map_or(f64::NAN, |r| r.in_ridge_ratio),
                    );
                }
            }
        }
    }
}
