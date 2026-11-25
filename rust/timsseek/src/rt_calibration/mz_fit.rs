use super::ErrorStats;
use serde::{
    Deserialize,
    Serialize,
};

pub(super) const MZ_SIGMA_MULTIPLIER: f64 = 3.0;
const MAX_ITERATIONS: u16 = 200;
const MIN_SAMPLES: usize = 32;
const CONVERGENCE_TOLERANCE: f64 = 1e-5;
const MAX_BACKGROUND_FRACTION: f64 = 0.98;
const SQRT_TWO_PI: f64 = 2.506_628_274_631_000_2;

/// Gaussian center and standard deviation, both in ppm.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MzFit {
    pub mu: f64,
    pub sigma: f64,
    pub background_fraction: f64,
    pub iterations: u16,
    pub converged: bool,
    pub n: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MzFitFailure {
    NoSharedPpmWindow,
    InvalidWindow,
    TooFewSamples,
    CollapsedSigma,
    BackgroundDominated,
    GaussianSpansWindow,
    NoEvidenceAgainstUniform,
    WindowDoesNotOverlap,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum MzFitOutcome {
    Fitted(MzFit),
    Fallback {
        reason: MzFitFailure,
        n: usize,
        #[serde(skip_serializing_if = "Option::is_none")]
        estimate: Option<MzFit>,
    },
}

impl MzFitOutcome {
    pub(super) fn fallback(reason: MzFitFailure, n: usize, estimate: Option<MzFit>) -> Self {
        Self::Fallback {
            reason,
            n,
            estimate,
        }
    }
}

pub(super) fn fit_mixture(samples: &[f32], stats: &ErrorStats, width: f64) -> MzFitOutcome {
    assert!(width.is_finite() && width > 0.0);
    let n = samples.len();
    if n < MIN_SAMPLES {
        return MzFitOutcome::fallback(MzFitFailure::TooFewSamples, n, None);
    }
    if samples.iter().all(|x| *x == samples[0]) {
        return MzFitOutcome::fallback(MzFitFailure::CollapsedSigma, n, None);
    }

    // Scale densities by the window width so the uniform density is one.
    // The sigma bound only guards the Gaussian density calculation.
    let min_sigma = (width * 1e-9).max(f64::MIN_POSITIVE);
    let mut fit = MzFit {
        mu: f64::from(stats.median),
        sigma: (1.4826 * f64::from(stats.mad)).max(min_sigma),
        background_fraction: 0.5,
        iterations: 0,
        converged: false,
        n,
    };
    assert!(fit.mu.is_finite());
    assert!(fit.sigma.is_finite());
    for iteration in 1..=MAX_ITERATIONS {
        let (signal_scale, background) = scaled_densities(&fit, width);
        let inverse_sigma = fit.sigma.recip();
        let mut signal_weight = 0.0;
        let mut background_weight = 0.0;
        let mut weighted_delta = 0.0;
        let mut weighted_delta_sq = 0.0;
        for &sample in samples {
            let delta = f64::from(sample) - fit.mu;
            let signal = scaled_gaussian(delta, inverse_sigma, signal_scale);
            let density = signal + background;
            let weight = signal / density;
            signal_weight += weight;
            background_weight += background / density;
            weighted_delta += weight * delta;
            weighted_delta_sq += weight * delta * delta;
        }
        if signal_weight <= f64::EPSILON {
            return MzFitOutcome::fallback(MzFitFailure::BackgroundDominated, n, Some(fit));
        }
        let shift = weighted_delta / signal_weight;
        let sigma = (weighted_delta_sq / signal_weight - shift * shift)
            .max(0.0)
            .sqrt()
            .max(min_sigma);
        let mu = fit.mu + shift;
        let background_fraction = background_weight / (signal_weight + background_weight);
        assert!(mu.is_finite());
        assert!(sigma.is_finite());
        assert!(background_fraction.is_finite());
        let converged = shift.abs() <= CONVERGENCE_TOLERANCE * fit.sigma
            && (sigma - fit.sigma).abs() <= CONVERGENCE_TOLERANCE * fit.sigma
            && (background_fraction - fit.background_fraction).abs() <= CONVERGENCE_TOLERANCE;
        fit = MzFit {
            mu,
            sigma,
            background_fraction,
            iterations: iteration,
            converged,
            n,
        };
        if converged {
            break;
        }
    }

    let reason = if fit.background_fraction >= MAX_BACKGROUND_FRACTION {
        MzFitFailure::BackgroundDominated
    } else if fit.sigma <= min_sigma {
        MzFitFailure::CollapsedSigma
    } else if 2.0 * MZ_SIGMA_MULTIPLIER * fit.sigma >= width {
        MzFitFailure::GaussianSpansWindow
    } else if !beats_uniform(samples, &fit, width) {
        MzFitFailure::NoEvidenceAgainstUniform
    } else {
        return MzFitOutcome::Fitted(fit);
    };
    MzFitOutcome::fallback(reason, n, Some(fit))
}

fn scaled_densities(fit: &MzFit, width: f64) -> (f64, f64) {
    let background = fit
        .background_fraction
        .clamp(f64::EPSILON, 1.0 - f64::EPSILON);
    (
        (1.0 - background) * (width / fit.sigma) / SQRT_TWO_PI,
        background,
    )
}

fn scaled_gaussian(delta: f64, inverse_sigma: f64, scale: f64) -> f64 {
    scale * (-0.5 * (delta * inverse_sigma).powi(2)).exp()
}

fn beats_uniform(samples: &[f32], fit: &MzFit, width: f64) -> bool {
    let (signal_scale, background) = scaled_densities(fit, width);
    let inverse_sigma = fit.sigma.recip();
    let log_likelihood_gain = samples
        .iter()
        .map(|&x| {
            (scaled_gaussian(f64::from(x) - fit.mu, inverse_sigma, signal_scale) + background).ln()
        })
        .sum::<f64>();
    // Require a BIC improvement over the known uniform: the mixture has three
    // free parameters and can otherwise fit a chance bump in uniform noise.
    2.0 * log_likelihood_gain > 3.0 * (fit.n as f64).ln()
}

#[cfg(test)]
pub(super) fn synthetic_samples(
    n: usize,
    mu: f64,
    sigma: f64,
    background: f64,
    window: (f64, f64),
    seed: u64,
) -> Vec<f32> {
    use rand::rngs::StdRng;
    use rand::{
        Rng,
        SeedableRng,
    };

    let mut rng = StdRng::seed_from_u64(seed);
    let mut samples = Vec::with_capacity(n);
    for _ in 0..n {
        let x = if rng.random::<f64>() < background {
            rng.random_range(window.0..window.1)
        } else {
            let r = (-2.0 * (1.0 - rng.random::<f64>()).ln()).sqrt();
            let angle = rng.random_range(0.0..std::f64::consts::TAU);
            mu + sigma * r * angle.cos()
        };
        assert!((window.0..=window.1).contains(&x));
        samples.push(x as f32);
    }
    samples
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_gaussian_center_spread_and_background() {
        let (mu, sigma) = (1.25, 1.6);
        let window = (-12.0, 18.0);
        for (seed, background) in [0.0, 0.25, 0.5, 0.8].into_iter().enumerate() {
            let samples = synthetic_samples(30_000, mu, sigma, background, window, seed as u64);
            let result = fit_mixture(&samples, &ErrorStats::from_slice(&samples), 30.0);
            let MzFitOutcome::Fitted(fit) = result else {
                panic!("background {background}: {result:?}");
            };
            assert!((fit.mu - mu).abs() < 0.12, "{fit:?}");
            assert!((fit.sigma - sigma).abs() < 0.12, "{fit:?}");
            assert!(
                (fit.background_fraction - background).abs() < 0.025,
                "{fit:?}"
            );
            assert_eq!(fit.n, samples.len());
            assert!(fit.converged, "{fit:?}");
            assert!(fit.iterations < MAX_ITERATIONS);
        }
    }

    #[test]
    fn rejects_too_few_samples_and_a_collapsed_core() {
        for (samples, expected) in [
            (&[1.0; MIN_SAMPLES - 1][..], MzFitFailure::TooFewSamples),
            (&[1.0; MIN_SAMPLES][..], MzFitFailure::CollapsedSigma),
        ] {
            let result = fit_mixture(samples, &ErrorStats::from_slice(samples), 30.0);
            assert_eq!(
                result,
                MzFitOutcome::fallback(expected, samples.len(), None)
            );
        }
    }

    #[test]
    fn uniform_background_falls_back_at_the_iteration_cap() {
        let samples = synthetic_samples(20_000, 0.0, 1.0, 1.0, (-15.0, 15.0), 4);
        let result = fit_mixture(&samples, &ErrorStats::from_slice(&samples), 30.0);
        let MzFitOutcome::Fallback {
            estimate: Some(fit),
            ..
        } = result
        else {
            panic!("{result:?}");
        };
        assert_eq!(fit.iterations, MAX_ITERATIONS);
        assert!(!fit.converged);
    }
}
