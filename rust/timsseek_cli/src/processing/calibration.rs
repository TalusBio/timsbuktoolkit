use indicatif::ProgressIterator;
use timsquery::models::FlatIdx;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
};
use timsquery::traits::QueryGeom;
use timsquery::{
    MzMobilityStatsCollector,
    SpectralCollector,
    Tolerance,
};
use timsseek::data_sources::reference_library::ReferenceLibrary;
use timsseek::rt_calibration::{
    CALIBRANT_WEIGHT,
    CalibRtError,
    CalibratedGrid,
    CalibrationResult,
    DerivationParams,
    DimensionErrors,
    ErrorStats,
    LibraryRT,
    Point,
    RTCalibration,
    RidgeSummary,
    ridge_half_width_interp,
};
use timsseek::scoring::offsets::MzMobilityOffsets;
use timsseek::scoring::pipeline::Scorer;
use timsseek::scoring::timings::make_progress_bar;
use timsseek::scoring::{
    CALIBRANT_HEAP_GROUPS,
    CALIBRANT_RT_BANDS,
    CalibrantCandidate,
    CalibrationConfig,
    StratifiedCalibrantHeaps,
};
use timsseek::{
    IonAnnot,
    ScorerQueriable,
};
use tracing::{
    info,
    warn,
};

/// Lookup from `(precursor m/z * 100, precursor charge)` to the main speclib's
/// entries in that bucket: `(library RT seconds, sorted fragment m/z * 100)`.
/// Used to re-anchor a separate calibration library onto the main library's
/// iRT axis via shared-fragment matching.
pub(super) type PrecursorFragmentLookup =
    std::collections::HashMap<(i64, u8), Vec<(f32, Vec<i64>)>>;

pub(super) fn build_precursor_fragment_lookup(
    speclib: &ReferenceLibrary,
) -> PrecursorFragmentLookup {
    let mut map: PrecursorFragmentLookup = std::collections::HashMap::new();
    for item in speclib.iter() {
        let mz_key = (item.mono_precursor_mz() * 100.0).round() as i64;
        let charge = item.precursor_charge();
        let mut frag_mzs: Vec<i64> = item
            .iter_fragments_refs()
            .map(|(_, mz)| (mz * 100.0).round() as i64)
            .collect();
        frag_mzs.sort_unstable();
        map.entry((mz_key, charge))
            .or_default()
            .push((item.rt_seconds(), frag_mzs));
    }
    info!(
        "Built precursor+fragment lookup with {} unique (mz, charge) buckets from main speclib",
        map.len()
    );
    map
}

fn library_rt_for_candidate(
    candidate: &CalibrantCandidate,
    phase1_lib: &ReferenceLibrary,
    main_lookup: Option<&PrecursorFragmentLookup>,
) -> Option<f64> {
    let Some(lookup) = main_lookup else {
        return Some(candidate.library_rt.0 as f64);
    };
    let calib_item = phase1_lib.item_at(candidate.speclib_index);
    let mz_key = (calib_item.mono_precursor_mz() * 100.0).round() as i64;
    let charge = calib_item.precursor_charge();
    let bucket = lookup.get(&(mz_key, charge))?;

    let mut calib_frags: Vec<i64> = calib_item
        .iter_fragments_refs()
        .map(|(_, mz)| (mz * 100.0).round() as i64)
        .collect();
    calib_frags.sort_unstable();

    let calib_rt = calib_item.rt_seconds();
    bucket
        .iter()
        .filter_map(|(main_rt, main_frags)| {
            let shared = count_shared_fragments(&calib_frags, main_frags);
            (shared >= MIN_SHARED_FRAGMENTS).then_some((
                shared,
                (main_rt - calib_rt).abs(),
                *main_rt,
            ))
        })
        .min_by(|a, b| {
            b.0.cmp(&a.0).then(
                a.1.partial_cmp(&b.1)
                    .expect("NaN RT residual in calibrant matching"),
            )
        })
        .map(|(_, _, rt)| rt as f64)
}

/// Fixed so a run and its Phase-1 heap history can be reproduced exactly.
const PHASE1_SHUFFLE_SEED: u64 = 0x05EE_D130;

fn capacity_per_rt_band(config: &CalibrationConfig) -> usize {
    config.n_calibrants.div_ceil(CALIBRANT_RT_BANDS).max(1)
}

const CURVE_SAMPLES: usize = 64;
const MAX_MEAN_CURVE_DELTA_FRACTION: f32 = 0.01;
const MAX_P95_CURVE_DELTA_FRACTION: f32 = 0.03;
const MAX_POINT_CURVE_DELTA_FRACTION: f32 = 0.10;
const MAX_RIDGE_WIDTH_FRACTION: f32 = 0.25;
const REQUIRED_CONVERGENCE_PASSES: u8 = 3;

struct CurveSummary {
    visited: usize,
    predictions: [f32; CURVE_SAMPLES],
    ridge_width_fraction: f32,
}

#[derive(Debug, Clone, Copy)]
struct CurveDistance {
    mean: f32,
    p95: f32,
    max: f32,
    compared: usize,
}

impl CurveDistance {
    fn between(left: &CurveSummary, right: &CurveSummary) -> Self {
        let mut errors = [0.0; CURVE_SAMPLES];
        let mut compared = 0;
        for (left, right) in left.predictions.iter().zip(right.predictions) {
            if left.is_finite() && right.is_finite() {
                errors[compared] = (left - right).abs();
                compared += 1;
            }
        }
        if compared == 0 {
            return Self {
                mean: f32::INFINITY,
                p95: f32::INFINITY,
                max: f32::INFINITY,
                compared: 0,
            };
        }
        let errors = &mut errors[..compared];
        errors.sort_unstable_by(f32::total_cmp);
        Self {
            mean: errors.iter().sum::<f32>() / compared as f32,
            p95: errors[((compared - 1) as f32 * 0.95).round() as usize],
            max: errors[compared - 1],
            compared,
        }
    }

    fn passes(self, acquisition_duration: f32) -> bool {
        self.compared >= CURVE_SAMPLES / 2
            && self.mean <= acquisition_duration * MAX_MEAN_CURVE_DELTA_FRACTION
            && self.p95 <= acquisition_duration * MAX_P95_CURVE_DELTA_FRACTION
            && self.max <= acquisition_duration * MAX_POINT_CURVE_DELTA_FRACTION
    }
}

struct Phase1Convergence {
    states: [CalibratedGrid; CALIBRANT_HEAP_GROUPS + 1],
    history: Vec<CurveSummary>,
    consecutive_passes: u8,
    library_rt_range: (f32, f32),
    acquisition_duration: f32,
    grid_size: usize,
    mapped_library_rts: std::collections::HashMap<FlatIdx, Option<f64>>,
}

impl Phase1Convergence {
    fn new(
        config: &CalibrationConfig,
        library_rt_range: (f32, f32),
        acquisition_rt_range: (f64, f64),
    ) -> Result<Self, CalibRtError> {
        let states = [
            CalibratedGrid::deferred(config.grid_size, config.dp_lookback)?,
            CalibratedGrid::deferred(config.grid_size, config.dp_lookback)?,
            CalibratedGrid::deferred(config.grid_size, config.dp_lookback)?,
        ];
        Ok(Self {
            states,
            history: Vec::new(),
            consecutive_passes: 0,
            library_rt_range,
            acquisition_duration: (acquisition_rt_range.1 - acquisition_rt_range.0) as f32,
            grid_size: config.grid_size,
            mapped_library_rts: std::collections::HashMap::new(),
        })
    }

    fn summarize<'a>(
        state: &mut CalibratedGrid,
        candidates: impl Iterator<Item = &'a CalibrantCandidate> + Clone,
        visited: usize,
        library_rt_range: (f32, f32),
        acquisition_duration: f32,
        grid_size: usize,
        mapped_library_rts: &std::collections::HashMap<FlatIdx, Option<f64>>,
    ) -> Option<CurveSummary> {
        state
            .refit(
                grid_size,
                candidates.clone().filter_map(|candidate| {
                    mapped_library_rts
                        .get(&candidate.speclib_index)
                        .copied()
                        .flatten()
                        .map(|library_rt| (library_rt, candidate.apex_rt.0 as f64))
                }),
            )
            .ok()?;
        let curve = state.curve()?;
        let predictions = sample_curve(curve, library_rt_range);
        let ridge = RidgeSummary::of(state.ridge_widths())?;
        let ridge_width_fraction = (2.0 * ridge.weighted_half_width as f32) / acquisition_duration;
        Some(CurveSummary {
            visited,
            predictions,
            ridge_width_fraction,
        })
    }

    fn closest(&self, target: usize) -> Option<&CurveSummary> {
        self.history
            .iter()
            .min_by_key(|checkpoint| checkpoint.visited.abs_diff(target))
    }

    fn observe(
        &mut self,
        heaps: &StratifiedCalibrantHeaps,
        visited: usize,
        phase1_lib: &ReferenceLibrary,
        main_lookup: Option<&PrecursorFragmentLookup>,
    ) -> bool {
        for candidate in heaps.iter() {
            self.mapped_library_rts
                .entry(candidate.speclib_index)
                .or_insert_with(|| library_rt_for_candidate(candidate, phase1_lib, main_lookup));
        }
        let [state_a, state_b, state_union] = &mut self.states;
        let Some(curve_a) = Self::summarize(
            state_a,
            heaps.group_iter(0),
            visited,
            self.library_rt_range,
            self.acquisition_duration,
            self.grid_size,
            &self.mapped_library_rts,
        ) else {
            self.consecutive_passes = 0;
            return false;
        };
        let Some(curve_b) = Self::summarize(
            state_b,
            heaps.group_iter(1),
            visited,
            self.library_rt_range,
            self.acquisition_duration,
            self.grid_size,
            &self.mapped_library_rts,
        ) else {
            self.consecutive_passes = 0;
            return false;
        };
        let Some(curve_union) = Self::summarize(
            state_union,
            heaps.iter(),
            visited,
            self.library_rt_range,
            self.acquisition_duration,
            self.grid_size,
            &self.mapped_library_rts,
        ) else {
            self.consecutive_passes = 0;
            return false;
        };

        let cross_group = CurveDistance::between(&curve_a, &curve_b);
        let half = self.closest(visited / 2);
        let quarter = self.closest(visited / 4);
        let multiscale = half.zip(quarter).filter(|(half, quarter)| {
            half.visited != quarter.visited && half.visited != curve_union.visited
        });
        let passes = multiscale.is_some_and(|(half, quarter)| {
            cross_group.passes(self.acquisition_duration)
                && CurveDistance::between(&curve_union, half).passes(self.acquisition_duration)
                && CurveDistance::between(half, quarter).passes(self.acquisition_duration)
                && curve_a.ridge_width_fraction <= MAX_RIDGE_WIDTH_FRACTION
                && curve_b.ridge_width_fraction <= MAX_RIDGE_WIDTH_FRACTION
                && curve_union.ridge_width_fraction <= MAX_RIDGE_WIDTH_FRACTION
        });

        self.history.push(curve_union);
        self.consecutive_passes = if passes {
            self.consecutive_passes.saturating_add(1)
        } else {
            0
        };
        if self.consecutive_passes >= REQUIRED_CONVERGENCE_PASSES {
            info!(
                visited,
                mean_delta_seconds = cross_group.mean,
                p95_delta_seconds = cross_group.p95,
                max_delta_seconds = cross_group.max,
                "Phase 1 calibration converged"
            );
            true
        } else {
            false
        }
    }
}

fn sample_curve(curve: &RTCalibration, library_rt_range: (f32, f32)) -> [f32; CURVE_SAMPLES] {
    std::array::from_fn(|index| {
        let fraction = index as f64 / (CURVE_SAMPLES - 1) as f64;
        let library_rt =
            library_rt_range.0 as f64 + fraction * (library_rt_range.1 - library_rt_range.0) as f64;
        curve
            .predict(LibraryRT(library_rt))
            .map(|observed| observed.0 as f32)
            .unwrap_or(f32::NAN)
    })
}

/// Everything the pipeline below needs from the dev-only calibration fit
/// dashboard. `arm` is the only `#[cfg]`ed part: it has a real and a no-op
/// version, exposing the same items with the same signatures. Every call site
/// therefore calls unconditionally and `#[cfg]` stays off the search path.
pub(super) mod calib_dash_hook {
    use super::{
        CalibrantCandidate,
        CalibrationConfig,
        CalibrationResult,
    };

    pub use arm::*;

    #[cfg(feature = "calib-dashboard")]
    mod arm {
        use super::{
            CalibrantCandidate,
            CalibrationConfig,
            CalibrationResult,
        };
        use std::hash::{
            DefaultHasher,
            Hash,
            Hasher,
        };
        use std::io::IsTerminal;
        use timsquery::models::FlatIdx;

        // Dashboard-private identity for frame-to-frame visualization. Never
        // leaves this optional module or participates in search decisions.
        fn identity_hash(slot: FlatIdx) -> u64 {
            let mut hasher = DefaultHasher::new();
            slot.hash(&mut hasher);
            hasher.finish()
        }

        /// The dashboard for a whole run. Stays `None` inside unless
        /// `TIMSSEEK_CALIB_DASHBOARD` asks for it, so an ordinary run of a
        /// dashboard-enabled binary allocates nothing -- and neither does
        /// anything else in this module that keys off it.
        pub struct Dash {
            inner: Option<calib_dash::CalibDash>,
        }

        fn dashboard_requested() -> bool {
            let Some(raw) = std::env::var_os("TIMSSEEK_CALIB_DASHBOARD") else {
                return false;
            };
            match raw.to_string_lossy().trim().to_ascii_lowercase().as_str() {
                "" => false,
                "1" | "true" | "yes" => true,
                _ => {
                    tracing::warn!(
                        value = ?raw,
                        "TIMSSEEK_CALIB_DASHBOARD is not one of 1/true/yes; leaving the \
                         calibration dashboard off"
                    );
                    false
                }
            }
        }

        /// `n_chunks` is the chunk count of the library Phase 1 actually
        /// walks, which under `--calib-lib` is not the library the rest of the
        /// pipeline scores.
        pub fn attach(n_chunks: usize, config: &CalibrationConfig) -> Dash {
            let inner = dashboard_requested().then(|| {
                if !std::io::stdout().is_terminal() {
                    tracing::warn!(
                        "TIMSSEEK_CALIB_DASHBOARD is set but stdout is not a terminal; the \
                         dashboard will skip every pause and Phase 1 will run unattended"
                    );
                }
                let budget_bytes = match std::env::var_os("CALIB_DASH_FRAME_BUDGET_MB") {
                    None => calib_dash::DEFAULT_RUN_BUDGET_BYTES,
                    Some(raw) => match raw
                        .to_str()
                        .and_then(|v| v.parse::<usize>().ok())
                        .filter(|mb| *mb > 0)
                    {
                        Some(mb) => mb.saturating_mul(1024 * 1024),
                        None => {
                            tracing::warn!(
                                value = ?raw,
                                "CALIB_DASH_FRAME_BUDGET_MB is not a positive whole number \
                                 of megabytes; falling back to the default budget"
                            );
                            calib_dash::DEFAULT_RUN_BUDGET_BYTES
                        }
                    },
                };
                calib_dash::CalibDash::new(
                    n_chunks,
                    config.n_calibrants,
                    config.grid_size,
                    config.dp_lookback,
                    budget_bytes,
                )
            });
            Dash { inner }
        }

        /// Exits the process on `Ctrl-C` at a pause: raw mode swallows `SIGINT`,
        /// and returning to the pipeline would only run Phases 2 and 3 against a
        /// truncated calibrant heap.
        pub fn on_batch<'a>(
            dash: &mut Dash,
            chunk: usize,
            calibrants: impl Iterator<Item = &'a CalibrantCandidate>,
        ) {
            let Some(d) = dash.inner.as_mut() else {
                return;
            };
            let flow = d.on_batch(
                chunk,
                calibrants.map(|c| calib_dash::CalibrantPoint {
                    library_rt: c.library_rt.0 as f64,
                    observed_rt: c.apex_rt.0 as f64,
                    identity: identity_hash(c.speclib_index),
                }),
            );
            if matches!(flow, calib_dash::Flow::Abort) {
                tracing::warn!("calibration dashboard aborted Phase 1 at chunk {chunk}; exiting");
                std::process::exit(130);
            }
        }

        pub fn finish(dash: &mut Dash) {
            if let Some(d) = dash.inner.as_mut() {
                d.finish();
            }
        }

        /// Wires the Phase 2 fit into the Tolerances tab, then pauses so the
        /// tabs and the batch scrubber are reachable before Phase 3 starts.
        pub fn show_final(dash: &mut Dash, calibration: &CalibrationResult) {
            let Some(d) = dash.inner.as_mut() else {
                return;
            };
            let mobility = calibration.mobility_tolerance();
            let rt_tolerance_seconds = calibration.rt_tolerance_minutes() as f64 * 60.0;
            d.show_final(
                calib_dash::FitRecording::from_state(calibration.state()),
                calib_dash::ToleranceSummary {
                    mz_ppm: calibration.mz_tolerance(),
                    mobility_pct: (mobility.0 as f64, mobility.1 as f64),
                    rt_seconds: rt_tolerance_seconds,
                    n_calibrants: calibration.errors().rt_seconds.n,
                },
            );
            d.present();
        }
    }

    #[cfg(not(feature = "calib-dashboard"))]
    mod arm {
        use super::{
            CalibrantCandidate,
            CalibrationConfig,
            CalibrationResult,
        };

        pub struct Dash;

        pub fn attach(_n_chunks: usize, _config: &CalibrationConfig) -> Dash {
            Dash
        }

        pub fn on_batch<'a>(
            _dash: &mut Dash,
            _chunk: usize,
            _calibrants: impl Iterator<Item = &'a CalibrantCandidate>,
        ) {
        }

        pub fn finish(_dash: &mut Dash) {}

        pub fn show_final(_dash: &mut Dash, _calibration: &CalibrationResult) {}
    }
}

/// Check that two speclibs are on a compatible RT scale.
/// Warns loudly if the RT ranges don't overlap, which would produce a useless calibration.
pub(super) fn check_rt_scale_compatibility(
    main_lib: &ReferenceLibrary,
    calib_lib: &ReferenceLibrary,
) {
    fn rt_range(lib: &ReferenceLibrary) -> (f32, f32) {
        let mut min_rt = f32::INFINITY;
        let mut max_rt = f32::NEG_INFINITY;
        for q in lib.iter() {
            let rt = q.rt_seconds();
            min_rt = min_rt.min(rt);
            max_rt = max_rt.max(rt);
        }
        (min_rt, max_rt)
    }

    let (main_min, main_max) = rt_range(main_lib);
    let (calib_min, calib_max) = rt_range(calib_lib);

    info!(
        "RT ranges -- main speclib: [{:.1}, {:.1}]s, calib lib: [{:.1}, {:.1}]s",
        main_min, main_max, calib_min, calib_max
    );

    // Check overlap
    let overlap_start = main_min.max(calib_min);
    let overlap_end = main_max.min(calib_max);

    if overlap_start >= overlap_end {
        warn!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        warn!("!! RT SCALE MISMATCH: main speclib and calibration library  !!");
        warn!("!! have NO overlapping RT range. The calibration will be    !!");
        warn!("!! meaningless. Ensure both libraries use the same iRT      !!");
        warn!("!! scale (e.g., both from the same prediction model).       !!");
        warn!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        return;
    }

    let main_span = main_max - main_min;
    let overlap_span = overlap_end - overlap_start;
    let overlap_pct = overlap_span / main_span * 100.0;

    if overlap_pct < 50.0 {
        warn!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        warn!(
            "!! RT SCALE WARNING: only {:.0}% overlap between main speclib !!",
            overlap_pct
        );
        warn!("!! and calibration library. Calibration may be unreliable.  !!");
        warn!("!! Ensure both libraries use the same iRT scale.            !!");
        warn!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    } else if overlap_pct < 80.0 {
        warn!(
            "RT overlap between main speclib and calib lib is {:.0}% -- may affect calibration at the extremes",
            overlap_pct
        );
    }
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub(super) fn phase1_prescore<I: ScorerQueriable>(
    speclib: &ReferenceLibrary,
    pipeline: &Scorer<I>,
    chunk_size: usize,
    config: &CalibrationConfig,
    main_lookup: Option<&PrecursorFragmentLookup>,
    mapped_rt_range: Option<(f32, f32)>,
) -> (
    Vec<CalibrantCandidate>,
    timsseek::scoring::PrescoreTimings,
    calib_dash_hook::Dash,
) {
    let total = speclib.len();
    let n_chunks = total.div_ceil(chunk_size);
    let pb = make_progress_bar(n_chunks as u64, "Phase 1");

    let mut dash = calib_dash_hook::attach(n_chunks, config);

    let mut timings = timsseek::scoring::PrescoreTimings::default();
    let shuffled = speclib.shuffled_flats(PHASE1_SHUFFLE_SEED);
    let rt_range = speclib.rt_range().unwrap_or((0.0, 0.0));
    let convergence_rt_range = mapped_rt_range.unwrap_or(rt_range);
    let capacity_per_band = capacity_per_rt_band(config);
    let acquisition_rt_range = pipeline.acquisition_rt_range_seconds();
    let mut heaps = StratifiedCalibrantHeaps::new(capacity_per_band, rt_range);
    let mut convergence =
        match Phase1Convergence::new(config, convergence_rt_range, acquisition_rt_range) {
            Ok(convergence) => Some(convergence),
            Err(error) => {
                warn!(
                    ?error,
                    "could not initialize Phase-1 convergence; running exhaustively"
                );
                None
            }
        };
    let mut visited = 0usize;

    // Shuffle individual flat indices, not source-order chunks: libraries are
    // commonly ordered by RT, sequence, or decoy variant, and a future prefix
    // stop must not inherit any of those structures.
    for (batch_idx, batch) in shuffled.chunks(chunk_size).enumerate().progress_with(pb) {
        let batch_heaps = pipeline.prescore_partitioned_batch(
            speclib,
            batch,
            capacity_per_band,
            rt_range,
            &mut timings,
        );
        heaps.merge_in_place(batch_heaps);
        visited += batch.len();

        calib_dash_hook::on_batch(&mut dash, batch_idx, heaps.iter());
        if visited < total
            && convergence
                .as_mut()
                .is_some_and(|tracker| tracker.observe(&heaps, visited, speclib, main_lookup))
        {
            break;
        }
    }
    calib_dash_hook::finish(&mut dash);

    (heaps.into_vec(), timings, dash)
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
/// Count shared fragment m/z values between two sorted lists (within 0.01 Da = 1 unit of quantized m/z).
fn count_shared_fragments(a: &[i64], b: &[i64]) -> usize {
    let mut i = 0;
    let mut j = 0;
    let mut count = 0;
    while i < a.len() && j < b.len() {
        let diff = a[i] - b[j];
        if diff.abs() <= 1 {
            count += 1;
            i += 1;
            j += 1;
        } else if diff < 0 {
            i += 1;
        } else {
            j += 1;
        }
    }
    count
}

const MIN_SHARED_FRAGMENTS: usize = 5;

/// Every calibrant was found inside the first-pass window, so a derived window
/// wider than it measures contamination of the per-calibrant offset (a mean
/// over a 50 ppm query box), not the instrument. On a 21-min HeLa run the
/// unclamped derivation gave 60 ppm; at 10-15 ppm Phase 3 found 34% more
/// peptides. Returns the bound when it clamped.
fn clamp_mz_window_to_first_pass(
    mz_ppm: &mut (f64, f64),
    first_pass: &Tolerance,
) -> Option<(f64, f64)> {
    let MzTolerance::Ppm((lo, hi)) = first_pass.ms else {
        return None;
    };
    if mz_ppm.0 <= lo && mz_ppm.1 <= hi {
        return None;
    }
    *mz_ppm = (mz_ppm.0.min(lo), mz_ppm.1.min(hi));
    Some((lo, hi))
}

#[cfg(test)]
mod clamp_tests {
    use super::*;

    fn first_pass(lo: f64, hi: f64) -> Tolerance {
        Tolerance {
            ms: MzTolerance::Ppm((lo, hi)),
            rt: RtTolerance::Unrestricted,
            mobility: MobilityTolerance::Unrestricted,
            quad: QuadTolerance::Absolute((0.1, 0.1)),
        }
    }

    #[test]
    fn clamps_each_side_to_the_first_pass_window() {
        let mut w = (59.1, 61.3);
        assert_eq!(
            clamp_mz_window_to_first_pass(&mut w, &first_pass(15.0, 15.0)),
            Some((15.0, 15.0))
        );
        assert_eq!(w, (15.0, 15.0));

        let mut w = (8.0, 20.0);
        assert!(clamp_mz_window_to_first_pass(&mut w, &first_pass(15.0, 15.0)).is_some());
        assert_eq!(w, (8.0, 15.0));
    }

    #[test]
    fn inside_the_first_pass_window_is_untouched() {
        let mut w = (8.9, 11.1);
        assert_eq!(
            clamp_mz_window_to_first_pass(&mut w, &first_pass(15.0, 15.0)),
            None
        );
        assert_eq!(w, (8.9, 11.1));
    }
}

pub(super) fn calibrate_from_phase1<I: ScorerQueriable>(
    candidates: Vec<CalibrantCandidate>,
    phase1_lib: &ReferenceLibrary,
    main_lookup: Option<&PrecursorFragmentLookup>,
    pipeline: &Scorer<I>,
    config: &CalibrationConfig,
) -> Result<CalibrationResult, CalibRtError> {
    // === Step A: Fit iRT -> RT curve ===
    // With a separate calib lib, the curve's x-axis is the main speclib's iRT
    // (matched via shared fragments). We keep this per-candidate so later
    // stages (ridge filter, residual stats) use the same x-axis the fit saw.
    let library_rt_for_candidate: Vec<Option<f64>> = candidates
        .iter()
        .map(|candidate| library_rt_for_candidate(candidate, phase1_lib, main_lookup))
        .collect();

    let points: Vec<Point> = candidates
        .iter()
        .zip(library_rt_for_candidate.iter())
        .filter_map(|(c, lib_rt)| {
            let lib_rt = (*lib_rt)?;
            Some(Point {
                library: lib_rt,
                observed: c.apex_rt.0 as f64,
                weight: CALIBRANT_WEIGHT,
            })
        })
        .collect();

    if main_lookup.is_some() {
        let matched = points.len();
        let total = candidates.len();
        let rate = if total > 0 {
            matched as f64 / total as f64
        } else {
            0.0
        };
        // Surface a collapsing cross-library match loudly -- otherwise it silently
        // becomes a ZeroRange grid -> identity fallback -> ~0 IDs, with no hint why.
        if matched < 2 || rate < 0.10 {
            warn!(
                "Calibration: only {}/{} calibrants ({:.1}%) matched the main speclib (>= {} shared fragments). \
                 Calibration will fail or degrade (likely ZeroRange -> fallback). Common causes: main speclib \
                 entries carry < {} usable fragments (e.g. over-aggressive fragment filtering on load), or the \
                 calib lib and main speclib m/z scales don't align.",
                matched,
                total,
                rate * 100.0,
                MIN_SHARED_FRAGMENTS,
                MIN_SHARED_FRAGMENTS,
            );
        } else {
            info!(
                "Calibration: {} of {} calibrants matched in main speclib (>={} shared fragments)",
                matched, total, MIN_SHARED_FRAGMENTS,
            );
        }
    }

    // Use CalibrationState for fitting + ridge width measurement
    let mut cal_state = CalibratedGrid::deferred(config.grid_size, config.dp_lookback)?;
    cal_state.refit(
        config.grid_size,
        points.iter().map(|p| (p.library, p.observed)),
    )?;
    let cal_curve = cal_state.curve().ok_or(CalibRtError::NoPoints)?;

    // Position-dependent RT tolerance comes from the ridge the fit measured.
    let ridge_widths = cal_state.ridge_widths();
    if let Some(s) = RidgeSummary::of(ridge_widths) {
        info!(
            "Ridge width: weighted avg {:.1}s across {} columns (min {:.1}s, max {:.1}s)",
            s.weighted_half_width, s.n_columns, s.min_half_width, s.max_half_width,
        );
    }

    // === Step B: Measure m/z and mobility errors at calibrant apexes ===
    let query_tolerance = Tolerance {
        ms: MzTolerance::Ppm((50.0, 50.0)),
        rt: RtTolerance::Minutes((
            config.calibration_query_rt_window_minutes,
            config.calibration_query_rt_window_minutes,
        )),
        mobility: MobilityTolerance::Pct((5.0, 5.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1)),
    };

    let mut mz_errors_ppm: Vec<f32> = Vec::with_capacity(candidates.len());
    let mut mobility_errors_pct: Vec<f32> = Vec::with_capacity(candidates.len());
    let mut rt_residuals_seconds: Vec<f32> = Vec::with_capacity(candidates.len());

    let mut n_off_ridge = 0usize;
    for (candidate, library_rt_opt) in candidates.iter().zip(library_rt_for_candidate.iter()) {
        let item = phase1_lib.item_at(candidate.speclib_index);

        // Skip calibrants that never matched in the main speclib.
        let Some(library_rt_s) = *library_rt_opt else {
            continue;
        };

        let predicted_rt = match cal_curve.predict(LibraryRT(library_rt_s)) {
            Ok(o) => o.0,
            Err(_) => continue,
        };
        let rt_residual_signed = candidate.apex_rt.0 as f64 - predicted_rt;
        let half_width = ridge_half_width_interp(ridge_widths, library_rt_s);
        let in_ridge = match half_width {
            Some(hw) => rt_residual_signed.abs() <= hw,
            None => true,
        };
        if !in_ridge {
            n_off_ridge += 1;
            continue;
        }

        // No clone: build collector from the flyweight (QueryGeom) + rt override.
        let mut agg: SpectralCollector<IonAnnot, MzMobilityStatsCollector> =
            SpectralCollector::new(&item);
        agg.reset_with_overrides(&item, Some(candidate.apex_rt.0), None);
        pipeline.index.add_query(&mut agg, &query_tolerance);

        let expected_mob = item.mobility_ook0() as f64;
        let offsets = MzMobilityOffsets::new(&agg, expected_mob);
        let Some((mz_err_ppm, mob_err_pct)) = offsets.weighted_ms1() else {
            continue;
        };

        mz_errors_ppm.push(mz_err_ppm);
        // Mobility is NaN for a non-searchable-axis run (mzML/no-IM lib); don't
        // let it poison the MAD-based mobility tolerance (unused there anyway).
        if mob_err_pct.is_finite() {
            mobility_errors_pct.push(mob_err_pct);
        }
        rt_residuals_seconds.push(rt_residual_signed as f32);
    }

    info!(
        "Ridge filter: kept {}/{} calibrants (dropped {} off-ridge)",
        mz_errors_ppm.len(),
        candidates.len(),
        n_off_ridge,
    );

    // === Step C: Derive tolerances from error distributions ===
    let errors = DimensionErrors {
        mz_ppm: ErrorStats::from_slice(&mz_errors_ppm),
        mobility_pct: ErrorStats::from_slice(&mobility_errors_pct),
        rt_seconds: ErrorStats::from_slice(&rt_residuals_seconds),
    };
    let mut derivation = DerivationParams::default();
    derivation.sigma.mz = config.mz_sigma;
    derivation.sigma.mobility = config.mobility_sigma;
    derivation.sigma.rt = config.rt_sigma_factor;
    derivation.floors.rt_minutes = config.min_rt_tolerance_minutes;

    let mut windows = errors.derive_windows(&derivation);
    if let Some((lo, hi)) =
        clamp_mz_window_to_first_pass(&mut windows.mz_ppm, &pipeline.broad_tolerance)
    {
        warn!(
            "derived m/z window exceeded the first-pass window ({lo:.1}, {hi:.1}) ppm; clamped. \
             Residual MAD {:.1} ppm is contamination of the box-mean offset, not the instrument.",
            errors.mz_ppm.mad
        );
    }
    info!(
        "RT residuals: MAD={:.1}s, n={}",
        errors.rt_seconds.mad, errors.rt_seconds.n
    );
    info!(
        "Calibration: RT tol={:.2} min, m/z tol=({:.1}, {:.1}) ppm, mob tol=({:.1}, {:.1}) %",
        windows.rt_minutes,
        windows.mz_ppm.0,
        windows.mz_ppm.1,
        windows.mobility_pct.0,
        windows.mobility_pct.1,
    );

    Ok(CalibrationResult::new(
        cal_state,
        windows.rt_minutes,
        windows.mz_ppm,
        windows.mobility_pct,
    )?
    .with_error_stats(errors)
    .with_derivation(derivation))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn curve_summary(visited: usize, value: f32) -> CurveSummary {
        CurveSummary {
            visited,
            predictions: [value; CURVE_SAMPLES],
            ridge_width_fraction: 0.1,
        }
    }

    #[test]
    fn curve_distance_uses_classic_error_metrics() {
        let left = curve_summary(20, 10.0);
        let mut right = curve_summary(40, 11.0);
        right.predictions[CURVE_SAMPLES - 1] = 20.0;

        let distance = CurveDistance::between(&left, &right);
        assert_eq!(distance.compared, CURVE_SAMPLES);
        assert!((distance.mean - 1.140625).abs() < f32::EPSILON);
        assert_eq!(distance.p95, 1.0);
        assert_eq!(distance.max, 10.0);
    }

    #[test]
    fn curve_distance_requires_half_the_sampling_domain() {
        let left = curve_summary(20, 10.0);
        let mut right = curve_summary(40, f32::NAN);
        right.predictions[..CURVE_SAMPLES / 2 - 1].fill(10.0);

        assert!(!CurveDistance::between(&left, &right).passes(1_000.0));
    }
}
