use super::config::OutputConfig;
use indicatif::ProgressIterator;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
};
use timsquery::traits::QueryGeom;
use timsquery::{
    IndexedTimstofPeaks,
    MzMobilityStatsCollector,
    SpectralCollector,
    Tolerance,
};
use timsseek::data_sources::reference_library::ReferenceLibrary;
use timsseek::errors::TimsSeekError;
use timsseek::ml::qvalues::{
    ThresholdCounts,
    report_qvalues_at_thresholds,
};
use timsseek::ml::{
    RescoreFeatureStats,
    RescoreModel,
    rescore_with,
};
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
    RidgeSummary,
    ridge_half_width_interp,
};
use timsseek::scoring::offsets::MzMobilityOffsets;
use timsseek::scoring::pipeline::Scorer;
use timsseek::scoring::timings::{
    TimedStep,
    make_progress_bar,
};
use timsseek::scoring::{
    CalibrantCandidate,
    CalibrantHeap,
    CalibrationConfig,
    CompetedCandidate,
    PipelineReport,
    ScoreTimings,
    ScoredCandidate,
    SkipCounts,
};
use timsseek::{
    IonAnnot,
    ScorerQueriable,
};
use tracing::{
    debug,
    info,
    warn,
};

/// Lookup from `(precursor m/z * 100, precursor charge)` to the main speclib's
/// entries in that bucket: `(library RT seconds, sorted fragment m/z * 100)`.
/// Used to re-anchor a separate calibration library onto the main library's
/// iRT axis via shared-fragment matching.
type PrecursorFragmentLookup = std::collections::HashMap<(i64, u8), Vec<(f32, Vec<i64>)>>;

/// Everything the pipeline below needs from the dev-only calibration fit
/// dashboard. `arm` is the only `#[cfg]`ed part: it has a real and a no-op
/// version, exposing the same items with the same signatures. Every call site
/// therefore calls unconditionally and `#[cfg]` stays off the search path.
mod calib_dash_hook {
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

        /// The dashboard compares calibrant identity across batches but never
        /// displays it, so it takes a hash of the arena handle rather than the
        /// id that handle resolves to. That keeps its point type `Copy` and,
        /// more importantly, keeps `size_of` the true cost of a point -- which is
        /// what its replay memory budget is computed from. Only meaningful
        /// within one process.
        ///
        /// Concrete in `FlatIdx` rather than `impl Hash`, which would also
        /// accept `score.to_bits()` -- compiles, means nothing.
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
fn check_rt_scale_compatibility(main_lib: &ReferenceLibrary, calib_lib: &ReferenceLibrary) {
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

pub const FEATURE_STATS_FILENAME: &str = "results.feature_stats.tsv";
pub const FEATURE_IMPORTANCE_FILENAME: &str = "results.feature_importance.tsv";

fn write_tsv(
    path: &std::path::Path,
    header: &str,
    mut write_rows: impl FnMut(&mut String),
    label: &str,
) -> std::io::Result<()> {
    use std::fmt::Write as _;
    let mut buf = String::new();
    writeln!(buf, "{header}").unwrap();
    write_rows(&mut buf);
    std::fs::write(path, buf)?;
    eprintln!("wrote {label}: {}", path.display());
    tracing::info!(path = %path.display(), "wrote {} tsv", label);
    Ok(())
}

fn write_feature_stats_sidecar(
    stats: &RescoreFeatureStats,
    parquet_path: &std::path::Path,
) -> std::io::Result<()> {
    use std::fmt::Write as _;

    write_tsv(
        &parquet_path.with_file_name(FEATURE_STATS_FILENAME),
        "name\tmean\tmissing\tfold",
        |buf| {
            for fold in stats.iter() {
                for fs in fold.feature_stats.iter() {
                    writeln!(
                        buf,
                        "{}\t{}\t{}\t{}",
                        fs.name, fs.mean, fs.nan_ratio, fold.fold
                    )
                    .unwrap();
                }
            }
        },
        "feature stats",
    )?;

    write_tsv(
        &parquet_path.with_file_name(FEATURE_IMPORTANCE_FILENAME),
        "name\tgain\tfold",
        |buf| {
            for fold in stats.iter() {
                for (name, gain) in fold.feature_importance.iter() {
                    writeln!(buf, "{}\t{}\t{}", name, gain, fold.fold).unwrap();
                }
            }
        },
        "feature importance",
    )?;

    Ok(())
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn execute_pipeline<I: ScorerQueriable>(
    speclib: &ReferenceLibrary,
    calib_lib: Option<&ReferenceLibrary>,
    pipeline: &Scorer<I>,
    chunk_size: usize,
    out_path: &OutputConfig,
    max_qvalue: f32,
    calib_config: &CalibrationConfig,
    no_feature_stats: bool,
    rescore_model: RescoreModel,
) -> std::result::Result<PipelineReport, TimsSeekError> {
    // === PHASE 1: Broad prescore -> collect top calibrants ===
    // Use calibration library if provided, otherwise fall back to main speclib
    let phase1_lib = calib_lib.unwrap_or(speclib);
    if let Some(clib) = calib_lib {
        info!(
            "Phase 1: Broad prescore using calibration library ({} entries)...",
            clib.len()
        );
        check_rt_scale_compatibility(speclib, clib);
    } else {
        info!("Phase 1: Broad prescore (unrestricted RT)...");
    }
    // Reset the `for_each_peak` filter-funnel counters so the per-phase
    // dump reflects just this phase's work. Both calls no-op when the
    // `query-instr` feature is disabled.
    timscentroid::indexing::reset_for_each_peak_funnel();
    let step = TimedStep::begin("Phase 1: Prescore");
    let (calibrants, phase1_timings, mut calib_dash_state) =
        phase1_prescore(phase1_lib, pipeline, chunk_size, calib_config);
    let phase1_ms = step
        .finish_with(format_args!(
            "{} calibrants/ {} candidates",
            calibrants.len(),
            phase1_timings.n_scored
        ))
        .as_millis() as u64;
    alloc_track::snap!("Phase 1: Prescore");
    timscentroid::indexing::dump_for_each_peak_funnel("Phase 1");
    timscentroid::indexing::reset_for_each_peak_funnel();
    info!(
        "Phase 1 detail: extraction {:?}, scoring {:?}, {} passed filter, {} scored",
        phase1_timings.extraction,
        phase1_timings.scoring,
        phase1_timings.n_passed_filter,
        phase1_timings.n_scored,
    );

    // === PHASE 2: Calibration (fit RT + measure errors + derive tolerances) ===
    // Build lookup from main speclib when using a separate calib lib.
    // Maps (quantized_precursor_mz, charge) -> Vec<(rt, sorted_fragment_mzs)>.
    // Matching requires same precursor (0.01 Da) + charge + at least 5 shared fragment masses.
    let main_lookup: Option<PrecursorFragmentLookup> = if calib_lib.is_some() {
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
        Some(map)
    } else {
        None
    };

    info!("Phase 2: Calibration...");
    let step = TimedStep::begin("Phase 2: Calibrate");
    let calibration = match calibrate_from_phase1(
        calibrants,
        phase1_lib,
        main_lookup.as_ref(),
        pipeline,
        calib_config,
    ) {
        Ok(calib) => {
            info!("Calibration succeeded");
            calib
        }
        Err(e) => {
            tracing::error!("Calibration failed: {:?}. Using fallback.", e);
            CalibrationResult::fallback(pipeline)
        }
    };
    let phase2_ms = step
        .finish_with(format_args!(
            "{} fit points → {} path nodes",
            calibration.state().fit_points().len(),
            calibration.ridge_width_summary().map_or(0, |s| s.n_columns),
        ))
        .as_millis() as u64;
    alloc_track::snap!("Phase 2: Calibrate");
    // Print tolerance summary
    if let Some(summary) = calibration.ridge_width_summary() {
        println!(
            "  RT tolerance (ridge): avg {:.0}s, min {:.0}s, max {:.0}s ({} cols, {:.0}% in-ridge)",
            summary.weighted_half_width,
            summary.min_half_width,
            summary.max_half_width,
            summary.n_columns,
            summary.in_ridge_ratio * 100.0,
        );
    }
    // The mobility tolerance is only meaningful for a searchable TIMS 1/K0 run;
    // for mzML/FAIMS it is fit from sentinel mobilities and never used (the
    // query gate unrestricts mobility), so report it as disabled rather than
    // print a misleading number.
    if pipeline.index.mobility_kind().is_scoreable() {
        println!(
            "  m/z: ({:.1}, {:.1}) ppm   mobility: ({:.1}, {:.1}) %",
            calibration.mz_tolerance().0,
            calibration.mz_tolerance().1,
            calibration.mobility_tolerance().0,
            calibration.mobility_tolerance().1,
        );
    } else {
        println!(
            "  m/z: ({:.1}, {:.1}) ppm   mobility: disabled (no searchable axis)",
            calibration.mz_tolerance().0,
            calibration.mz_tolerance().1,
        );
    }

    // Save the calibration for the viewer to load.
    if !calibration.is_fallback() {
        let (rt_lo_ms, rt_hi_ms) = pipeline.index.ms1_cycle_mapping().range_milis();
        let rt_lo = rt_lo_ms as f64 / 1000.0;
        let rt_hi = rt_hi_ms as f64 / 1000.0;
        let cal_json_path = std::path::Path::new(&out_path.uri).join("calibration.json");
        if let Err(e) = calibration.save_json([rt_lo, rt_hi], phase1_lib.len(), &cal_json_path) {
            tracing::warn!("Failed to save calibration: {}", e);
        } else {
            info!("Saved calibration to {:?}", cal_json_path);
        }
    }

    calib_dash_hook::show_final(&mut calib_dash_state, &calibration);

    // === PHASE 3: Narrow scoring with calibrated tolerances ===
    info!("Phase 3: Scoring with calibrated extraction...");
    timscentroid::indexing::reset_for_each_peak_funnel();
    let step = TimedStep::begin("Phase 3: Score");
    let mut phase3_timings = ScoreTimings::default();
    let (results, phase3_skips) = phase3_score(
        speclib,
        pipeline,
        &calibration,
        chunk_size,
        &mut phase3_timings,
    );
    step.finish_with(format_args!("{} peptides", results.len()));
    alloc_track::snap!("Phase 3: Score");
    timscentroid::indexing::dump_for_each_peak_funnel("Phase 3");
    timscentroid::indexing::reset_for_each_peak_funnel();

    let total_scored = results.len();

    // === PHASE 4: Target-decoy competition ===
    let step = TimedStep::begin("Phase 4: Compete");
    let mut competed = target_decoy_compete(results);
    competed.sort_unstable_by(|x, y| {
        y.scoring
            .primary
            .main_score
            .partial_cmp(&x.scoring.primary.main_score)
            .expect("NaN main_score should have been filtered during Phase 3 scoring")
    });
    let total_after_competition = competed.len();
    let phase4_ms = step
        .finish_with(format_args!("{} candidates", total_after_competition))
        .as_millis() as u64;
    alloc_track::snap!("Phase 4: Compete");

    // === PHASE 5: Rescore ===
    let step = TimedStep::begin("Phase 5: Rescore");
    // The CLI flag overrides the configured model; MLP is the default.
    info!("Phase 5 rescore model: {rescore_model:?}");
    let (data, feature_stats) = rescore_with(rescore_model, competed)?;
    let phase5_ms = step.finish().as_millis() as u64;
    alloc_track::snap!("Phase 5: Rescore");

    // Collect q-value threshold counts -- full report to log, key result to stdout
    let qval_report = report_qvalues_at_thresholds(&data, &[0.01, 0.05, 0.1, 0.5, 1.0]);
    let mut targets_at_1pct_qval = 0usize;
    let mut targets_at_5pct_qval = 0usize;
    let mut targets_at_10pct_qval = 0usize;
    for &ThresholdCounts { q, targets, decoys } in &qval_report {
        info!(
            "q-value threshold {:.2}: {} targets, {} decoys ({} total)",
            q,
            targets,
            decoys,
            targets + decoys
        );
        if (q - 0.01).abs() < 1e-6 {
            targets_at_1pct_qval = targets;
        } else if (q - 0.05).abs() < 1e-6 {
            targets_at_5pct_qval = targets;
        } else if (q - 0.10).abs() < 1e-6 {
            targets_at_10pct_qval = targets;
        }
    }

    // Built BEFORE Phase 6, because the writer below consumes `data`. So the
    // precompute -- including its ~1 GB feature matrix -- runs while nothing is
    // on disk yet; only the blocking TUI is after the write.
    #[cfg(feature = "dashboard")]
    let dashboard = crate::dashboard::build(&data, &feature_stats, &qval_report);

    // === PHASE 6: Write Parquet output ===
    let step = TimedStep::begin("Phase 6: Write output");
    // ARTIFACT-LIST (per-sample): keep in sync with validate_inputs in main.rs.
    let out_path_pq = std::path::Path::new(&out_path.uri).join("results.parquet");
    let mut pq_writer = timsseek::scoring::parquet_writer::ResultParquetWriter::new(
        &out_path_pq,
        20_000,
        speclib.parsable_sequences(),
        &speclib.geom,
    )
    .map_err(|e| TimsSeekError::Io {
        path: out_path_pq.clone().into(),
        source: e,
    })?;
    for res in data.into_iter() {
        if res.qvalue <= max_qvalue {
            pq_writer.add(res).map_err(|e| TimsSeekError::Io {
                path: out_path_pq.clone().into(),
                source: e,
            })?;
        }
    }
    pq_writer.close().map_err(|e| TimsSeekError::Io {
        path: out_path_pq.clone().into(),
        source: e,
    })?;
    let phase6_ms = step.finish().as_millis() as u64;
    alloc_track::snap!("Phase 6: Write output");
    info!("Wrote final results to {:?}", out_path_pq);

    if !no_feature_stats && let Err(e) = write_feature_stats_sidecar(&feature_stats, &out_path_pq) {
        // Non-fatal: log and continue.
        tracing::warn!("Failed to write feature_stats sidecar: {}", e);
    }

    // After Phase 6, so a dashboard left open overnight -- or killed -- still has
    // its results written.
    // No `snap!` around this: it blocks until the user quits, so any measurement
    // taken across it reports how long they looked at the screen.
    #[cfg(feature = "dashboard")]
    if let Some(dash) = dashboard
        && let Err(e) = rescore_dash::run(dash)
    {
        tracing::warn!("rescore dashboard exited with an error: {e}");
    }

    // Key result to stdout. The final output URI is printed by main.rs
    // per-file footer -- out_path_pq here is the local working path (which
    // is a tempdir for remote destinations), not the eventual location.
    println!();
    println!("{} targets at 1% FDR", targets_at_1pct_qval);

    Ok(PipelineReport {
        load_index_ms: 0, // set by caller after return
        phase1_prescore_ms: phase1_ms,
        phase1_detail: phase1_timings,
        phase2_calibration_ms: phase2_ms,
        phase3_extraction_thread_ms: phase3_timings.extraction.as_millis() as u64,
        phase3_scoring_thread_ms: phase3_timings.scoring.as_millis() as u64,
        phase3_spectral_query_thread_ms: phase3_timings.spectral_query.as_millis() as u64,
        phase3_assembly_thread_ms: phase3_timings.assembly.as_millis() as u64,
        phase4_competition_ms: phase4_ms,
        phase5_rescore_ms: phase5_ms,
        phase6_output_ms: phase6_ms,
        total_scored,
        total_after_competition,
        targets_at_1pct_qval,
        targets_at_5pct_qval,
        targets_at_10pct_qval,
        phase3_skips,
    })
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
fn phase1_prescore<I: ScorerQueriable>(
    speclib: &ReferenceLibrary,
    pipeline: &Scorer<I>,
    chunk_size: usize,
    config: &CalibrationConfig,
) -> (
    Vec<CalibrantCandidate>,
    timsseek::scoring::PrescoreTimings,
    calib_dash_hook::Dash,
) {
    let total = speclib.len();
    let n_chunks = total.div_ceil(chunk_size);
    let pb = make_progress_bar(n_chunks as u64, "Phase 1");

    let mut dash = calib_dash_hook::attach(n_chunks, config);

    let mut global_heap = CalibrantHeap::new(config.n_calibrants);
    let mut timings = timsseek::scoring::PrescoreTimings::default();

    // Chunk the flat index space `0..len` (no materialized slice on the lazy
    // arm); each flat index is also the global speclib index.
    for (batch_idx, batch) in speclib.chunks(chunk_size).enumerate().progress_with(pb) {
        let chunk_heap = pipeline.prescore_batch(speclib, &batch, config, &mut timings);
        global_heap = global_heap.merge(chunk_heap);

        calib_dash_hook::on_batch(&mut dash, batch_idx, global_heap.iter());
    }
    calib_dash_hook::finish(&mut dash);

    (global_heap.into_vec(), timings, dash)
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

fn calibrate_from_phase1<I: ScorerQueriable>(
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
        .map(|c| {
            let calib_item = phase1_lib.item_at(c.speclib_index);
            match main_lookup {
                Some(lookup) => {
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
                            if shared >= MIN_SHARED_FRAGMENTS {
                                Some((shared, (main_rt - calib_rt).abs(), *main_rt))
                            } else {
                                None
                            }
                        })
                        .min_by(|a, b| {
                            b.0.cmp(&a.0).then(
                                a.1.partial_cmp(&b.1)
                                    .expect("NaN RT residual in calibrant matching"),
                            )
                        })
                        .map(|(_, _, rt)| rt as f64)
                }
                None => Some(calib_item.rt_seconds() as f64),
            }
        })
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

    let windows = errors.derive_windows(&derivation);
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

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
fn phase3_score<I: ScorerQueriable>(
    speclib: &ReferenceLibrary,
    pipeline: &Scorer<I>,
    calibration: &CalibrationResult,
    chunk_size: usize,
    timings: &mut ScoreTimings,
) -> (Vec<ScoredCandidate>, SkipCounts) {
    let total_peptides = speclib.len();
    let n_chunks = total_peptides.div_ceil(chunk_size);
    let pb = make_progress_bar(n_chunks as u64, "Phase 3");

    let mut results = Vec::new();
    let mut skips = SkipCounts::default();

    // Chunk the flat index space `0..len` and hand each range to the batch
    // scorer, which drives the lazy flyweight by index.
    for batch in speclib.chunks(chunk_size).progress_with(pb) {
        let (batch_results, batch_timings, batch_skips) =
            pipeline.score_calibrated_batch(speclib, &batch, calibration);
        *timings += batch_timings;
        skips += batch_skips;
        results.extend(batch_results);
    }

    let skipped = total_peptides - results.len();
    if skipped > total_peptides / 20 {
        warn!(
            "{}/{} peptides produced no Phase 3 result (>{:.0}%): {}. \
             If this is unexpected, check calibration quality.",
            skipped,
            total_peptides,
            (skipped as f64 / total_peptides as f64) * 100.0,
            skips,
        );
    }

    (results, skips)
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
/// Sort by `(sequence, main_score desc, target-first, precursor_mz)` then
/// collapse exact `(sequence, charge, precursor_mz)` duplicates, keeping the
/// first. The trailing `precursor_mz` tiebreak makes the sort a TOTAL order:
/// within a shared sequence the target and its ±decoys have distinct precursor
/// m/z, so the survivor is deterministic regardless of the input vec's order
/// (an unstable sort would otherwise leave a `(seq, score, is_target)` tie in
/// arbitrary relative order).
fn dedup_by_sequence(results: &mut Vec<ScoredCandidate>) {
    results.sort_unstable_by(|x, y| {
        let seq_ord = x
            .scoring
            .identity
            .peptide
            .as_str()
            .cmp(y.scoring.identity.peptide.as_str());
        // Then sort descending by main_score
        // NOTE: same sequences should always have the same score EXCEPT when we apply a mass shift
        // to some of them to make a "decoy"
        let score_ord = y
            .scoring
            .primary
            .main_score
            .partial_cmp(&x.scoring.primary.main_score)
            .expect("NaN main_score should have been filtered during Phase 3 scoring");
        let ord = seq_ord.then(score_ord);

        if ord == std::cmp::Ordering::Equal {
            // Move to the first position the target
            match (x.scoring.identity.is_target, y.scoring.identity.is_target) {
                (true, false) => std::cmp::Ordering::Less,
                (false, true) => std::cmp::Ordering::Greater,
                // Total order: within a shared sequence the target and its
                // ±decoys have distinct precursor m/z (mono / mono±shift/z),
                // so this breaks the tie deterministically.
                _ => x
                    .scoring
                    .identity
                    .precursor_mz
                    .total_cmp(&y.scoring.identity.precursor_mz),
            }
        } else {
            ord
        }
    });
    results.dedup_by(|x, y| {
        (x.scoring.identity.peptide.as_str() == y.scoring.identity.peptide.as_str())
            && (x.scoring.identity.precursor_charge == y.scoring.identity.precursor_charge)
            && (x.scoring.identity.precursor_mz == y.scoring.identity.precursor_mz)
    });
}

fn target_decoy_compete(mut results: Vec<ScoredCandidate>) -> Vec<CompetedCandidate> {
    // TODO: re-implement so we dont drop results but instead just flag them as rejected (maybe
    // a slice where we push rejected results to the end and keep the trailing slice as the "active")

    fn glimpse_result_head(results: &[ScoredCandidate]) -> Vec<String> {
        results[..10.min(results.len())]
            .iter()
            .map(|x| {
                format!(
                    "{} {} {} {}",
                    x.scoring.identity.peptide.as_str(),
                    x.scoring.identity.precursor_charge,
                    x.scoring.identity.precursor_mz,
                    x.scoring.primary.main_score
                )
            })
            .collect::<Vec<_>>()
    }
    // Deduplicate by sequence, keeping the best scoring target
    // This is meant to remove instances where reversing a target creates another target.
    debug!(
        "First 10 result before deduplication for seq+charge+mz: {:#?}",
        glimpse_result_head(&results)
    );
    dedup_by_sequence(&mut results);
    debug!(
        "First 10 result after deduplication for seq+charge+mz: {:#?}",
        glimpse_result_head(&results)
    );

    // Sort competing results adjacent, best first within each competition.
    results.sort_unstable_by(|x, y| {
        x.scoring
            .identity
            .competition_key()
            .cmp(&y.scoring.identity.competition_key())
            .then_with(|| {
                x.scoring
                    .primary
                    .main_score
                    .partial_cmp(&y.scoring.primary.main_score)
                    .expect("NaN main_score should have been filtered during Phase 3 scoring")
                    .reverse()
            })
    });
    info!(
        "Number of results before t/d competition: {}",
        results.len()
    );

    // Each run is one competition, best first. Only the winner survives, and
    // its separation feature is its margin over the runner-up -- NaN when it
    // ran alone, since there is nothing to separate from and NaN is the model's
    // missing marker.
    let group_lens: Vec<usize> = results
        .chunk_by(|a, b| a.scoring.identity.competes_with(&b.scoring.identity))
        .map(<[_]>::len)
        .collect();

    info!(
        "Number of results after t/d competition: {}",
        group_lens.len()
    );

    let mut members = results.into_iter();
    group_lens
        .into_iter()
        .map(|len| {
            let winner = members.next().expect("chunk_by yields non-empty runs");
            let ln1p_winner = winner.scoring.primary.main_score.ln_1p();
            // Guard on the run length before taking: a lone winner must not
            // consume the next competition's winner as its runner-up.
            let (diff, ratio) = if len >= 2 {
                let runner_up = members.next().expect("a run of >=2 has a second member");
                let ln1p_runner_up = runner_up.scoring.primary.main_score.ln_1p();
                (ln1p_winner - ln1p_runner_up, ln1p_runner_up / ln1p_winner)
            } else {
                (f32::NAN, f32::NAN)
            };
            // Everyone below the runner-up loses and carries no feature.
            members.by_ref().take(len.saturating_sub(2)).for_each(drop);
            winner.into_competed(diff, ratio)
        })
        .collect()
}

pub fn run_pipeline(
    speclib: &ReferenceLibrary,
    calib_lib: Option<&ReferenceLibrary>,
    pipeline: &Scorer<IndexedTimstofPeaks>,
    chunk_size: usize,
    output: &OutputConfig,
    max_qvalue: f32,
    load_index_ms: u64,
    calib_config: &CalibrationConfig,
    no_feature_stats: bool,
    rescore_model: RescoreModel,
) -> std::result::Result<PipelineReport, TimsSeekError> {
    // ARTIFACT-LIST (per-sample): keep in sync with validate_inputs in main.rs.
    let performance_report_path = std::path::Path::new(&output.uri).join("performance_report.json");

    let mut timings = execute_pipeline(
        speclib,
        calib_lib,
        pipeline,
        chunk_size,
        output,
        max_qvalue,
        calib_config,
        no_feature_stats,
        rescore_model,
    )?;
    timings.load_index_ms = load_index_ms;
    // Write per-file report
    let perf_report =
        serde_json::to_string_pretty(&timings).map_err(|e| TimsSeekError::ParseError {
            msg: format!("Error serializing performance report to JSON: {}", e),
        })?;
    std::fs::write(&performance_report_path, perf_report).map_err(|e| TimsSeekError::Io {
        path: performance_report_path.into(),
        source: e,
    })?;
    Ok(timings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use timsquery::models::test_handles;
    use timsseek::models::DecoyMarking;
    use timsseek::models::sequence::Peptide;
    use timsseek::scoring::results::ScoringFields;

    /// `group` names the competition group: candidates sharing one compete.
    /// Production reads it off the arena as an opaque code; a fixture only needs
    /// the codes to differ where the test wants different groups.
    fn candidate(seq: &str, mz: f64, is_target: bool, group: u32) -> ScoredCandidate {
        let decoy = if is_target {
            DecoyMarking::Target
        } else {
            DecoyMarking::MassShiftedDecoy
        };
        let peptide = Peptide {
            raw: Arc::from(seq),
            decoy,
            sequence_features: false,
        };
        let mut scoring = ScoringFields::sample(peptide);
        scoring.identity.precursor_mz = mz;
        scoring.identity.is_target = is_target;
        scoring.identity.group = test_handles::group(group);
        // All fixtures tie on score so the (seq, score, is_target) tiebreak arm
        // is what orders them.
        scoring.primary.main_score = 5.0;
        ScoredCandidate { scoring }
    }

    fn survivors(mut results: Vec<ScoredCandidate>) -> Vec<(String, u64, bool)> {
        dedup_by_sequence(&mut results);
        results
            .iter()
            .map(|c| {
                (
                    c.scoring.identity.peptide.as_str().to_string(),
                    c.scoring.identity.precursor_mz.to_bits(),
                    c.scoring.identity.is_target,
                )
            })
            .collect()
    }

    #[test]
    fn dedup_total_order() {
        // One shared sequence: a target and its ±decoys, all tied on main_score
        // and charge, distinct precursor m/z. Two input orderings differing ONLY
        // in the tie order must produce identical (ordered) survivors.
        let order_a = vec![
            candidate("PEPTIDEK", 501.0, true, 0),
            candidate("PEPTIDEK", 500.0, false, 0),
            candidate("PEPTIDEK", 502.0, false, 0),
        ];
        let order_b = vec![
            candidate("PEPTIDEK", 502.0, false, 0),
            candidate("PEPTIDEK", 500.0, false, 0),
            candidate("PEPTIDEK", 501.0, true, 0),
        ];

        let sa = survivors(order_a);
        let sb = survivors(order_b);
        assert_eq!(sa, sb, "dedup survivors must be order-independent");
        // Target first (target-first arm), then decoys ascending by m/z.
        assert_eq!(
            sa,
            vec![
                ("PEPTIDEK".to_string(), 501.0f64.to_bits(), true),
                ("PEPTIDEK".to_string(), 500.0f64.to_bits(), false),
                ("PEPTIDEK".to_string(), 502.0f64.to_bits(), false),
            ]
        );
    }

    #[test]
    fn competition_features_match_their_ln1p_names() {
        let mut best = candidate("BEST", 501.0, true, 7);
        best.scoring.primary.main_score = 8.0;
        let mut runner_up = candidate("RUNNER", 502.0, false, 7);
        runner_up.scoring.primary.main_score = 3.0;

        let competed = target_decoy_compete(vec![runner_up, best]);
        assert_eq!(competed.len(), 1);
        let winner = &competed[0];
        let best_ln1p = 8.0f32.ln_1p();
        let runner_up_ln1p = 3.0f32.ln_1p();

        assert!((winner.delta_group_ln1p_diff - (best_ln1p - runner_up_ln1p)).abs() < 1e-6);
        assert!((winner.delta_group_ln1p_ratio - runner_up_ln1p / best_ln1p).abs() < 1e-6);
    }

    /// `MassShift` makes every group three members, so a
    /// two-member test never exercised the case that mattered: the margin has
    /// to come from the runner-up, not from the worst member.
    #[test]
    fn a_three_member_group_separates_the_winner_from_the_runner_up() {
        let mut best = candidate("BEST", 501.0, true, 7);
        best.scoring.primary.main_score = 8.0;
        let mut middle = candidate("MIDDLE", 502.0, false, 7);
        middle.scoring.primary.main_score = 3.0;
        let mut worst = candidate("WORST", 503.0, false, 7);
        worst.scoring.primary.main_score = 1.0;

        let competed = target_decoy_compete(vec![worst, middle, best]);
        assert_eq!(competed.len(), 1, "one group, one survivor");

        let (b, m) = (8.0f32.ln_1p(), 3.0f32.ln_1p());
        assert!((competed[0].delta_group_ln1p_diff - (b - m)).abs() < 1e-6);
        assert!((competed[0].delta_group_ln1p_ratio - m / b).abs() < 1e-6);
    }

    /// Nothing to separate from, so the feature does not apply. NaN is the
    /// model's missing marker; 0 would claim a tie with a rival that never
    /// existed.
    #[test]
    fn a_lone_group_member_reports_no_separation() {
        let mut only = candidate("ALONE", 501.0, true, 7);
        only.scoring.primary.main_score = 8.0;

        let competed = target_decoy_compete(vec![only]);
        assert_eq!(competed.len(), 1);
        assert!(competed[0].delta_group_ln1p_diff.is_nan());
        assert!(competed[0].delta_group_ln1p_ratio.is_nan());
    }

    /// Groups are found by walking runs of the sort order, so a result that
    /// belongs to a different competition must not be folded into its
    /// neighbour's.
    #[test]
    fn separate_groups_do_not_bleed_into_each_other() {
        let mut a = candidate("A", 501.0, true, 1);
        a.scoring.primary.main_score = 8.0;
        let mut b = candidate("B", 502.0, true, 2);
        b.scoring.primary.main_score = 3.0;

        let competed = target_decoy_compete(vec![a, b]);
        assert_eq!(competed.len(), 2, "two groups, two survivors");
        assert!(competed.iter().all(|c| c.delta_group_ln1p_diff.is_nan()));
    }
}
