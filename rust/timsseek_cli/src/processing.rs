use super::config::OutputConfig;
use crate::artifacts::{
    FEATURE_IMPORTANCE_TSV,
    FEATURE_STATS_TSV,
    PERFORMANCE_REPORT,
    RESULTS_PARQUET,
};
use indicatif::ProgressIterator;
use timsquery::IndexedTimstofPeaks;
use timsseek::ScorerQueriable;
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
use timsseek::rt_calibration::CalibrationResult;
use timsseek::scoring::pipeline::Scorer;
use timsseek::scoring::timings::{
    TimedStep,
    make_progress_bar,
};
use timsseek::scoring::{
    CalibrationConfig,
    CompetedCandidate,
    PipelineReport,
    ScoreTimings,
    ScoredCandidate,
    SkipCounts,
};
use tracing::{
    debug,
    info,
    warn,
};

mod calibration;

use calibration::{
    build_precursor_fragment_lookup,
    calib_dash_hook,
    calibrate_from_phase1,
    check_rt_scale_compatibility,
    phase1_prescore,
};

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
        &parquet_path.with_file_name(FEATURE_STATS_TSV),
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
        &parquet_path.with_file_name(FEATURE_IMPORTANCE_TSV),
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

/// Decoy-free libraries use full-RT spectrum scoring and a schema without FDR fields.
fn execute_raw_pipeline<I: ScorerQueriable>(
    library: &ReferenceLibrary,
    pipeline: &Scorer<I>,
    options: &PipelineOptions<'_>,
) -> Result<PipelineReport, TimsSeekError> {
    info!(
        "No decoys: raw spectrum scoring over full acquisition RT; calibration, competition, rescoring and q-value filtering disabled"
    );
    let path = std::path::Path::new(&options.output.uri).join(RESULTS_PARQUET);
    let io_error = |source| TimsSeekError::Io {
        source,
        path: Some(path.clone()),
    };
    let mut writer =
        timsseek::scoring::parquet_writer::ResultParquetWriter::raw(&path, 20_000, library)
            .map_err(io_error)?;
    let mut report = PipelineReport {
        raw_scores: true,
        ..Default::default()
    };
    let mut timings = ScoreTimings::default();
    for batch in library.chunks(options.chunk_size) {
        let (rows, timing, skips) = pipeline.score_raw_batch(library, &batch);
        timings += timing;
        report.phase3_skips += skips;
        report.total_scored += rows.len();
        for row in rows {
            writer.add_raw(row).map_err(io_error)?;
        }
    }
    writer.close().map_err(io_error)?;
    report.phase3_extraction_thread_ms = timings.extraction.as_millis() as u64;
    report.phase3_scoring_thread_ms = timings.scoring.as_millis() as u64;
    report.phase3_spectral_query_thread_ms = timings.spectral_query.as_millis() as u64;
    report.phase3_assembly_thread_ms = timings.assembly.as_millis() as u64;
    println!(
        "{} raw scored candidates; no FDR estimates",
        report.total_scored
    );
    Ok(report)
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub struct PipelineOptions<'a> {
    pub chunk_size: usize,
    pub output: &'a OutputConfig,
    pub max_qvalue: f32,
    pub calibration: &'a CalibrationConfig,
    pub no_feature_stats: bool,
    pub rescore_model: RescoreModel,
}

pub fn execute_pipeline<I: ScorerQueriable>(
    speclib: &ReferenceLibrary,
    calib_lib: Option<&ReferenceLibrary>,
    pipeline: &Scorer<I>,
    options: &PipelineOptions<'_>,
) -> std::result::Result<PipelineReport, TimsSeekError> {
    if !matches!(
        speclib.geometry().capabilities().decoys,
        timsquery::models::capabilities::DecoyStrategy::MassShift { .. }
    ) && speclib.geometry().n_stored_decoys() == 0
    {
        return execute_raw_pipeline(speclib, pipeline, options);
    }
    let PipelineOptions {
        chunk_size,
        output: out_path,
        max_qvalue,
        calibration: calib_config,
        no_feature_stats,
        rescore_model,
    } = *options;
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
    let main_lookup = calib_lib.map(|_| build_precursor_fragment_lookup(speclib));
    let mapped_rt_range = calib_lib.and_then(|_| speclib.rt_range());
    // Reset the `for_each_peak` filter-funnel counters so the per-phase
    // dump reflects just this phase's work. Both calls no-op when the
    // `query-instr` feature is disabled.
    timscentroid::indexing::reset_for_each_peak_funnel();
    let step = TimedStep::begin("Phase 1: Prescore");
    let (calibrants, phase1_timings, mut calib_dash_state) = phase1_prescore(
        phase1_lib,
        pipeline,
        chunk_size,
        calib_config,
        main_lookup.as_ref(),
        mapped_rt_range,
    );
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
    info!("Phase 2: Calibration...");
    let step = TimedStep::begin("Phase 2: Calibrate");
    let calibration = {
        match calibrate_from_phase1(
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
        }
    };
    let calibration = calibration.with_rt_axis(speclib.geometry().rt_axis());
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
    println!(
        "  m/z: {:?}   mobility: {}   RT fit: {}",
        calibration.mz_tolerance(),
        if pipeline.index.mobility_kind().is_scoreable() {
            format!("{:?}", calibration.mobility_tolerance())
        } else {
            "disabled (no searchable axis)".to_owned()
        },
        calibration.has_rt_calibration()
    );

    // Save the calibration for the viewer to load.
    {
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
    let mut competed = target_decoy_compete(results, speclib);
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
    let (data, feature_stats) = rescore_with(rescore_model, competed, speclib)?;
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
    let dashboard = crate::dashboard::build(speclib, &data, &feature_stats, &qval_report);

    // === PHASE 6: Write Parquet output ===
    let step = TimedStep::begin("Phase 6: Write output");
    let out_path_pq = std::path::Path::new(&out_path.uri).join(RESULTS_PARQUET);
    let mut pq_writer =
        timsseek::scoring::parquet_writer::ResultParquetWriter::new(&out_path_pq, 20_000, speclib)
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
        raw_scores: false,
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
            "{}/{} peptides ({:.1}%) produced no Phase 3 result: {}. \
             Check library/acquisition m/z coverage for range skips; check RT calibration only \
             for retention-time skips.",
            skipped,
            total_peptides,
            (skipped as f64 / total_peptides as f64) * 100.0,
            skips.percentage_summary(total_peptides),
        );
    }

    (results, skips)
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
/// Group complete peptide structure + charge + m/z, then retain the best score,
/// preferring a target on ties. Incomparable chemistry stays row-specific.
fn dedup_by_sequence(results: &mut Vec<ScoredCandidate>, library: &ReferenceLibrary) {
    let chemistry_cmp = |a: &ScoredCandidate, b: &ScoredCandidate| {
        let a = &a.scoring.identity;
        let b = &b.scoring.identity;
        let pa = library
            .geometry()
            .analyte(a.row)
            .peptide
            .known()
            .filter(|p| p.modifications.known().is_some());
        let pb = library
            .geometry()
            .analyte(b.row)
            .peptide
            .known()
            .filter(|p| p.modifications.known().is_some());
        match (pa, pb) {
            (Some(pa), Some(pb)) => pa.residues.cmp(pb.residues).then_with(|| {
                pa.modifications
                    .known()
                    .unwrap()
                    .iter()
                    .map(|(site, m)| (site, m.annotation()))
                    .cmp(
                        pb.modifications
                            .known()
                            .unwrap()
                            .iter()
                            .map(|(site, m)| (site, m.annotation())),
                    )
            }),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => a.row.cmp(&b.row),
        }
    };
    // Sort by the entire duplicate key before ranking survivors. Otherwise a
    // different charge/mz can separate equal keys and escape deduplication.
    let key_cmp = |a: &ScoredCandidate, b: &ScoredCandidate| {
        chemistry_cmp(a, b)
            .then(
                a.scoring
                    .identity
                    .precursor_charge
                    .cmp(&b.scoring.identity.precursor_charge),
            )
            .then(
                a.scoring
                    .identity
                    .precursor_mz
                    .total_cmp(&b.scoring.identity.precursor_mz),
            )
    };
    results.sort_unstable_by(|a, b| {
        key_cmp(a, b)
            .then(
                b.scoring
                    .primary
                    .main_score
                    .total_cmp(&a.scoring.primary.main_score),
            )
            .then(
                b.scoring
                    .identity
                    .is_target
                    .cmp(&a.scoring.identity.is_target),
            )
            .then(
                a.scoring
                    .identity
                    .source_id
                    .cmp(&b.scoring.identity.source_id),
            )
            .then(a.scoring.identity.row.cmp(&b.scoring.identity.row))
    });
    results.dedup_by(|a, b| key_cmp(a, b).is_eq());
}

/// Targets discarded because another target won their competition, which is the
/// search-side symptom of a library whose groups hold more than one target per
/// charge.
///
/// A target beaten by a *decoy* is ordinary target/decoy competition and is not
/// counted -- that is the mechanism working. Every target past the first in one
/// `(group, charge)` is counted whichever member won, since one of the two is
/// gone either way.
///
/// Reads the vec the competition sort left behind, so competitors are adjacent
/// and this is one pass with no allocation. It has to run after
/// `dedup_by_sequence`: the rows that removes are one precursor listed twice,
/// which is a duplicate rather than a mis-declared group, and counting them here
/// would report a defect the library does not have.
fn targets_beaten_by_targets(sorted: &[ScoredCandidate]) -> usize {
    sorted
        .chunk_by(|a, b| a.scoring.identity.competes_with(&b.scoring.identity))
        .map(|run| {
            run.iter()
                .filter(|c| c.scoring.identity.is_target)
                .count()
                .saturating_sub(1)
        })
        .sum()
}

fn target_decoy_compete(
    mut results: Vec<ScoredCandidate>,
    library: &ReferenceLibrary,
) -> Vec<CompetedCandidate> {
    // TODO: re-implement so we dont drop results but instead just flag them as rejected (maybe
    // a slice where we push rejected results to the end and keep the trailing slice as the "active")

    fn glimpse_result_head(results: &[ScoredCandidate]) -> Vec<String> {
        results[..10.min(results.len())]
            .iter()
            .map(|x| {
                format!(
                    "{} {} {} {}",
                    x.scoring.identity.source_id,
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
    dedup_by_sequence(&mut results, library);
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
    let outcompeted_targets = targets_beaten_by_targets(&results);
    if outcompeted_targets > 0 {
        warn!(
            "{outcompeted_targets} targets were discarded because another target won their \
             competition, which is a real identification thrown away to keep another: the \
             library declares competition groups holding more than one target per charge."
        );
    }

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
    load_index_ms: u64,
    options: &PipelineOptions<'_>,
) -> std::result::Result<PipelineReport, TimsSeekError> {
    let performance_report_path =
        std::path::Path::new(&options.output.uri).join(PERFORMANCE_REPORT);

    let mut timings = execute_pipeline(speclib, calib_lib, pipeline, options)?;
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
    fn fixture_library(results: &mut [ScoredCandidate]) -> ReferenceLibrary {
        use timsquery::models::{
            Row,
            TargetColumnsBuilder,
        };
        use timsquery::serde::TargetTable;
        let mut builder = TargetColumnsBuilder::with_capabilities(
            timsquery::models::TargetCapabilities::default_diann(),
        );
        let frags = [(timsquery::ion::IonAnnot::try_from("y1").unwrap(), 300.0)];
        for (i, c) in results.iter_mut().enumerate() {
            c.scoring.identity.row = test_handles::row(i as u32);
            let analyte = timsquery::chemistry::analyte::Analyte::from_sequence(
                &c.scoring.identity.source_id.to_string(),
            );
            builder.push_row(Row {
                analyte: analyte.as_input(),
                charge: 2,
                precursor_mz: 500.0,
                frags: &frags,
                ..Default::default()
            });
        }
        ReferenceLibrary::try_from(TargetTable::Mzpaf {
            geom: builder.seal(Default::default()).unwrap(),
            frag_intens: Some(vec![1.0; results.len()]),
        })
        .unwrap()
    }
    fn compete(mut results: Vec<ScoredCandidate>) -> Vec<CompetedCandidate> {
        let library = fixture_library(&mut results);
        target_decoy_compete(results, &library)
    }

    use super::*;
    use timsquery::models::test_handles;
    use timsseek::scoring::results::ScoringFields;

    /// `group` names the competition group: candidates sharing one compete.
    /// Production reads it off the arena as an opaque code; a fixture only needs
    /// the codes to differ where the test wants different groups.
    fn candidate(seq: &str, mz: f64, is_target: bool, group: u32) -> ScoredCandidate {
        let mut scoring = ScoringFields::sample_default();
        scoring.identity.precursor_mz = mz;
        scoring.identity.source_id = seq.into();
        scoring.identity.is_target = is_target;
        scoring.identity.group = test_handles::group(group);
        // All fixtures tie on score so the (seq, score, is_target) tiebreak arm
        // is what orders them.
        scoring.primary.main_score = 5.0;
        ScoredCandidate { scoring }
    }

    fn survivors(mut results: Vec<ScoredCandidate>) -> Vec<(String, u64, bool)> {
        let library = fixture_library(&mut results);
        dedup_by_sequence(&mut results, &library);
        results
            .iter()
            .map(|c| {
                (
                    c.scoring.identity.source_id.to_string(),
                    c.scoring.identity.precursor_mz.to_bits(),
                    c.scoring.identity.is_target,
                )
            })
            .collect()
    }

    #[test]
    fn structural_dedup_groups_charge_and_mz_before_ranking() {
        let mut a = candidate("AC(UniMod:4)K", 500.0, false, 0);
        let mut b = candidate("AC[UNIMOD:4]K", 500.0, true, 0);
        let mut different_charge = candidate("AC[UNIMOD:4]K", 500.0, true, 1);
        different_charge.scoring.identity.precursor_charge = 3;
        a.scoring.primary.main_score = 5.0;
        b.scoring.primary.main_score = 5.0;
        different_charge.scoring.primary.main_score = 6.0;
        let mut rows = vec![a, different_charge, b];
        let library = fixture_library(&mut rows);
        dedup_by_sequence(&mut rows, &library);
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().all(|r| r.scoring.identity.is_target));
        assert_eq!(rows[0].scoring.identity.precursor_charge, 2);
        assert_eq!(rows[1].scoring.identity.precursor_charge, 3);
    }

    #[test]
    fn missing_and_unresolved_chemistry_do_not_form_one_duplicate_key() {
        let mut rows = vec![
            candidate("", 500.0, true, 0),
            candidate("", 500.0, false, 1),
            candidate("unknown!", 500.0, true, 2),
            candidate("unknown!", 500.0, false, 3),
            candidate("PEPTIDE", 500.0, true, 4),
        ];
        let library = fixture_library(&mut rows);
        dedup_by_sequence(&mut rows, &library);
        assert_eq!(rows.len(), 5);
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
        // Duplicate keys are adjacent; distinct m/z values sort ascending.
        assert_eq!(
            sa,
            vec![
                ("PEPTIDEK".to_string(), 500.0f64.to_bits(), false),
                ("PEPTIDEK".to_string(), 501.0f64.to_bits(), true),
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

        let competed = compete(vec![runner_up, best]);
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

        let competed = compete(vec![worst, middle, best]);
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

        let competed = compete(vec![only]);
        assert_eq!(competed.len(), 1);
        assert!(competed[0].delta_group_ln1p_diff.is_nan());
        assert!(competed[0].delta_group_ln1p_ratio.is_nan());
    }

    /// A target beaten by a decoy is what competition is for, so it counts as
    /// nothing lost.
    #[test]
    fn a_target_losing_to_a_decoy_is_not_counted_as_an_outcompeted_target() {
        let mut target = candidate("PEPTIDEK", 501.0, true, 7);
        target.scoring.primary.main_score = 3.0;
        let mut decoy = candidate("KEDITPEP", 502.0, false, 7);
        decoy.scoring.primary.main_score = 8.0;

        assert_eq!(targets_beaten_by_targets(&[decoy, target]), 0);
    }

    /// Two targets competing at one charge means one of them is discarded to
    /// keep the other, which is the search-side symptom of a mis-declared
    /// group and the number the log has to name.
    #[test]
    fn a_target_losing_to_another_target_is_counted() {
        let mut winner = candidate("PEPTIDEK", 501.0, true, 7);
        winner.scoring.primary.main_score = 8.0;
        let mut loser = candidate("SAMPLERK", 502.0, true, 7);
        loser.scoring.primary.main_score = 3.0;

        assert_eq!(targets_beaten_by_targets(&[winner, loser]), 1);
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

        let competed = compete(vec![a, b]);
        assert_eq!(competed.len(), 2, "two groups, two survivors");
        assert!(competed.iter().all(|c| c.delta_group_ln1p_diff.is_nan()));
    }
}
