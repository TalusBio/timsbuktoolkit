use super::config::OutputConfig;
use indicatif::{
    ProgressIterator,
    ProgressStyle,
};
use std::path::Path;
use std::time::Instant;
use timsquery::IndexedTimstofPeaks;
use timsquery::MzMobilityStatsCollector;
use timsquery::SpectralCollector;
use timsquery::Tolerance;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
};
use timsseek::data_sources::speclib::Speclib;
use timsseek::errors::TimsSeekError;
use timsseek::ml::qvalues::report_qvalues_at_thresholds;
use timsseek::ml::rescore;
use timsseek::rt_calibration::{
    CalibRtError,
    CalibrationResult,
    Point,
    calibrate_with_ranges,
};
use timsseek::scoring::{
    CalibrantCandidate,
    CalibrantHeap,
    CalibrationConfig,
    CompetedCandidate,
    PipelineTimings,
    ScoredCandidate,
    ScoreTimings,
};
use timsseek::scoring::pipeline::Scorer;
use timsseek::{
    DecoyStrategy,
    IonAnnot,
    ScorerQueriable,
};
use tracing::{
    debug,
    info,
    warn,
};

/// Check that two speclibs are on a compatible RT scale.
/// Warns loudly if the RT ranges don't overlap, which would produce a useless calibration.
fn check_rt_scale_compatibility(main_lib: &Speclib, calib_lib: &Speclib) {
    fn rt_range(lib: &Speclib) -> (f32, f32) {
        let mut min_rt = f32::INFINITY;
        let mut max_rt = f32::NEG_INFINITY;
        for item in lib.as_slice() {
            let rt = item.query.rt_seconds();
            min_rt = min_rt.min(rt);
            max_rt = max_rt.max(rt);
        }
        (min_rt, max_rt)
    }

    let (main_min, main_max) = rt_range(main_lib);
    let (calib_min, calib_max) = rt_range(calib_lib);

    info!(
        "RT ranges — main speclib: [{:.1}, {:.1}]s, calib lib: [{:.1}, {:.1}]s",
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
            "!! RT SCALE WARNING: only {:.0}% overlap between main speclib !!", overlap_pct
        );
        warn!("!! and calibration library. Calibration may be unreliable.  !!");
        warn!("!! Ensure both libraries use the same iRT scale.            !!");
        warn!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    } else if overlap_pct < 80.0 {
        warn!(
            "RT overlap between main speclib and calib lib is {:.0}% — may affect calibration at the extremes",
            overlap_pct
        );
    }
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn main_loop<I: ScorerQueriable>(
    speclib: Speclib,
    calib_lib: Option<Speclib>,
    pipeline: &Scorer<I>,
    chunk_size: usize,
    _out_path: &OutputConfig,
) -> std::result::Result<PipelineTimings, TimsSeekError> {
    let calib_config = CalibrationConfig::default();

    // === PHASE 1: Broad prescore -> collect top calibrants ===
    // Use calibration library if provided, otherwise fall back to main speclib
    let phase1_lib = calib_lib.as_ref().unwrap_or(&speclib);
    if let Some(ref clib) = calib_lib {
        info!(
            "Phase 1: Broad prescore using calibration library ({} entries)...",
            clib.len()
        );
        check_rt_scale_compatibility(&speclib, clib);
    } else {
        info!("Phase 1: Broad prescore (unrestricted RT)...");
    }
    let phase1_start = Instant::now();
    let calibrants = phase1_prescore(phase1_lib, pipeline, chunk_size, &calib_config);
    let phase1_ms = phase1_start.elapsed().as_millis() as u64;
    info!(
        "Phase 1 complete: {} calibrant candidates in {}ms",
        calibrants.len(),
        phase1_ms
    );

    // === PHASE 2: Calibration (fit RT + measure errors + derive tolerances) ===
    // Build lookup from main speclib when using a separate calib lib.
    // Maps (quantized_precursor_mz, charge) -> Vec<(rt, sorted_fragment_mzs)>.
    // Matching requires same precursor (0.01 Da) + charge + at least 5 shared fragment masses.
    let main_lookup: Option<
        std::collections::HashMap<(i64, u8), Vec<(f32, Vec<i64>)>>,
    > = if calib_lib.is_some() {
        let mut map: std::collections::HashMap<(i64, u8), Vec<(f32, Vec<i64>)>> =
            std::collections::HashMap::new();
        for item in speclib.as_slice() {
            let mz_key = (item.query.mono_precursor_mz() * 100.0).round() as i64;
            let charge = item.query.precursor_charge();
            let mut frag_mzs: Vec<i64> = item
                .query
                .iter_fragments()
                .map(|(_, mz)| (mz * 100.0).round() as i64)
                .collect();
            frag_mzs.sort_unstable();
            map.entry((mz_key, charge))
                .or_default()
                .push((item.query.rt_seconds(), frag_mzs));
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
    let phase2_start = Instant::now();
    let calibration = match calibrate_from_phase1(
        calibrants,
        phase1_lib,
        main_lookup.as_ref(),
        pipeline,
        &calib_config,
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
    let phase2_ms = phase2_start.elapsed().as_millis() as u64;

    // === PHASE 3: Narrow scoring with calibrated tolerances ===
    info!("Phase 3: Scoring with calibrated extraction...");
    let phase3_start = Instant::now();
    let mut phase3_timings = ScoreTimings::default();
    let results = phase3_score(
        &speclib,
        pipeline,
        &calibration,
        chunk_size,
        &mut phase3_timings,
    );
    info!(
        "Phase 3 complete: {} scored peptides in {:?}",
        results.len(),
        phase3_start.elapsed()
    );

    // === Post-processing ===
    let mut competed = target_decoy_compete(results);
    competed.sort_unstable_by(|x, y| {
        y.scoring.main_score.partial_cmp(&x.scoring.main_score).unwrap()
    });

    let data = rescore(competed);
    for val in report_qvalues_at_thresholds(&data, &[0.01, 0.05, 0.1, 0.5, 1.0]) {
        let (thresh, n_below_thresh, n_targets, n_decoys) = val;
        println!(
            "Found {} targets and {} decoys at q-value threshold {:.2} ({} total)",
            n_targets, n_decoys, thresh, n_below_thresh
        );
    }

    // TODO: Task 7 — manual Parquet writer for FinalResult
    let _ = &data;

    Ok(PipelineTimings {
        phase1_prescore_ms: phase1_ms,
        phase2_calibration_ms: phase2_ms,
        phase3_prescore_ms: phase3_timings.prescore.as_millis() as u64,
        phase3_localize_ms: phase3_timings.localize.as_millis() as u64,
        phase3_secondary_query_ms: phase3_timings.secondary_query.as_millis() as u64,
        phase3_finalization_ms: phase3_timings.finalization.as_millis() as u64,
    })
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
fn phase1_prescore<I: ScorerQueriable>(
    speclib: &Speclib,
    pipeline: &Scorer<I>,
    chunk_size: usize,
    config: &CalibrationConfig,
) -> Vec<CalibrantCandidate> {
    let style = ProgressStyle::with_template(
        "{spinner:.green} Phase 1 [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
    )
    .unwrap();

    let mut global_heap = CalibrantHeap::new(config.n_calibrants);
    let mut offset = 0usize;

    for chunk in speclib.as_slice().chunks(chunk_size).progress_with_style(style) {
        let chunk_heap = pipeline.prescore_batch(chunk, offset, config);
        global_heap = global_heap.merge(chunk_heap);
        offset += chunk.len();
    }

    global_heap.into_vec()
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
    phase1_lib: &Speclib,
    main_lookup: Option<&std::collections::HashMap<(i64, u8), Vec<(f32, Vec<i64>)>>>,
    pipeline: &Scorer<I>,
    config: &CalibrationConfig,
) -> Result<CalibrationResult, CalibRtError> {
    // === Step A: Fit iRT -> RT curve ===
    // When a separate calib lib is used, we need the main speclib's iRT as x.
    // The calibration curve must map main_speclib_irt -> observed_rt.
    let points: Vec<Point> = candidates
        .iter()
        .filter_map(|c| {
            let calib_item = &phase1_lib.as_slice()[c.speclib_index];

            let irt_for_curve = match main_lookup {
                Some(lookup) => {
                    let mz_key =
                        (calib_item.query.mono_precursor_mz() * 100.0).round() as i64;
                    let charge = calib_item.query.precursor_charge();
                    let bucket = lookup.get(&(mz_key, charge))?;

                    // Build sorted fragment m/z list for the calib entry
                    let mut calib_frags: Vec<i64> = calib_item
                        .query
                        .iter_fragments()
                        .map(|(_, mz)| (mz * 100.0).round() as i64)
                        .collect();
                    calib_frags.sort_unstable();

                    // Find best match: most shared fragments, break ties by closest RT
                    let calib_rt = calib_item.query.rt_seconds();
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
                        // Best = most shared fragments, then closest RT
                        .min_by(|a, b| b.0.cmp(&a.0).then(a.1.partial_cmp(&b.1).unwrap()))
                        .map(|(_, _, rt)| rt)?
                }
                None => calib_item.query.rt_seconds(),
            };

            Some(Point {
                x: irt_for_curve as f64,
                y: c.apex_rt_seconds as f64,
                weight: 1.0,
            })
        })
        .collect();

    if main_lookup.is_some() {
        info!(
            "Calibration: {} of {} calibrants matched in main speclib (>={} shared fragments)",
            points.len(),
            candidates.len(),
            MIN_SHARED_FRAGMENTS,
        );
    }

    let (min_x, max_x, min_y, max_y) = points.iter().fold(
        (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ),
        |(mnx, mxx, mny, mxy), p| (mnx.min(p.x), mxx.max(p.x), mny.min(p.y), mxy.max(p.y)),
    );

    let cal_curve =
        calibrate_with_ranges(&points, (min_x, max_x), (min_y, max_y), config.grid_size)?;

    // === Step B: Measure m/z and mobility errors at calibrant apexes ===
    let query_tolerance = Tolerance {
        ms: MzTolerance::Ppm((10.0, 10.0)),
        rt: RtTolerance::Minutes((
            config.calibration_query_rt_window_minutes,
            config.calibration_query_rt_window_minutes,
        )),
        mobility: MobilityTolerance::Pct((5.0, 5.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1)),
    };

    let mut mz_errors_ppm: Vec<f32> = Vec::with_capacity(candidates.len());
    let mut mobility_errors_pct: Vec<f32> = Vec::with_capacity(candidates.len());

    for candidate in &candidates {
        let item = &phase1_lib.as_slice()[candidate.speclib_index];
        let query_at_apex = item
            .query
            .clone()
            .with_rt_seconds(candidate.apex_rt_seconds);
        let mut agg: SpectralCollector<IonAnnot, MzMobilityStatsCollector> =
            SpectralCollector::new(query_at_apex);
        pipeline.index.add_query(&mut agg, &query_tolerance);

        for ((_key, expected_mz), stats) in agg.iter_precursors() {
            if let (Ok(obs_mz), Ok(obs_mob)) = (stats.mean_mz(), stats.mean_mobility()) {
                let mz_err = (obs_mz - expected_mz) / expected_mz * 1e6;
                mz_errors_ppm.push(mz_err as f32);

                let expected_mob = item.query.mobility_ook0() as f64;
                let mob_err = (obs_mob - expected_mob) / expected_mob * 100.0;
                mobility_errors_pct.push(mob_err as f32);
                break;
            }
        }
    }

    // === Step C: Derive tolerances from error distributions ===
    let rt_tolerance_minutes = {
        let mut abs_residuals: Vec<f64> = points
            .iter()
            .map(|p| {
                let predicted = cal_curve.predict(p.x).unwrap_or(p.y);
                (p.y - predicted).abs()
            })
            .collect();
        abs_residuals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mad_seconds = abs_residuals
            .get(abs_residuals.len() / 2)
            .copied()
            .unwrap_or(0.0);
        info!(
            "RT residuals: MAD={:.1}s, n={}",
            mad_seconds, abs_residuals.len()
        );
        (config.rt_sigma_factor * mad_seconds as f32 / 60.0).max(config.min_rt_tolerance_minutes)
    };

    let mz_tolerance_ppm = {
        let (l, r) = asymmetric_tolerance(&mz_errors_ppm, config.mz_sigma, 0.1);
        (l as f64, r as f64)
    };

    let mobility_tolerance_pct =
        asymmetric_tolerance(&mobility_errors_pct, config.mobility_sigma, 0.1);

    info!(
        "Calibration: RT tol={:.2} min, m/z tol=({:.1}, {:.1}) ppm, mob tol=({:.1}, {:.1}) %",
        rt_tolerance_minutes,
        mz_tolerance_ppm.0,
        mz_tolerance_ppm.1,
        mobility_tolerance_pct.0,
        mobility_tolerance_pct.1,
    );

    Ok(CalibrationResult::new(
        cal_curve,
        rt_tolerance_minutes,
        mz_tolerance_ppm,
        mobility_tolerance_pct,
    ))
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
fn phase3_score<I: ScorerQueriable>(
    speclib: &Speclib,
    pipeline: &Scorer<I>,
    calibration: &CalibrationResult,
    chunk_size: usize,
    timings: &mut ScoreTimings,
) -> Vec<ScoredCandidate> {
    let style = ProgressStyle::with_template(
        "{spinner:.green} Phase 3 [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
    )
    .unwrap();

    let total_peptides = speclib.as_slice().len();
    let mut results = Vec::new();

    for chunk in speclib
        .as_slice()
        .chunks(chunk_size)
        .progress_with_style(style)
    {
        let (batch_results, batch_timings) =
            pipeline.score_calibrated_batch(chunk, calibration);
        *timings += batch_timings;
        results.extend(batch_results);
    }

    let skipped = total_peptides - results.len();
    if skipped > total_peptides / 20 {
        warn!(
            "{}/{} peptides produced no Phase 3 result (>{:.0}%). \
             If this is unexpected, check calibration quality.",
            skipped,
            total_peptides,
            (skipped as f64 / total_peptides as f64) * 100.0
        );
    }

    results
}

/// Derive asymmetric tolerance from error distribution.
#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
fn asymmetric_tolerance(errors: &[f32], n_sigma: f32, min_val: f32) -> (f32, f32) {
    if errors.is_empty() {
        return (min_val, min_val);
    }
    let mean = errors.iter().sum::<f32>() / errors.len() as f32;
    let variance = errors.iter().map(|e| (e - mean).powi(2)).sum::<f32>() / errors.len() as f32;
    let std = variance.sqrt();
    let left = (-(mean - n_sigma * std)).max(min_val);
    let right = (mean + n_sigma * std).max(min_val);
    (left, right)
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
fn target_decoy_compete(mut results: Vec<ScoredCandidate>) -> Vec<CompetedCandidate> {
    // TODO: re-implement so we dont drop results but instead just flag them as rejected (maybe
    // a slice where we push rejected results to the end and keep the trailing slice as the "active")

    fn glimpse_result_head(results: &[ScoredCandidate]) -> Vec<String> {
        results[..10.min(results.len())]
            .iter()
            .map(|x| {
                format!(
                    "{} {} {} {}",
                    x.scoring.sequence, x.scoring.precursor_charge, x.scoring.precursor_mz, x.scoring.main_score
                )
            })
            .collect::<Vec<_>>()
    }
    // Deduplicate by sequence, keeping the best scoring target
    // This is meant to remove instances where reversing a target creates another target.
    results.sort_unstable_by(|x, y| {
        let seq_ord = x.scoring.sequence.cmp(&y.scoring.sequence);
        // Then sort descending by main_score
        // NOTE: same sequences should always have the same score EXCEPT when we apply a mass shift
        // to some of them to make a "decoy"
        let score_ord = y.scoring.main_score.partial_cmp(&x.scoring.main_score).unwrap();
        let ord = seq_ord.then(score_ord);

        if ord == std::cmp::Ordering::Equal {
            // Move to the first position the target
            match (x.scoring.is_target, y.scoring.is_target) {
                (true, false) => std::cmp::Ordering::Less,
                (false, true) => std::cmp::Ordering::Greater,
                _ => std::cmp::Ordering::Equal,
            }
        } else {
            seq_ord
        }
    });
    // As debug lets print the first and last results after deduplication
    debug!(
        "First 10 result before deduplication for seq+charge+mz: {:#?}",
        glimpse_result_head(&results)
    );
    results.dedup_by(|x, y| {
        (x.scoring.sequence == y.scoring.sequence)
            && (x.scoring.precursor_charge == y.scoring.precursor_charge)
            && (x.scoring.precursor_mz == y.scoring.precursor_mz)
    });
    debug!(
        "First 10 result after deduplication for seq+charge+mz: {:#?}",
        glimpse_result_head(&results)
    );

    // Compete target-decoy pairs at precursor level
    results.sort_unstable_by(|x, y| {
        x.scoring.decoy_group_id
            .cmp(&y.scoring.decoy_group_id)
            .then_with(|| x.scoring.precursor_charge.cmp(&y.scoring.precursor_charge))
            .then_with(|| x.scoring.main_score.partial_cmp(&y.scoring.main_score).unwrap().reverse())
    });
    info!(
        "Number of results before t/d competition: {}",
        results.len()
    );

    // Calculate delta scores between consecutive target/decoy pairs
    // Results are sorted by (decoy_group_id, precursor_charge, score desc)
    // We store (group_id, charge, index, main_score) and the computed deltas per index.
    let mut delta_map: Vec<(f32, f32)> = vec![(f32::NAN, f32::NAN); results.len()];
    let mut previous: Option<(u32, u8, usize, f32)> = None;

    for i in 0..results.len() {
        let current = &results[i];
        let current_key = (current.scoring.decoy_group_id, current.scoring.precursor_charge);

        if let Some((prev_group_id, prev_charge, prev_index, prev_score)) = previous {
            let prev_key = (prev_group_id, prev_charge);

            if current_key == prev_key {
                // This is the second item in a target/decoy pair
                let delta_score = current.scoring.main_score - prev_score;
                let delta_ratio = current.scoring.main_score / prev_score;

                delta_map[prev_index] = (-delta_score, delta_ratio);

                // Skip updating previous - we only compare first two items per group
                continue;
            }
        }

        // Start of a new group or first item overall
        previous = Some((
            current.scoring.decoy_group_id,
            current.scoring.precursor_charge,
            i,
            current.scoring.main_score,
        ));
    }

    // Dedup by (decoy_group_id, charge) — keep the first (best scoring)
    // We need indices to grab the right deltas, so collect the deduped indices first.
    let mut kept_indices: Vec<usize> = Vec::with_capacity(results.len());
    {
        let mut last_key: Option<(u32, u8)> = None;
        for i in 0..results.len() {
            let key = (results[i].scoring.decoy_group_id, results[i].scoring.precursor_charge);
            if last_key == Some(key) {
                continue; // duplicate in same group
            }
            last_key = Some(key);
            kept_indices.push(i);
        }
    }

    info!("Number of results after t/d competition: {}", kept_indices.len());

    // Build CompetedCandidate vec from the kept indices.
    // We need to pull elements out of `results` by index, but they are non-Copy.
    // Convert the whole Vec into an indexed form we can drain.
    let mut results_opt: Vec<Option<ScoredCandidate>> = results.into_iter().map(Some).collect();
    let competed: Vec<CompetedCandidate> = kept_indices
        .into_iter()
        .map(|i| {
            let (dg, dgr) = delta_map[i];
            results_opt[i]
                .take()
                .expect("index should be unique")
                .into_competed(dg, dgr)
        })
        .collect();

    competed
}

pub fn process_speclib(
    path: &Path,
    calib_lib_path: Option<&Path>,
    pipeline: &Scorer<IndexedTimstofPeaks>,
    chunk_size: usize,
    output: &OutputConfig,
    decoy_strategy: DecoyStrategy,
) -> std::result::Result<(), TimsSeekError> {
    info!("Building database from speclib file {:?}", path);
    info!("Decoy generation strategy: {}", decoy_strategy);

    let st = std::time::Instant::now();
    let performance_report_path = output.directory.join("performance_report.json");
    let speclib = Speclib::from_file(path, decoy_strategy)?;
    let elap_time = st.elapsed();
    info!(
        "Loading speclib of length {} took: {:?} for {}",
        speclib.len(),
        elap_time,
        path.display()
    );

    let calib_lib = match calib_lib_path {
        Some(p) => {
            info!("Loading calibration library from {:?}", p);
            let st = std::time::Instant::now();
            let lib = Speclib::from_file(p, decoy_strategy)?;
            info!(
                "Loaded calibration library of length {} in {:?}",
                lib.len(),
                st.elapsed()
            );
            Some(lib)
        }
        None => None,
    };

    let timings = main_loop(speclib, calib_lib, pipeline, chunk_size, output)?;
    let perf_report =
        serde_json::to_string_pretty(&timings).map_err(|e| TimsSeekError::ParseError {
            msg: format!("Error serializing performance report to JSON: {}", e),
        })?;
    std::fs::write(&performance_report_path, perf_report).map_err(|e| TimsSeekError::Io {
        path: performance_report_path.into(),
        source: e,
    })?;
    Ok(())
}
