use super::config::OutputConfig;
use indicatif::{
    ProgressIterator,
    ProgressStyle,
};
use std::path::Path;
use std::time::Instant;
use timsquery::IndexedTimstofPeaks;
use timsseek::data_sources::speclib::Speclib;
use timsseek::errors::TimsSeekError;
use timsseek::ml::qvalues::report_qvalues_at_thresholds;
use timsseek::ml::rescore;
use timsseek::rt_calibration::{
    recalibrate_results,
    recalibrate_speclib,
};
use timsseek::scoring::ScoreTimings;
use timsseek::scoring::pipeline::ScoringPipeline;
use timsseek::scoring::search_results::{
    IonSearchResults,
    ResultParquetWriter,
};
use timsseek::{
    DecoyStrategy,
    ScorerQueriable,
};
use tracing::{
    debug,
    info,
};

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn main_loop<I: ScorerQueriable>(
    // query_iterator: impl ExactSizeIterator<Item = QueryItemToScore>,
    // # I would like this to be streaming
    mut query_iterator: Speclib,
    pipeline: &ScoringPipeline<I>,
    chunk_size: usize,
    out_path: &OutputConfig,
) -> std::result::Result<ScoreTimings, TimsSeekError> {
    let total = query_iterator.len();
    let mut chunk_num = 0;
    let mut nqueried = 0;
    let mut nwritten = 0;
    let start = Instant::now();

    let mut all_timings = ScoreTimings::default();
    let style = ProgressStyle::with_template(
        "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
    )
    .unwrap();
    let mut results: Vec<IonSearchResults> = Vec::new();
    query_iterator
        .as_slice()
        .chunks(chunk_size)
        .progress_with_style(style)
        .for_each(|chunk| {
            nqueried += chunk.len();
            // Parallelism happens here within the process_batch function
            let (mut out, timings): (Vec<IonSearchResults>, ScoreTimings) = pipeline.process_batch(chunk);
            all_timings += timings;
            nwritten += out.len();
            out.sort_unstable_by(|x, y| x.main_score.partial_cmp(&y.main_score).unwrap());
            debug!("Worst score in chunk: {:#?}", out[0]);
            if let Some(last) = out.last() {
                debug!("Best Score in chunk: {:#?}", last);
            }
            results.extend(out.iter().cloned());

            chunk_num += 1;
            let pct_done = (nqueried as f64 / total as f64) * 100.0;
            let estimated_total = start.elapsed().as_secs_f64() * (100.0 / pct_done);
            let eta = (estimated_total - start.elapsed().as_secs_f64()).max(0.0);
            let eta_duration = std::time::Duration::from_secs_f64(eta);
            info!(
                "Processed chunk {}, total queries {}, total written {}, elapsed {:?}, {:.2}% done, ETA {:?}",
                chunk_num, nqueried, nwritten, start.elapsed(), pct_done, eta_duration
            );
        });

    info!("Processed {} queries, found {} results", nqueried, nwritten);
    info!(
        "Finished processing {} chunks in {:?}",
        chunk_num,
        start.elapsed()
    );

    let mut results = target_decoy_compete(results);

    // Sort in descending order of score
    results.sort_unstable_by(|x, y| y.main_score.partial_cmp(&x.main_score).unwrap());
    assert!(results.first().unwrap().main_score >= results.last().unwrap().main_score);

    match recalibrate_speclib(&mut query_iterator, &results) {
        Ok(calib) => {
            info!("Recalibrated speclib retention times based on search results");
            recalibrate_results(&calib, results.as_mut_slice());
        }
        Err(e) => {
            tracing::error!("Error recalibrating speclib retention times: {:?}", e);
        }
    };

    let data = rescore(results);
    for val in report_qvalues_at_thresholds(&data, &[0.01, 0.05, 0.1, 0.5, 1.0]) {
        let (thresh, n_below_thresh, n_targets, n_decoys) = val;
        println!(
            "Found {} targets and {} decoys at q-value threshold {:.2} ({} total)",
            n_targets, n_decoys, thresh, n_below_thresh
        );
    }
    let out_path_pq = out_path.directory.join("results.parquet");
    let mut pq_writer = ResultParquetWriter::new(out_path_pq.clone(), 20_000).map_err(|e| {
        tracing::error!(
            "Error creating parquet writer for path {:?}: {}",
            out_path_pq.clone(),
            e
        );
        TimsSeekError::Io {
            path: out_path_pq.clone().into(),
            source: e,
        }
    })?;
    for res in data.into_iter() {
        pq_writer.add(res);
    }
    pq_writer.close();
    info!("Wrote final results to {:?}", out_path_pq);

    Ok(all_timings)
}

fn target_decoy_compete(mut results: Vec<IonSearchResults>) -> Vec<IonSearchResults> {
    // TODO: re-implement so we dont drop results but instead just flag them as rejected (maybe
    // a slice where we push rejected results to the end and keep the trailing slice as the "active")

    // I KNOW this is an ugly place for a function...
    fn glimpse_result_head(results: &[IonSearchResults]) -> Vec<String> {
        results[..10.min(results.len())]
            .iter()
            .map(|x| {
                format!(
                    "{} {} {} {}",
                    x.sequence, x.precursor_charge, x.precursor_mz, x.main_score
                )
            })
            .collect::<Vec<_>>()
    }
    // Deduplicate by sequence, keeping the best scoring target
    // This is meant to remove instances where reversing a target creates another target.
    results.sort_unstable_by(|x, y| {
        let seq_ord = x.sequence.cmp(&y.sequence);
        // Then sort descending by main_score
        // NOTE: same sequences should always have the same score EXCEPT when we apply a mass shift
        // to some of them to make a "decoy"
        let score_ord = y.main_score.partial_cmp(&x.main_score).unwrap();
        let ord = seq_ord.then(score_ord);

        if ord == std::cmp::Ordering::Equal {
            // Move to the first position the target
            match (x.is_target, y.is_target) {
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
        (x.sequence == y.sequence)
            && (x.precursor_charge == y.precursor_charge)
            && (x.precursor_mz == y.precursor_mz)
    });
    debug!(
        "First 10 result after deduplication for seq+charge+mz: {:#?}",
        glimpse_result_head(&results)
    );

    // Compete target-decoy pairs at precursor level
    results.sort_unstable_by(|x, y| {
        x.decoy_group_id
            .cmp(&y.decoy_group_id)
            .then_with(|| x.precursor_charge.cmp(&y.precursor_charge))
            .then_with(|| x.main_score.partial_cmp(&y.main_score).unwrap().reverse())
    });
    info!(
        "Number of results before t/d competition: {}",
        results.len()
    );
    debug!(
        "First 10 result before target-decoy in decoy group id: {:#?}",
        glimpse_result_head(&results)
    );

    // Calculate delta scores between consecutive target/decoy pairs
    // Results are sorted by (decoy_group_id, precursor_charge, score desc)
    let mut previous: Option<(u32, u8, usize, f32)> = None;

    for i in 0..results.len() {
        let current = &results[i];
        let current_key = (current.decoy_group_id, current.precursor_charge);

        if let Some((prev_group_id, prev_charge, prev_index, prev_score)) = previous {
            let prev_key = (prev_group_id, prev_charge);

            if current_key == prev_key {
                // This is the second item in a target/decoy pair
                let delta_score = current.main_score - prev_score;
                let delta_ratio = current.main_score / prev_score;

                results[prev_index].delta_group = -delta_score;
                results[prev_index].delta_group_ratio = delta_ratio;

                // Skip updating previous - we only compare first two items per group
                continue;
            }
        }

        // Start of a new group or first item overall
        previous = Some((
            current.decoy_group_id,
            current.precursor_charge,
            i,
            current.main_score,
        ));
    }

    results.dedup_by(|x, y| {
        (x.decoy_group_id == y.decoy_group_id) & (x.precursor_charge == y.precursor_charge)
    });
    debug!(
        "First 10 result after deduplication for decoy_group_id+charge: {:#?}",
        glimpse_result_head(&results)
    );
    info!("Number of results after t/d competition: {}", results.len());
    debug!(
        "First 10 result after target-decoy competition: {:#?}",
        glimpse_result_head(&results)
    );
    results
}

pub fn process_speclib(
    path: &Path,
    pipeline: &ScoringPipeline<IndexedTimstofPeaks>,
    chunk_size: usize,
    output: &OutputConfig,
    decoy_strategy: DecoyStrategy,
) -> std::result::Result<(), TimsSeekError> {
    // TODO: I should probably "inline" this function with the main loop
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
    let timings = main_loop(speclib, pipeline, chunk_size, output)?;
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
