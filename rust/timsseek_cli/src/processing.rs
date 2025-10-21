use super::config::OutputConfig;
use indicatif::{
    ProgressIterator,
    ProgressStyle,
};
use std::path::PathBuf;
use std::time::Instant;
use timsquery::{
    GenerallyQueriable,
    IndexedTimstofPeaks,
};
use timsseek::IonAnnot;
use timsseek::data_sources::speclib::Speclib;
use timsseek::errors::TimsSeekError;
use timsseek::ml::qvalues::report_qvalues_at_thresholds;
use timsseek::ml::rescore;
use timsseek::rt_calibration::recalibrate_speclib;
use timsseek::scoring::scorer::{
    ScoreTimings,
    Scorer,
};
use timsseek::scoring::search_results::{
    IonSearchResults,
    ResultParquetWriter,
};
use tracing::{
    debug,
    info,
};

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn main_loop<I: GenerallyQueriable<IonAnnot>>(
    // query_iterator: impl ExactSizeIterator<Item = QueryItemToScore>,
    // # I would like this to be streaming
    mut query_iterator: Speclib,
    scorer: &Scorer<I>,
    chunk_size: usize,
    out_path: &OutputConfig,
) -> std::result::Result<ScoreTimings, TimsSeekError> {
    let total = query_iterator.len();
    let mut chunk_num = 0;
    let mut nqueried = 0;
    let mut nwritten = 0;
    let start = Instant::now();

    let out_path_pq = out_path.directory.join("results.parquet");
    let mut pq_writer = ResultParquetWriter::new(out_path_pq.clone(), 20_000).map_err(|e| {
        tracing::error!(
            "Error creating parquet writer for path {:?}: {}",
            out_path_pq,
            e
        );
        TimsSeekError::Io {
            path: out_path_pq.into(),
            source: e,
        }
    })?;
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
            // Parallelism happens here within the score_iter function
            let (mut out, timings): (Vec<IonSearchResults>, ScoreTimings) = scorer.score_iter(chunk);
            all_timings += timings;
            nwritten += out.len();
            out.sort_unstable_by(|x, y| x.main_score.partial_cmp(&y.main_score).unwrap());
            debug!("Worst score in chunk: {:#?}", out[0]);
            if let Some(last) = out.last() {
                debug!("Best Score in chunk: {:#?}", last);
            }
            results.extend(out.iter().cloned());

            for x in out.into_iter() {
                pq_writer.add(x);
            }
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

    pq_writer.close();
    println!("Processed {} queries, wrote {} results", nqueried, nwritten);
    println!(
        "Finished processing {} chunks in {:?}",
        chunk_num,
        start.elapsed()
    );

    // This is a really dirty deduplication step ... since some targets can be decoys as well
    // we just sort by sequence and keep the best scoring one
    results.sort_unstable_by(|x, y| {
        let seq_ord = x.sequence.cmp(&y.sequence);
        if seq_ord == std::cmp::Ordering::Equal {
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
    results.dedup_by(|x, y| x.sequence == y.sequence);

    // Sort in descending order of score
    results.sort_unstable_by(|x, y| y.main_score.partial_cmp(&x.main_score).unwrap());
    assert!(results.first().unwrap().main_score >= results.last().unwrap().main_score);

    match recalibrate_speclib(
        &mut query_iterator,
        &results[..(100_000.min(results.len()))],
    ) {
        Ok(_) => info!("Recalibrated speclib retention times based on search results"),
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
    let out_path_pq = out_path.directory.join("results_rescored.parquet");
    let mut pq_writer = ResultParquetWriter::new(out_path_pq.clone(), 20_000).map_err(|e| {
        tracing::error!(
            "Error creating parquet writer for path {:?}: {}",
            out_path_pq,
            e
        );
        TimsSeekError::Io {
            path: out_path_pq.into(),
            source: e,
        }
    })?;
    for res in data.into_iter() {
        pq_writer.add(res);
    }
    pq_writer.close();

    Ok(all_timings)
}

pub fn process_speclib(
    path: PathBuf,
    scorer: &Scorer<IndexedTimstofPeaks>,
    chunk_size: usize,
    output: &OutputConfig,
) -> std::result::Result<(), TimsSeekError> {
    // TODO: I should probably "inline" this function with the main loop
    info!("Building database from speclib file {:?}", path);
    let st = std::time::Instant::now();
    let performance_report_path = output.directory.join("performance_report.json");
    let speclib = Speclib::from_file(&path)?;
    let elap_time = st.elapsed();
    info!(
        "Loading speclib of length {} took: {:?} for {}",
        speclib.len(),
        elap_time,
        path.display()
    );
    let timings = main_loop(speclib, scorer, chunk_size, output)?;
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
