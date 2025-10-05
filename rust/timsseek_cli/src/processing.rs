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

pub fn main_loop<I: GenerallyQueriable<IonAnnot>>(
    // query_iterator: impl ExactSizeIterator<Item = QueryItemToScore>,
    // # I would like this to be streaming
    query_iterator: Speclib,
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
