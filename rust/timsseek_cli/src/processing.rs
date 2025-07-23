use super::config::OutputConfig;
use indicatif::{
    ProgressIterator,
    ProgressStyle,
};
use std::path::PathBuf;
use std::time::Instant;
use timsquery::GenerallyQueriable;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsseek::IonAnnot;
use timsseek::data_sources::speclib::Speclib;
use timsseek::errors::TimsSeekError;
use timsseek::scoring::scorer::Scorer;
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
) -> std::result::Result<(), TimsSeekError> {
    let mut chunk_num = 0;
    let mut nqueried = 0;
    let mut nwritten = 0;
    let start = Instant::now();

    let out_path_pq = out_path.directory.join("results.parquet");
    let mut pq_writer = ResultParquetWriter::new(out_path_pq, 20_000).unwrap();
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
            let mut out: Vec<IonSearchResults> = scorer.score_iter(chunk);
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
        });

    pq_writer.close();
    println!("Processed {} queries, wrote {} results", nqueried, nwritten);
    println!(
        "Finished processing {} chunks in {:?}",
        chunk_num,
        start.elapsed()
    );
    Ok(())
}

pub fn process_speclib(
    path: PathBuf,
    scorer: &Scorer<QuadSplittedTransposedIndex>,
    chunk_size: usize,
    output: &OutputConfig,
) -> std::result::Result<(), TimsSeekError> {
    // TODO: I should probably "inline" this function with the main loop
    info!("Building database from speclib file {:?}", path);
    let st = std::time::Instant::now();
    let speclib = Speclib::from_file(&path)?;
    let elap_time = st.elapsed();
    info!(
        "Loading speclib of length {} took: {:?} for {}",
        speclib.len(),
        elap_time,
        path.display()
    );

    main_loop(speclib, scorer, chunk_size, output)?;
    Ok(())
}
