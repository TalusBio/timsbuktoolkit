use super::config::{
    AnalysisConfig,
    DigestionConfig,
    OutputConfig,
};
use indicatif::{
    ProgressIterator,
    ProgressStyle,
};
use rayon::prelude::*;
use timsseek::fragment_mass::IonAnnot;
use timsseek::scoring::scorer::Scorer;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{
    Duration,
    Instant,
};
use timsquery::models::aggregators::EGCAggregator;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsquery::{
    GenerallyQueriable, QueriableData, Tolerance
};
use timsseek::data_sources::speclib::Speclib;
use timsseek::digest::digestion::{
    DigestionEnd,
    DigestionParameters,
    DigestionPattern,
};
use timsseek::errors::TimsSeekError;
use timsseek::fragment_mass::elution_group_converter::SequenceToElutionGroupConverter;
use timsseek::models::{
    DigestSlice,
    DigestedSequenceIterator,
    NamedQueryChunk,
    deduplicate_digests,
};
use timsseek::protein::fasta::ProteinSequenceCollection;
use timsseek::scoring::calculate_scores::{
    IntensityArrays,
    LocalizedPreScore,
    LongitudinalMainScoreElements,
    PreScore,
};
use timsseek::scoring::full_results::FullQueryResult;
use timsseek::scoring::search_results::{
    IonSearchResults,
    SearchResultBuilder,
    write_results_to_parquet,
};
use tracing::info;


pub fn main_loop<'a, I: GenerallyQueriable<IonAnnot>>(
    chunked_query_iterator: impl ExactSizeIterator<Item = NamedQueryChunk>,
    scorer: &'a Scorer<I>,
    tolerance: &'a Tolerance,
    out_path: &OutputConfig,
) -> std::result::Result<(), TimsSeekError> {
    let mut chunk_num = 0;
    let start = Instant::now();

    let style = ProgressStyle::with_template(
        "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
    )
    .unwrap();
    chunked_query_iterator
        .progress_with_style(style)
        .for_each(|chunk| {
            // Here I am branching whether I want the 'full output', which is
            // which writes to disk both the search results and the extractions
            // to disk. If that is not requested, the (hypothetically) more
            // efficient version is used.
            let out: (Vec<IonSearchResults>, RuntimeMetrics) = scorer.score_iter(chunk);
            let out_path_pq = &out_path
                .directory
                .join(format!("chunk_{}.parquet", chunk_num));
            write_results_to_parquet(&out.0, out_path_pq).unwrap();
            chunk_num += 1;
        });
    println!(
        "Finished processing {} chunks in {:?}",
        chunk_num,
        start.elapsed()
    );
    Ok(())
}

pub fn process_speclib(
    path: PathBuf,
    index: &QuadSplittedTransposedIndex,
    analysis: &AnalysisConfig,
    output: &OutputConfig,
) -> std::result::Result<(), TimsSeekError> {
    tracing::info!("Building database from speclib file {:?}", path);
    let st = std::time::Instant::now();
    let speclib = Speclib::from_ndjson_file(&path)?;
    let elap_time = st.elapsed();
    println!(
        "Loading speclib of length {} took: {:?} for {}",
        speclib.len(),
        elap_time,
        path.display()
    );
    let speclib_iter = speclib.as_iterator(analysis.chunk_size);

    main_loop(speclib_iter, index, &analysis.tolerance, output)?;
    Ok(())
}
