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
use std::path::{
    Path,
    PathBuf,
};
use std::sync::Arc;
use std::time::Instant;
use timsquery::models::aggregators::MultiCMGStatsFactory;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsquery::queriable_tims_data::queriable_tims_data::query_multi_group;
use timsquery::traits::tolerance::DefaultTolerance;
use timsseek::data_sources::speclib::Speclib;
use timsseek::digest::digestion::{
    DigestionEnd,
    DigestionParameters,
    DigestionPattern,
};
use timsseek::errors::TimsSeekError;
use timsseek::fragment_mass::elution_group_converter::SequenceToElutionGroupConverter;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;
use timsseek::models::{
    DigestSlice,
    DigestedSequenceIterator,
    NamedQueryChunk,
    deduplicate_digests,
};
use timsseek::protein::fasta::ProteinSequenceCollection;
use timsseek::scoring::calculate_scores::PreScore;
use timsseek::scoring::search_results::{
    IonSearchResults,
    SearchResultBuilder,
    write_results_to_csv,
};
use tracing::info;

pub fn process_chunk<'a>(
    queries: NamedQueryChunk,
    index: &'a QuadSplittedTransposedIndex,
    factory: &'a MultiCMGStatsFactory<SafePosition>,
    tolerance: &'a DefaultTolerance,
    ref_time_ms: Arc<[u32]>,
) -> Vec<IonSearchResults> {
    let start = Instant::now();
    let num_queries = queries.len();
    let res = query_multi_group(index, tolerance, &queries.queries, &|x| {
        factory.build_with_elution_group(x)
    });
    let elap_time = start.elapsed();
    info!("Querying took {:?}", elap_time);

    let start = Instant::now();

    let tmp: Vec<(IonSearchResults, f32)> = res
        .into_par_iter()
        .zip(queries.into_zip_par_iter())
        .map(
            |(res_elem, (expect_inten, (eg_elem, (digest, charge_elem))))| {
                let builder = SearchResultBuilder::default();
                let prescore = PreScore {
                    charge: charge_elem,
                    digest: &digest,
                    reference: &eg_elem,
                    expected_intensities: &expect_inten,
                    query_values: &res_elem,
                    ref_time_ms: ref_time_ms.clone(),
                };

                let loc = prescore.localize();
                if loc.is_err() {
                    // TODO: Implement filtering out queries that cannot match the data
                    // So we dont get here, to a point where queries can be empty bc no data
                    // can match them.
                    tracing::debug!(
                        "Error localizing pre score: id={:} eg: {:#?}, because of: {:#?}",
                        res_elem.id,
                        eg_elem,
                        loc
                    );
                    return None;
                }
                let loc = loc.unwrap();
                let res = builder.with_localized_pre_score(&loc).finalize();
                if res.is_err() {
                    tracing::error!(
                        "Error creating search result for Digest: {:#?} \nElutionGroup: {:#?}\n Error: {:?}",
                        digest,
                        eg_elem,
                        res,
                    );
                    return None;
                }
                let res = res.unwrap();
                let main_score = res.main_score;
                Some((res, main_score))
            },
        )
        .flatten()
        .collect();

    if tmp.is_empty() {
        // TODO: Remove this and check the error elsewhere.
        panic!("No results found");
    }

    let (out, main_scores): (Vec<IonSearchResults>, Vec<f32>) = tmp.into_iter().unzip();

    let avg_main_scores = main_scores.iter().sum::<f32>() / main_scores.len() as f32;

    assert!(!avg_main_scores.is_nan());
    let elapsed = start.elapsed();
    tracing::info!(
        "Bundling took {:?} for {} elution_groups",
        elapsed,
        num_queries,
    );
    tracing::info!("Avg main score: {:?}", avg_main_scores);

    out
}

pub fn main_loop<'a>(
    chunked_query_iterator: impl ExactSizeIterator<Item = NamedQueryChunk>,
    index: &'a QuadSplittedTransposedIndex,
    factory: &'a MultiCMGStatsFactory<SafePosition>,
    tolerance: &'a DefaultTolerance,
    ref_time_ms: Arc<[u32]>,
    out_path: &Path,
) -> std::result::Result<(), TimsSeekError> {
    let mut chunk_num = 0;
    let mut nqueries = 0;
    let start = Instant::now();

    let style = ProgressStyle::with_template(
        "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
    )
    .unwrap();
    chunked_query_iterator
        .progress_with_style(style)
        .for_each(|chunk| {
            let out = process_chunk(chunk, index, factory, tolerance, ref_time_ms.clone());
            nqueries += out.len();
            let out_path = out_path.join(format!("chunk_{}.csv", chunk_num));
            write_results_to_csv(&out, out_path).unwrap();
            chunk_num += 1;
        });
    let elap_time = start.elapsed();
    println!("Querying took {:?} for {} queries", elap_time, nqueries);
    Ok(())
}

pub fn process_fasta(
    path: PathBuf,
    index: &QuadSplittedTransposedIndex, // TODO: Make generic
    ref_time_ms: Arc<[u32]>,
    factory: &MultiCMGStatsFactory<SafePosition>,
    digestion: DigestionConfig,
    analysis: &AnalysisConfig,
    output: &OutputConfig,
) -> std::result::Result<(), TimsSeekError> {
    let digestion_params = DigestionParameters {
        min_length: digestion.min_length as usize,
        max_length: digestion.max_length as usize,
        pattern: DigestionPattern::trypsin(),
        digestion_end: DigestionEnd::CTerm,
        max_missed_cleavages: digestion.max_missed_cleavages as usize,
    };

    println!(
        "Digesting {} with parameters: \n {:?}",
        path.display(),
        digestion_params
    );

    let fasta_proteins = match ProteinSequenceCollection::from_fasta_file(&path) {
        Ok(x) => x,
        Err(e) => {
            return Err(TimsSeekError::Io {
                source: e,
                path: Some(path),
            });
        }
    };
    let sequences: Vec<Arc<str>> = fasta_proteins
        .sequences
        .iter()
        .map(|x| x.sequence.clone())
        .collect();

    let digest_sequences: Vec<DigestSlice> =
        deduplicate_digests(digestion_params.digest_multiple(&sequences));

    let def_converter = SequenceToElutionGroupConverter::default();
    let chunked_query_iterator = DigestedSequenceIterator::new(
        digest_sequences,
        analysis.chunk_size,
        def_converter,
        digestion.build_decoys,
    );

    main_loop(
        chunked_query_iterator,
        index,
        factory,
        &analysis.tolerance,
        ref_time_ms.clone(),
        &output.directory,
    )?;
    Ok(())
}

pub fn process_speclib(
    path: PathBuf,
    index: &QuadSplittedTransposedIndex,
    ref_time_ms: Arc<[u32]>,
    factory: &MultiCMGStatsFactory<SafePosition>,
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

    main_loop(
        speclib_iter,
        index,
        factory,
        &analysis.tolerance,
        ref_time_ms.clone(),
        &output.directory,
    )?;
    Ok(())
}
