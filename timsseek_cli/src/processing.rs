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
use std::time::{
    Duration,
    Instant,
};
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
use timsseek::scoring::calculate_scores::{
    LocalizedPreScore,
    PreScore,
};
use timsseek::scoring::search_results::{
    IonSearchResults,
    SearchResultBuilder,
    write_results_to_csv,
};
use tracing::info;

#[derive(Debug, Clone, Default, Copy)]
pub struct RuntimeMetrics {
    pub num_queries: usize,
    pub num_skipped: usize,
    pub query_time: Duration,
    pub total_time: Duration,
    pub time_localizing: Duration,
    pub time_scoring: Duration,
}

pub struct RuntimeThroughput {
    pub query_throughput_persec: f64,
    pub total_throughput_persec: f64,
    pub localizing_throughput_persec: f64,
    pub scoring_throughput_persec: f64,
}

impl RuntimeMetrics {
    pub fn throughput(&self) -> RuntimeThroughput {
        RuntimeThroughput {
            query_throughput_persec: self.num_queries as f64 / self.query_time.as_secs_f64(),
            total_throughput_persec: self.num_queries as f64 / self.total_time.as_secs_f64(),
            localizing_throughput_persec: self.num_queries as f64
                / self.time_localizing.as_secs_f64(),
            scoring_throughput_persec: self.num_queries as f64 / self.time_scoring.as_secs_f64(),
        }
    }

    fn fold(self, other: Self) -> Self {
        Self {
            num_queries: self.num_queries + other.num_queries,
            num_skipped: self.num_skipped + other.num_skipped,
            query_time: self.query_time + other.query_time,
            total_time: self.total_time + other.total_time,
            time_localizing: self.time_localizing + other.time_localizing,
            time_scoring: self.time_scoring + other.time_scoring,
        }
    }
}

impl std::fmt::Display for RuntimeThroughput {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Runtime throughput: ")?;
        write!(f, "Queries/sec: {:#?}", self.query_throughput_persec)?;
        write!(f, " Total/sec: {:#?}", self.total_throughput_persec)?;
        write!(
            f,
            " Localizing/sec: {:#?} ",
            self.localizing_throughput_persec
        )?;
        write!(f, " Scoring/sec: {:#?}", self.scoring_throughput_persec)?;
        Ok(())
    }
}

impl std::fmt::Display for RuntimeMetrics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Runtime metrics: ")?;
        write!(f, "Queries: {}", self.num_queries)?;
        write!(f, " Skipped: {}", self.num_skipped)?;
        write!(f, " Query time: {:#?}", self.query_time)?;
        write!(f, " Total time: {:#?}", self.total_time)?;
        write!(f, " Localizing time: {:#?}", self.time_localizing)?;
        write!(f, " Scoring time: {:#?}", self.time_scoring)?;
        Ok(())
    }
}

pub fn process_chunk<'a>(
    queries: NamedQueryChunk,
    index: &'a QuadSplittedTransposedIndex,
    factory: &'a MultiCMGStatsFactory<SafePosition>,
    tolerance: &'a DefaultTolerance,
    ref_time_ms: Arc<[u32]>,
) -> (Vec<IonSearchResults>, RuntimeMetrics) {
    let start = Instant::now();
    let num_queries = queries.len();
    let res = query_multi_group(index, tolerance, &queries.queries, &|x| {
        factory.build_with_elution_group(x)
    });
    let elap_time_querying = start.elapsed();
    info!("Querying took {:?}", elap_time_querying);

    let loc_start = Instant::now();
    let loc_scores: Vec<LocalizedPreScore> = res
        .into_par_iter()
        .zip(queries.into_zip_par_iter())
        .filter_map(
            |(res_elem, (expect_inten, (eg_elem, (digest, charge_elem))))| {
                let id = res_elem.id;
                let prescore = PreScore {
                    charge: charge_elem,
                    digest: digest,
                    reference: eg_elem,
                    expected_intensities: expect_inten,
                    query_values: res_elem,
                    ref_time_ms: ref_time_ms.clone(),
                };

                let loc = prescore.localize();
                match loc {
                    Ok(loc) => Some(loc),
                    Err(e) => {
                        tracing::debug!(
                            "Error localizing pre score: id={:}, because of: {:#?}",
                            id,
                            e
                        );
                        return None;
                    }
                }
            },
        )
        .collect();
    let elap_localizing = loc_start.elapsed();

    let start_scoring = Instant::now();
    let tmp: Vec<(IonSearchResults, f32)> = loc_scores
        .into_par_iter()
        .map(|loc| {
                let builder = SearchResultBuilder::default();
                let res = builder.with_localized_pre_score(&loc).finalize();
                if res.is_err() {
                    tracing::error!(
                        "Error creating search result for Digest: {:#?} \nElutionGroup: {:#?}\n Error: {:?}",
                        loc.pre_score.digest,
                        loc.pre_score.reference,
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

    let num_skipped = num_queries - tmp.len();
    let elap_scoring = start_scoring.elapsed();

    if tmp.is_empty() {
        // TODO: Remove this and check the error elsewhere.
        panic!("No results found");
    }

    let (out, main_scores): (Vec<IonSearchResults>, Vec<f32>) = tmp.into_iter().unzip();

    let avg_main_scores = main_scores.iter().sum::<f32>() / main_scores.len() as f32;

    assert!(!avg_main_scores.is_nan());
    let elapsed = start.elapsed();
    let metrics = RuntimeMetrics {
        num_queries,
        num_skipped,
        query_time: elap_time_querying,
        total_time: elapsed,
        time_localizing: elap_localizing,
        time_scoring: elap_scoring,
    };
    tracing::info!("{} {}", metrics, metrics.throughput());
    tracing::info!("Avg main score: {:?}", avg_main_scores);

    (out, metrics)
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
    let start = Instant::now();
    let mut metrics = RuntimeMetrics::default();

    let style = ProgressStyle::with_template(
        "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
    )
    .unwrap();
    chunked_query_iterator
        .progress_with_style(style)
        .for_each(|chunk| {
            let out = process_chunk(chunk, index, factory, tolerance, ref_time_ms.clone());
            metrics = metrics.fold(out.1);
            let out_path = out_path.join(format!("chunk_{}.csv", chunk_num));
            write_results_to_csv(&out.0, out_path).unwrap();
            chunk_num += 1;
        });
    println!(
        "Finished processing {} chunks in {:?}",
        chunk_num,
        start.elapsed()
    );
    println!("{} {}", metrics, metrics.throughput());
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
