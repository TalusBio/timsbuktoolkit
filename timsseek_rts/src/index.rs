use std::collections::HashMap;
use std::sync::Arc;

use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::NaturalFinalizedMultiCMGArrays;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsquery::models::aggregators::MultiCMGStatsFactory;
use timsquery::queriable_tims_data::queriable_tims_data::query_multi_group;
use timsquery::traits::tolerance::DefaultTolerance;
use timsseek::data_sources::speclib::ExpectedIntensities;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;
use timsseek::models::DigestSlice;
use timsseek::scoring::calculate_scores::{IntensityArrays, LongitudinalMainScoreElements, PreScore};
use timsseek::scoring::search_results::{IonSearchResults, SearchResultBuilder};
use timsseek::utils::tdf::get_ms1_frame_times_ms;
use timsquery::ElutionGroup;
use timsseek::errors::Result;
use serde::{Deserialize, Serialize};

pub struct BundledDotDIndex {
    index: QuadSplittedTransposedIndex,
    ref_time_ms: Arc<[u32]>,
    factory: MultiCMGStatsFactory<SafePosition>,
    tolerance: DefaultTolerance,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputQuery {
    pub sequence: String,
    pub charge: u8,
    pub elution_group: ElutionGroup<SafePosition>,
    pub expected_intensities: ExpectedIntensities,
}

impl InputQuery {
    pub fn sample() -> Self {
        Self {
            sequence: "PEPTIDE".to_string(),
            charge: 2,
            elution_group: ElutionGroup {
                id: 0,
                mobility: 0.75,
                rt_seconds: 0.0,
                precursor_mzs: vec![450.0, 450.5, 451.0, 451.5],
                fragment_mzs: HashMap::from_iter([
                    (SafePosition::from_str("a1").unwrap(), 450.0),
                    (SafePosition::from_str("a2").unwrap(), 450.5),
                    (SafePosition::from_str("a3").unwrap(), 451.0),
                    (SafePosition::from_str("a4").unwrap(), 451.5),
                ].iter().cloned()),
            },
            expected_intensities: ExpectedIntensities {
                precursor_intensities: vec![1.0, 1.0, 1.0, 1.0],
                fragment_intensities: HashMap::from_iter([
                    (SafePosition::from_str("a1").unwrap(), 1.0),
                    (SafePosition::from_str("a2").unwrap(), 1.0),
                    (SafePosition::from_str("a3").unwrap(), 1.0),
                    (SafePosition::from_str("a4").unwrap(), 1.0),
                ].iter().cloned()),
            }
        }
    }
}

impl From<InputQuery> for NamedQuery {
    fn from(value: InputQuery) -> Self {
        Self {
            digest: DigestSlice::from_string(value.sequence.clone(), false),
            charge: value.charge,
            elution_group: value.elution_group,
            expected_intensities: value.expected_intensities,
        }
    }
}

#[derive(Debug, Clone)]
pub struct NamedQuery {
    pub digest: DigestSlice,
    pub charge: u8,
    pub elution_group: ElutionGroup<SafePosition>,
    pub expected_intensities: ExpectedIntensities,
}

#[derive(Debug, Clone, Serialize)]
pub struct QueryResult {
    pub main_score_elements: LongitudinalMainScoreElements,
    pub extractions: NaturalFinalizedMultiCMGArrays<SafePosition>,
    pub search_results: IonSearchResults,
}

impl BundledDotDIndex {
    pub fn new(dotd_file_location: std::path::PathBuf, tolerance: DefaultTolerance) -> Result<BundledDotDIndex> {
        let index = QuadSplittedTransposedIndex::from_path_centroided(
            dotd_file_location
                .clone()
                .to_str()
                .expect("Path is not convertable to string"),
        )?;

        let tdf_path = &dotd_file_location.clone().join("analysis.tdf");
        let ref_time_ms = get_ms1_frame_times_ms(tdf_path.to_str().unwrap()).unwrap();

        let factory = MultiCMGStatsFactory {
            converters: (index.mz_converter, index.im_converter),
            _phantom: std::marker::PhantomData::<SafePosition>,
        };
        
        Ok(BundledDotDIndex {
            index,
            ref_time_ms,
            factory,
            tolerance,
        })
    }

    pub fn query(&self, queries: NamedQuery) -> Result<QueryResult> {
        let res = query_multi_group(&self.index, &self.tolerance, &[queries.elution_group.clone()], &|x| {
            self.factory.build_with_elution_group(x)
        });
        let builder = SearchResultBuilder::default();
        let prescore = PreScore {
            charge: queries.charge.clone(),
            digest: &queries.digest,
            reference: &queries.elution_group,
            expected_intensities: &queries.expected_intensities,
            query_values: &res[0],
            ref_time_ms: self.ref_time_ms.clone(),
        };

        let int_arrs = IntensityArrays::new(&res[0], &queries.expected_intensities);
        let longitudinal_main_score_elements = LongitudinalMainScoreElements::new(
            &int_arrs,
            self.ref_time_ms.clone(),
            &res[0].ms1_arrays.retention_time_miliseconds,
            &res[0].ms2_arrays.retention_time_miliseconds,
        );

        let res2 = builder.with_localized_pre_score(&prescore.localize()).finalize()?;

        Ok(QueryResult {
            main_score_elements: longitudinal_main_score_elements,
            extractions: res[0].clone(),
            search_results: res2,
        })
    }
}