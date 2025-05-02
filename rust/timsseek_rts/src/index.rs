use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;

use serde::{
    Deserialize,
    Serialize,
};
use timsquery::models::aggregators::EGCAggregator;
use timsquery::models::indices::ExpandedRawFrameIndex;
use timsquery::{
    ElutionGroup,
    QueriableData,
    Tolerance,
};
use timsseek::data_sources::speclib::ExpectedIntensities;
use timsseek::errors::Result;
use timsseek::fragment_mass::IonAnnot;
use timsseek::models::DigestSlice;
use timsseek::scoring::calculate_scores::{
    IntensityArrays,
    LongitudinalMainScoreElements,
    PreScore,
};
use timsseek::scoring::full_results::FullQueryResult;
use timsseek::scoring::search_results::SearchResultBuilder;
use timsseek::utils::tdf::get_ms1_frame_times_ms;
use timsseek::scoring::{QueryItemToScore, Scorer};

// TODO: replace with a trait ... This works for now though
type IndexUse = ExpandedRawFrameIndex;

pub struct BundledDotDIndex {
    index: IndexUse,
    tolerance: Tolerance,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputQuery {
    pub sequence: String,
    pub charge: u8,
    pub elution_group: Arc<ElutionGroup<IonAnnot>>,
    pub expected_intensities: ExpectedIntensities,
}

impl InputQuery {
    pub fn sample() -> Self {
        Self {
            sequence: "PEPTIDE".to_string(),
            charge: 2,
            elution_group: Arc::new(ElutionGroup {
                id: 0,
                mobility: 0.75,
                rt_seconds: 0.0,
                precursors: vec![(-1, 450.0), (0, 450.5), (1, 451.0), (2, 451.5)].into(),
                fragments: vec![
                    (IonAnnot::try_from("a1").unwrap(), 450.0),
                    (IonAnnot::try_from("a2").unwrap(), 450.5),
                    (IonAnnot::try_from("a3").unwrap(), 451.0),
                    (IonAnnot::try_from("a4").unwrap(), 451.5),
                ]
                .into(),
            }),
            expected_intensities: ExpectedIntensities {
                precursor_intensities: vec![1.0, 1.0, 1.0, 1.0],
                fragment_intensities: HashMap::from_iter(
                    [
                        (IonAnnot::try_from("a1").unwrap(), 1.0),
                        (IonAnnot::try_from("a2").unwrap(), 1.0),
                        (IonAnnot::try_from("a3").unwrap(), 1.0),
                        (IonAnnot::try_from("a4").unwrap(), 1.0),
                    ]
                    .iter()
                    .cloned(),
                ),
            },
        }
    }
}

impl From<InputQuery> for QueryItemToScore {
    fn from(value: InputQuery) -> Self {
        Self {
            digest: DigestSlice::from_string(value.sequence.clone(), false),
            charge: value.charge,
            query: value.elution_group,
            expected_intensity: value.expected_intensities,
        }
    }
}

#[derive(Debug, Clone)]
pub struct NamedQuery {
    pub digest: DigestSlice,
    pub charge: u8,
    pub elution_group: ElutionGroup<IonAnnot>,
    pub expected_intensities: ExpectedIntensities,
}

pub fn new_index(
    dotd_file_location: std::path::PathBuf,
    tolerance: Tolerance,
) -> Result<Scorer<IndexUse>> {
    let st = Instant::now();
    // Can use centroided for faster queries ...
    let index = ExpandedRawFrameIndex::from_path(
        dotd_file_location
            .clone()
            .to_str()
            .expect("Path is not convertable to string"),
    )?;
    let elap_time = st.elapsed();
    println!(
        "Loading index took: {:?} for {}",
        elap_time,
        dotd_file_location.display()
    );

    let tdf_path = &dotd_file_location.clone().join("analysis.tdf");
    let ref_time_ms = get_ms1_frame_times_ms(tdf_path.to_str().unwrap()).unwrap();

    Ok(Scorer {
        index_cycle_rt_ms: ref_time_ms,
        index,
        tolerance,
    })
}

impl BundledDotDIndex {

    pub fn query(&self, queries: NamedQuery) -> Result<FullQueryResult> {
        let mut res = EGCAggregator::new(
            Arc::new(queries.elution_group.clone()),
            self.index.cycle_rt_ms.clone(),
        )
        .unwrap();
        self.index.add_query(&mut res, &self.tolerance);
        let builder = SearchResultBuilder::default();
        let int_arrs = IntensityArrays::new(&res, &queries.expected_intensities)?;
        let prescore = PreScore {
            charge: queries.charge,
            digest: queries.digest,
            expected_intensities: queries.expected_intensities,
            query_values: res.clone(),
        };

        let longitudinal_main_score_elements = LongitudinalMainScoreElements::new(&int_arrs)?;

        let res2 = builder
            .with_localized_pre_score(&prescore.localize()?)
            .finalize()?;
        let longitudinal_main_score = longitudinal_main_score_elements.main_score_iter().collect();

        Ok(FullQueryResult {
            main_score_elements: longitudinal_main_score_elements,
            longitudinal_main_score,
            extractions: res.clone(),
            search_results: res2,
        })
    }
}
