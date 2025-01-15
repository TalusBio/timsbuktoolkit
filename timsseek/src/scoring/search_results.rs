use super::calculate_scores::{
    LocalizedPreScore,
    MainScore,
    PreScore,
    SortedErrors,
};
use crate::errors::Result;
use crate::fragment_mass::fragment_mass_builder::SafePosition;
use crate::models::{
    DecoyMarking,
    DigestSlice,
};
use csv::WriterBuilder;
use serde::Serialize;
use std::path::Path;
use std::time::Instant;
use timsquery::ElutionGroup;

#[derive(Debug, Default)]
pub struct SearchResultBuilder<'q> {
    digest_slice: Option<&'q DigestSlice>,
    ref_eg: Option<&'q ElutionGroup<SafePosition>>,
    decoy_marking: Option<DecoyMarking>,
    charge: Option<u8>,
    nqueries: Option<u8>,

    main_score: Option<f32>,
    delta_next: Option<f32>,
    rt_seconds: Option<f32>,
    observed_mobility: Option<f32>,

    npeaks: Option<u8>,
    lazyerscore: Option<f32>,
    lazyerscore_vs_baseline: Option<f32>,
    norm_lazyerscore_vs_baseline: Option<f32>,
    ms2_cosine_ref_similarity: Option<f32>,
    ms2_coelution_score: Option<f32>,
    ms2_summed_transition_intensity: Option<f32>,

    ms2_mz_errors: Option<[f32; 7]>,
    ms2_mobility_errors: Option<[f32; 7]>,

    ms1_cosine_ref_similarity: Option<f32>,
    ms1_coelution_score: Option<f32>,
    ms1_summed_precursor_intensity: Option<f32>,

    ms1_mz_errors: Option<[f32; 3]>,
    ms1_mobility_errors: Option<[f32; 3]>,
}

impl<'q> SearchResultBuilder<'q> {
    pub fn with_localized_pre_score(self, pre_score: &'q LocalizedPreScore) -> Self {
        self.with_pre_score(&pre_score.pre_score)
            .with_sorted_errors(&pre_score.inten_sorted_errors())
            .with_main_score(pre_score.main_score)
    }

    fn with_pre_score(mut self, pre_score: &'q PreScore) -> Self {
        self.digest_slice = Some(pre_score.digest);
        self.ref_eg = Some(pre_score.reference);
        self.nqueries = Some(pre_score.reference.fragment_mzs.len() as u8);
        self.decoy_marking = Some(pre_score.digest.decoy);
        self.charge = Some(pre_score.charge);
        self
    }

    fn with_sorted_errors(mut self, sorted_errors: &SortedErrors) -> Self {
        self.ms1_mz_errors = Some(sorted_errors.ms1_mz_errors);
        self.ms1_mobility_errors = Some(sorted_errors.ms1_mobility_errors);
        self.ms2_mz_errors = Some(sorted_errors.ms2_mz_errors);
        self.ms2_mobility_errors = Some(sorted_errors.ms2_mobility_errors);
        self
    }

    fn with_main_score(mut self, main_score: MainScore) -> Self {
        // TODO use exhaustive unpacking to make more explicit any values
        // I Might be ignoring.
        let MainScore {
            score,
            delta_next,
            observed_mobility,
            ms2_cosine_ref_sim,
            ms2_coelution_score,
            ms1_coelution_score,
            ms1_cosine_ref_sim,
            lazyscore,
            lazyscore_vs_baseline,
            lazyscore_z,
            npeaks,
            ms1_summed_intensity,
            ms2_summed_intensity,
            retention_time_ms,
            ..
        } = main_score;
        {
            self.main_score = Some(score);
            self.delta_next = Some(delta_next);
            self.observed_mobility = Some(observed_mobility);
            self.rt_seconds = Some(retention_time_ms as f32 / 1000.0);
            self.ms2_cosine_ref_similarity = Some(ms2_cosine_ref_sim);
            self.ms2_coelution_score = Some(ms2_coelution_score);
            self.ms1_coelution_score = Some(ms1_coelution_score);
            self.ms1_cosine_ref_similarity = Some(ms1_cosine_ref_sim);
            self.lazyerscore = Some(lazyscore);
            self.lazyerscore_vs_baseline = Some(lazyscore_vs_baseline);
            self.norm_lazyerscore_vs_baseline = Some(lazyscore_z);
            self.npeaks = Some(npeaks);
            self.ms1_summed_precursor_intensity = Some(ms1_summed_intensity);
            self.ms2_summed_transition_intensity = Some(ms2_summed_intensity);
        }

        self
    }

    pub fn finalize(self) -> Result<IonSearchResults> {
        let [mz1_e0, mz1_e1, mz1_e2] = self.ms1_mz_errors.unwrap();
        let [mz2_e0, mz2_e1, mz2_e2, mz2_e3, mz2_e4, mz2_e5, mz2_e6] = self.ms2_mz_errors.unwrap();

        let [mob1_e0, mob1_e1, mob1_e2] = self.ms1_mobility_errors.unwrap();
        let [
            mob2_e0,
            mob2_e1,
            mob2_e2,
            mob2_e3,
            mob2_e4,
            mob2_e5,
            mob2_e6,
        ] = self.ms2_mobility_errors.unwrap();

        let ref_eg = self.ref_eg.unwrap();
        // TODO replace this with exhaustive unpacking.

        let results = IonSearchResults {
            sequence: String::from(self.digest_slice.unwrap().clone()),
            precursor_mz: ref_eg.precursor_mzs[1],
            precursor_charge: self.charge.unwrap(),
            precursor_mobility_query: ref_eg.mobility,
            precursor_rt_query_seconds: ref_eg.rt_seconds,
            nqueries: self.nqueries.unwrap(),
            is_target: self.decoy_marking.unwrap().is_target(),
            main_score: self.main_score.unwrap(),
            delta_next: self.delta_next.unwrap(),
            obs_rt_seconds: self.rt_seconds.unwrap(),
            obs_mobility: self.observed_mobility.unwrap(),
            npeaks: self.npeaks.unwrap(),
            lazyerscore: self.lazyerscore.unwrap(),
            lazyerscore_vs_baseline: self.lazyerscore_vs_baseline.unwrap(),
            norm_lazyerscore_vs_baseline: self.norm_lazyerscore_vs_baseline.unwrap(),
            ms2_cosine_ref_similarity: self.ms2_cosine_ref_similarity.unwrap(),
            ms2_summed_transition_intensity: self.ms2_summed_transition_intensity.unwrap(),
            ms2_mz_error_0: mz2_e0,
            ms2_mz_error_1: mz2_e1,
            ms2_mz_error_2: mz2_e2,
            ms2_mz_error_3: mz2_e3,
            ms2_mz_error_4: mz2_e4,
            ms2_mz_error_5: mz2_e5,
            ms2_mz_error_6: mz2_e6,
            ms2_mobility_error_0: mob2_e0,
            ms2_mobility_error_1: mob2_e1,
            ms2_mobility_error_2: mob2_e2,
            ms2_mobility_error_3: mob2_e3,
            ms2_mobility_error_4: mob2_e4,
            ms2_mobility_error_5: mob2_e5,
            ms2_mobility_error_6: mob2_e6,
            ms2_coelution_score: self.ms2_coelution_score.unwrap(),
            ms1_cosine_ref_similarity: self.ms1_cosine_ref_similarity.unwrap(),
            ms1_summed_precursor_intensity: self.ms1_summed_precursor_intensity.unwrap(),
            ms1_mz_error_0: mz1_e0,
            ms1_mz_error_1: mz1_e1,
            ms1_mz_error_2: mz1_e2,
            ms1_coelution_score: self.ms1_coelution_score.unwrap(),
            ms1_mobility_error_0: mob1_e0,
            ms1_mobility_error_1: mob1_e1,
            ms1_mobility_error_2: mob1_e2,
        };

        Ok(results)
    }
}

#[derive(Debug, Serialize, Clone)]
pub struct IonSearchResults {
    sequence: String,
    precursor_mz: f64,
    precursor_charge: u8,
    precursor_mobility_query: f32,
    precursor_rt_query_seconds: f32,
    nqueries: u8,
    is_target: bool,

    // Combined
    pub main_score: f32,
    pub delta_next: f32,
    obs_rt_seconds: f32,
    obs_mobility: f32,

    // MS2
    npeaks: u8,
    lazyerscore: f32,
    lazyerscore_vs_baseline: f32,
    norm_lazyerscore_vs_baseline: f32,
    ms2_cosine_ref_similarity: f32,
    ms2_coelution_score: f32,
    ms2_summed_transition_intensity: f32,

    // MS2 - Split
    // Flattening manually bc serde(flatten)
    // is not supported by csv ...
    // https://github.com/BurntSushi/rust-csv/pull/223
    ms2_mz_error_0: f32,
    ms2_mz_error_1: f32,
    ms2_mz_error_2: f32,
    ms2_mz_error_3: f32,
    ms2_mz_error_4: f32,
    ms2_mz_error_5: f32,
    ms2_mz_error_6: f32,
    ms2_mobility_error_0: f32,
    ms2_mobility_error_1: f32,
    ms2_mobility_error_2: f32,
    ms2_mobility_error_3: f32,
    ms2_mobility_error_4: f32,
    ms2_mobility_error_5: f32,
    ms2_mobility_error_6: f32,

    // MS1
    ms1_cosine_ref_similarity: f32,
    ms1_coelution_score: f32,
    ms1_summed_precursor_intensity: f32,

    // MS1 Split
    ms1_mz_error_0: f32,
    ms1_mz_error_1: f32,
    ms1_mz_error_2: f32,
    ms1_mobility_error_0: f32,
    ms1_mobility_error_1: f32,
    ms1_mobility_error_2: f32,
}

pub fn write_results_to_csv<P: AsRef<Path>>(
    results: &[IonSearchResults],
    out_path: P,
) -> std::result::Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let mut writer = WriterBuilder::default()
        .has_headers(true)
        .from_path(out_path.as_ref())?;

    for result in results {
        writer.serialize(result).unwrap();
    }
    writer.flush()?;
    tracing::info!(
        "Writing took {:?} -> {:?}",
        start.elapsed(),
        out_path.as_ref()
    );
    Ok(())
}
