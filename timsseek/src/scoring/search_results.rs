use super::calculate_scores::{
    LocalizedPreScore,
    MainScore,
    PreScore,
    RelativeIntensities,
    SortedErrors,
};
use crate::errors::DataProcessingError;
use crate::fragment_mass::fragment_mass_builder::SafePosition;
use crate::models::{
    DecoyMarking,
    DigestSlice,
};
use csv::WriterBuilder;
use parquet::file::writer::SerializedFileWriter;
use parquet::record::RecordWriter;
use serde::Serialize;
use std::fs::File;
use std::path::Path;
use std::time::Instant;
use timsquery::ElutionGroup;

#[derive(Debug, Default)]
pub struct SearchResultBuilder<'q> {
    digest_slice: SetField<&'q DigestSlice>,
    ref_eg: SetField<&'q ElutionGroup<SafePosition>>,
    decoy_marking: SetField<DecoyMarking>,
    charge: SetField<u8>,
    nqueries: SetField<u8>,

    main_score: SetField<f32>,
    delta_next: SetField<f32>,
    rt_seconds: SetField<f32>,
    observed_mobility: SetField<f32>,
    delta_ms1_ms2_mobility: SetField<f32>,
    // ms1_ms2_correlation: SetField<f32>,
    npeaks: SetField<u8>,
    lazyerscore: SetField<f32>,
    lazyerscore_vs_baseline: SetField<f32>,
    norm_lazyerscore_vs_baseline: SetField<f32>,
    ms2_cosine_ref_similarity: SetField<f32>,
    ms2_coelution_score: SetField<f32>,
    ms2_summed_transition_intensity: SetField<f32>,

    ms2_mz_errors: SetField<[f32; 7]>,
    ms2_mobility_errors: SetField<[f32; 7]>,

    ms1_cosine_ref_similarity: SetField<f32>,
    ms1_coelution_score: SetField<f32>,
    ms1_summed_precursor_intensity: SetField<f32>,

    ms1_mz_errors: SetField<[f32; 3]>,
    ms1_mobility_errors: SetField<[f32; 3]>,

    relative_intensities: SetField<RelativeIntensities>,
}

#[derive(Debug, Clone, Copy, Default)]
pub enum SetField<T> {
    Some(T),
    #[default]
    None,
}

impl<T> SetField<T> {
    pub fn is_some(&self) -> bool {
        matches!(self, Self::Some(_))
    }

    pub fn expect_some(
        self,
        field_name: impl ToString,
        msg: impl ToString,
    ) -> Result<T, DataProcessingError> {
        match self {
            Self::Some(v) => Ok(v),
            Self::None => Err(DataProcessingError::ExpectedSetField {
                field: field_name.to_string(),
                context: msg.to_string(),
            }),
        }
    }
}

impl<'q> SearchResultBuilder<'q> {
    pub fn with_localized_pre_score(self, pre_score: &'q LocalizedPreScore) -> Self {
        self.with_pre_score(&pre_score.pre_score)
            .with_sorted_errors(pre_score.inten_sorted_errors())
            .with_relative_intensities(pre_score.relative_intensities())
            .with_main_score(pre_score.main_score)
    }

    fn with_pre_score(mut self, pre_score: &'q PreScore) -> Self {
        self.digest_slice = SetField::Some(&pre_score.digest);
        self.ref_eg = SetField::Some(&pre_score.reference);
        self.nqueries = SetField::Some(pre_score.reference.fragment_mzs.len() as u8);
        self.decoy_marking = SetField::Some(pre_score.digest.decoy);
        self.charge = SetField::Some(pre_score.charge);
        self
    }

    fn with_sorted_errors(mut self, sorted_errors: SortedErrors) -> Self {
        self.ms1_mz_errors = SetField::Some(sorted_errors.ms1_mz_errors);
        self.ms1_mobility_errors = SetField::Some(sorted_errors.ms1_mobility_errors);
        self.ms2_mz_errors = SetField::Some(sorted_errors.ms2_mz_errors);
        self.ms2_mobility_errors = SetField::Some(sorted_errors.ms2_mobility_errors);
        self
    }

    fn with_relative_intensities(mut self, relative_intensities: RelativeIntensities) -> Self {
        self.relative_intensities = SetField::Some(relative_intensities);
        self
    }

    fn with_main_score(mut self, main_score: MainScore) -> Self {
        // TODO use exhaustive unpacking to make more explicit any values
        // I Might be ignoring.
        let MainScore {
            score,
            delta_next,
            observed_mobility,
            observed_mobility_ms1,
            observed_mobility_ms2,
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
            let delta_ms1_ms2_mobility = observed_mobility_ms1 - observed_mobility_ms2;
            self.main_score = SetField::Some(score);
            self.delta_next = SetField::Some(delta_next);
            self.observed_mobility = SetField::Some(observed_mobility);
            self.delta_ms1_ms2_mobility = SetField::Some(delta_ms1_ms2_mobility);

            self.rt_seconds = SetField::Some(retention_time_ms as f32 / 1000.0);
            self.ms2_cosine_ref_similarity = SetField::Some(ms2_cosine_ref_sim);
            self.ms2_coelution_score = SetField::Some(ms2_coelution_score);
            self.ms1_coelution_score = SetField::Some(ms1_coelution_score);
            self.ms1_cosine_ref_similarity = SetField::Some(ms1_cosine_ref_sim);
            self.lazyerscore = SetField::Some(lazyscore);
            self.lazyerscore_vs_baseline = SetField::Some(lazyscore_vs_baseline);
            self.norm_lazyerscore_vs_baseline = SetField::Some(lazyscore_z);
            self.npeaks = SetField::Some(npeaks);
            self.ms1_summed_precursor_intensity = SetField::Some(ms1_summed_intensity);
            self.ms2_summed_transition_intensity = SetField::Some(ms2_summed_intensity);
        }

        self
    }

    pub fn finalize(self) -> Result<IonSearchResults, DataProcessingError> {
        let [mz1_e0, mz1_e1, mz1_e2] = self
            .ms1_mz_errors
            .expect_some("ms1_mz_errors", "ms1_mz_errors")?;
        let [mz2_e0, mz2_e1, mz2_e2, mz2_e3, mz2_e4, mz2_e5, mz2_e6] = self
            .ms2_mz_errors
            .expect_some("ms2_mz_errors", "ms2_mz_errors")?;

        let [mob1_e0, mob1_e1, mob1_e2] = self
            .ms1_mobility_errors
            .expect_some("ms1_mobility_errors", "ms1_mobility_errors")?;
        let [
            mob2_e0,
            mob2_e1,
            mob2_e2,
            mob2_e3,
            mob2_e4,
            mob2_e5,
            mob2_e6,
        ] = self
            .ms2_mobility_errors
            .expect_some("ms2_mobility_errors", "ms2_mobility_errors")?;

        let [int1_e0, int1_e1, int1_e2] = self
            .relative_intensities
            .expect_some("ms1_intensity_errors", "ms1_intensity_errors")?
            .ms1;
        let [
            int2_e0,
            int2_e1,
            int2_e2,
            int2_e3,
            int2_e4,
            int2_e5,
            int2_e6,
        ] = self
            .relative_intensities
            .expect_some("ms2_intensity_errors", "ms2_intensity_errors")?
            .ms2;

        let ref_eg = self.ref_eg.expect_some("ref_eg", "ref_eg")?;
        // TODO replace this with exhaustive unpacking.
        let obs_rt_seconds = self.rt_seconds.expect_some("rt_seconds", "rt_seconds")?;
        let delta_theo_rt = obs_rt_seconds - ref_eg.rt_seconds;
        let sq_delta_theo_rt = delta_theo_rt * delta_theo_rt;

        let delta_ms1_ms2_mobility = self
            .delta_ms1_ms2_mobility
            .expect_some("delta_ms1_ms2_mobility", "delta_ms1_ms2_mobility")?;
        let sq_delta_ms1_ms2_mobility = delta_ms1_ms2_mobility * delta_ms1_ms2_mobility;

        let results = IonSearchResults {
            sequence: String::from(
                self.digest_slice
                    .expect_some("digest_slice", "digest_slice")?
                    .clone(),
            ),
            precursor_mz: ref_eg.precursor_mzs[1],
            precursor_charge: self.charge.expect_some("charge", "charge")?,
            precursor_mobility_query: ref_eg.mobility,
            precursor_rt_query_seconds: ref_eg.rt_seconds,
            nqueries: self.nqueries.expect_some("nqueries", "nqueries")?,
            is_target: self
                .decoy_marking
                .expect_some("decoy_marking", "decoy_marking")?
                .is_target(),
            main_score: self.main_score.expect_some("main_score", "main_score")?,
            delta_next: self.delta_next.expect_some("delta_next", "delta_next")?,
            delta_theo_rt,
            sq_delta_theo_rt,
            obs_rt_seconds,
            obs_mobility: self
                .observed_mobility
                .expect_some("observed_mobility", "observed_mobility")?,
            delta_ms1_ms2_mobility,
            sq_delta_ms1_ms2_mobility,
            npeaks: self.npeaks.expect_some("npeaks", "npeaks")?,
            lazyerscore: self.lazyerscore.expect_some("lazyerscore", "lazyerscore")?,
            lazyerscore_vs_baseline: self
                .lazyerscore_vs_baseline
                .expect_some("lazyerscore_vs_baseline", "lazyerscore_vs_baseline")?,
            // ms1_ms2_correlation: self
            //     .ms1_ms2_correlation
            //     .expect_some("ms1_ms2_correlation", "ms1_ms2_correlation")?,
            norm_lazyerscore_vs_baseline: self.norm_lazyerscore_vs_baseline.expect_some(
                "norm_lazyerscore_vs_baseline",
                "norm_lazyerscore_vs_baseline",
            )?,
            ms2_cosine_ref_similarity: self
                .ms2_cosine_ref_similarity
                .expect_some("ms2_cosine_ref_similarity", "ms2_cosine_ref_similarity")?,
            ms2_summed_transition_intensity: self.ms2_summed_transition_intensity.expect_some(
                "ms2_summed_transition_intensity",
                "ms2_summed_transition_intensity",
            )?,
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
            ms2_coelution_score: self
                .ms2_coelution_score
                .expect_some("ms2_coelution_score", "ms2_coelution_score")?,
            ms1_cosine_ref_similarity: self
                .ms1_cosine_ref_similarity
                .expect_some("ms1_cosine_ref_similarity", "ms1_cosine_ref_similarity")?,
            ms1_summed_precursor_intensity: self.ms1_summed_precursor_intensity.expect_some(
                "ms1_summed_precursor_intensity",
                "ms1_summed_precursor_intensity",
            )?,
            ms1_mz_error_0: mz1_e0,
            ms1_mz_error_1: mz1_e1,
            ms1_mz_error_2: mz1_e2,
            ms1_coelution_score: self
                .ms1_coelution_score
                .expect_some("ms1_coelution_score", "ms1_coelution_score")?,
            ms1_mobility_error_0: mob1_e0,
            ms1_mobility_error_1: mob1_e1,
            ms1_mobility_error_2: mob1_e2,

            ms1_inten_ratio_0: int1_e0,
            ms1_inten_ratio_1: int1_e1,
            ms1_inten_ratio_2: int1_e2,

            ms2_inten_ratio_0: int2_e0,
            ms2_inten_ratio_1: int2_e1,
            ms2_inten_ratio_2: int2_e2,
            ms2_inten_ratio_3: int2_e3,
            ms2_inten_ratio_4: int2_e4,
            ms2_inten_ratio_5: int2_e5,
            ms2_inten_ratio_6: int2_e6,
        };

        Ok(results)
    }
}

#[derive(Debug, Clone, Serialize, ParquetRecordWriter)]
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
    pub obs_rt_seconds: f32,
    obs_mobility: f32,
    delta_theo_rt: f32,
    sq_delta_theo_rt: f32,
    delta_ms1_ms2_mobility: f32,
    // ms1_ms2_correlation: f32,
    sq_delta_ms1_ms2_mobility: f32,

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

    // Relative Intensities
    ms1_inten_ratio_0: f32,
    ms1_inten_ratio_1: f32,
    ms1_inten_ratio_2: f32,

    ms2_inten_ratio_0: f32,
    ms2_inten_ratio_1: f32,
    ms2_inten_ratio_2: f32,
    ms2_inten_ratio_3: f32,
    ms2_inten_ratio_4: f32,
    ms2_inten_ratio_5: f32,
    ms2_inten_ratio_6: f32,
}

pub fn write_results_to_parquet<P: AsRef<Path> + Clone>(
    results: &[IonSearchResults],
    out_path: P,
) -> std::result::Result<(), Box<dyn std::error::Error>> {
    // TODO: Implement multi chunk accumulator
    let file = match File::create_new(out_path.clone()) {
        Ok(file) => file,
        Err(err) => {
            tracing::error!(
                "Failed to open file {:?} with error: {}",
                out_path.as_ref(),
                err
            );
            return Err(Box::new(err));
        }
    };
    let schema = results.schema().unwrap();
    let mut writer = SerializedFileWriter::new(file, schema, Default::default()).unwrap();
    let mut row_group = writer.next_row_group().unwrap();
    results.write_to_row_group(&mut row_group).unwrap();
    row_group.close().unwrap();
    writer.close().unwrap();

    Ok(())
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
