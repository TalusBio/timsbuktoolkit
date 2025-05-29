use super::calculate_scores::{
    MainScore,
    PreScore,
    RelativeIntensities,
};
use super::offsets::MzMobilityOffsets;
use super::scorer::SecondaryLazyScores;
use super::{
    NUM_MS1_IONS,
    NUM_MS2_IONS,
};
use crate::IonAnnot;
use crate::errors::DataProcessingError;
use crate::models::{
    DecoyMarking,
    DigestSlice,
};
use parquet::file::writer::SerializedFileWriter;
use parquet::record::RecordWriter;
use serde::Serialize;
use std::fs::File;
use std::path::Path;
use timsquery::ElutionGroup;
use tracing::debug;

#[derive(Debug, Default)]
pub struct SearchResultBuilder<'q> {
    digest_slice: SetField<&'q DigestSlice>,
    ref_eg: SetField<&'q ElutionGroup<IonAnnot>>,
    decoy_marking: SetField<DecoyMarking>,
    charge: SetField<u8>,
    nqueries: SetField<u8>,

    main_score: SetField<f32>,
    delta_next: SetField<f32>,
    delta_second_next: SetField<f32>,
    rt_seconds: SetField<f32>,
    observed_mobility: SetField<f32>,
    delta_ms1_ms2_mobility: SetField<f32>,
    // ms1_ms2_correlation: SetField<f32>,
    npeaks: SetField<u8>,
    apex_lazyerscore: SetField<f32>,
    apex_lazyerscore_vs_baseline: SetField<f32>,
    apex_norm_lazyerscore_vs_baseline: SetField<f32>,
    ms2_cosine_ref_similarity: SetField<f32>,
    ms2_coelution_score: SetField<f32>,
    ms2_summed_transition_intensity: SetField<f32>,
    ms2_corr_v_gauss: SetField<f32>,
    ms2_lazyerscore: SetField<f32>,
    ms2_isotope_lazyerscore: SetField<f32>,
    ms2_isotope_lazyerscore_ratio: SetField<f32>,

    ms2_mz_errors: SetField<[f32; NUM_MS2_IONS]>,
    ms2_mobility_errors: SetField<[f32; NUM_MS2_IONS]>,

    ms1_cosine_ref_similarity: SetField<f32>,
    ms1_coelution_score: SetField<f32>,
    ms1_summed_precursor_intensity: SetField<f32>,
    ms1_corr_v_gauss: SetField<f32>,

    ms1_mz_errors: SetField<[f32; NUM_MS1_IONS]>,
    ms1_mobility_errors: SetField<[f32; NUM_MS1_IONS]>,

    relative_intensities: SetField<RelativeIntensities>,
    raising_cycles: SetField<u8>,
    falling_cycles: SetField<u8>,
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
        field_name: &'static str,
        // msg: impl ToString,
    ) -> Result<T, DataProcessingError> {
        match self {
            Self::Some(v) => Ok(v),
            Self::None => Err(DataProcessingError::ExpectedSetField {
                field: field_name,
                context: "".into(),
            }),
        }
    }
}

impl<'q> SearchResultBuilder<'q> {
    pub fn with_pre_score(mut self, pre_score: &'q PreScore) -> Self {
        self.digest_slice = SetField::Some(&pre_score.digest);
        self.ref_eg = SetField::Some(&pre_score.query_values.eg);
        self.nqueries = SetField::Some(pre_score.query_values.fragments.num_ions() as u8);
        self.decoy_marking = SetField::Some(pre_score.digest.decoy);
        self.charge = SetField::Some(pre_score.charge);
        self
    }

    pub fn with_sorted_offsets(mut self, offsets: &MzMobilityOffsets) -> Self {
        self.ms1_mz_errors = SetField::Some(offsets.ms1_mz_errors());
        self.ms1_mobility_errors = SetField::Some(offsets.ms1_mobility_errors());
        self.ms2_mz_errors = SetField::Some(offsets.ms2_mz_errors());
        self.ms2_mobility_errors = SetField::Some(offsets.ms2_mobility_errors());

        let mob_errors = offsets.avg_delta_mobs();
        let cum_err = mob_errors.0 + mob_errors.1;
        let obs_mob = (offsets.ref_mobility + cum_err.mean_mobility().unwrap_or(f64::NAN)) as f32;
        let d_err = match (mob_errors.0.mean_mobility(), mob_errors.1.mean_mobility()) {
            (Ok(mz), Ok(mob)) => mz - mob,
            _ => f64::NAN,
        };
        self.delta_ms1_ms2_mobility = SetField::Some(d_err as f32);
        self.observed_mobility = SetField::Some(obs_mob);
        self
    }

    pub fn with_secondary_lazyscores(mut self, lazyscores: SecondaryLazyScores) -> Self {
        self.ms2_lazyerscore = SetField::Some(lazyscores.lazyscore);
        self.ms2_isotope_lazyerscore = SetField::Some(lazyscores.iso_lazyscore);
        self.ms2_isotope_lazyerscore_ratio = SetField::Some(lazyscores.ratio);
        self
    }

    pub fn with_relative_intensities(mut self, relative_intensities: RelativeIntensities) -> Self {
        self.relative_intensities = SetField::Some(relative_intensities);
        self
    }

    pub fn with_main_score(mut self, main_score: MainScore) -> Self {
        // TODO use exhaustive unpacking to make more explicit any values
        // I Might be ignoring.
        let MainScore {
            score,
            delta_next,
            delta_second_next,
            ms2_cosine_ref_sim,
            ms2_coelution_score,
            ms1_coelution_score,
            ms1_cosine_ref_sim,
            lazyscore,
            lazyscore_vs_baseline,
            lazyscore_z,
            ms2_corr_v_gauss,
            ms1_corr_v_gauss,
            npeaks,
            ms1_summed_intensity,
            ms2_summed_intensity,
            retention_time_ms,
            raising_cycles,
            falling_cycles,
            // top_ions: _,
        } = main_score;
        {
            self.main_score = SetField::Some(score);
            self.delta_next = SetField::Some(delta_next);
            self.delta_second_next = SetField::Some(delta_second_next);

            self.rt_seconds = SetField::Some(retention_time_ms as f32 / 1000.0);
            self.ms2_cosine_ref_similarity = SetField::Some(ms2_cosine_ref_sim);
            self.ms2_coelution_score = SetField::Some(ms2_coelution_score);
            self.ms1_coelution_score = SetField::Some(ms1_coelution_score);
            self.ms1_cosine_ref_similarity = SetField::Some(ms1_cosine_ref_sim);
            self.apex_lazyerscore = SetField::Some(lazyscore);
            self.apex_lazyerscore_vs_baseline = SetField::Some(lazyscore_vs_baseline);
            self.apex_norm_lazyerscore_vs_baseline = SetField::Some(lazyscore_z);
            self.npeaks = SetField::Some(npeaks);
            self.ms1_summed_precursor_intensity = SetField::Some(ms1_summed_intensity);
            self.ms2_summed_transition_intensity = SetField::Some(ms2_summed_intensity);
            self.ms2_corr_v_gauss = SetField::Some(ms2_corr_v_gauss);
            self.ms1_corr_v_gauss = SetField::Some(ms1_corr_v_gauss);
            self.raising_cycles = SetField::Some(raising_cycles as u8);
            self.falling_cycles = SetField::Some(falling_cycles as u8);
        }

        self
    }

    pub fn finalize(self) -> Result<IonSearchResults, DataProcessingError> {
        macro_rules! expect_some {
            ($field:ident) => {
                self.$field.expect_some(stringify!($field))?
            };
        }

        let [mz1_e0, mz1_e1, mz1_e2] = expect_some!(ms1_mz_errors);
        let [mz2_e0, mz2_e1, mz2_e2, mz2_e3, mz2_e4, mz2_e5, mz2_e6] = expect_some!(ms2_mz_errors);

        let [mob1_e0, mob1_e1, mob1_e2] = expect_some!(ms1_mobility_errors);
        let [
            mob2_e0,
            mob2_e1,
            mob2_e2,
            mob2_e3,
            mob2_e4,
            mob2_e5,
            mob2_e6,
        ] = expect_some!(ms2_mobility_errors);

        let relints = expect_some!(relative_intensities);
        let [int1_e0, int1_e1, int1_e2] = relints.ms1.get_values();
        let [
            int2_e0,
            int2_e1,
            int2_e2,
            int2_e3,
            int2_e4,
            int2_e5,
            int2_e6,
        ] = relints.ms2.get_values();

        let ref_eg = expect_some!(ref_eg);
        // TODO replace this with exhaustive unpacking.
        let obs_rt_seconds = expect_some!(rt_seconds);
        let delta_theo_rt = obs_rt_seconds - ref_eg.rt_seconds;
        let sq_delta_theo_rt = delta_theo_rt * delta_theo_rt;

        let delta_ms1_ms2_mobility = expect_some!(delta_ms1_ms2_mobility);
        let sq_delta_ms1_ms2_mobility = delta_ms1_ms2_mobility * delta_ms1_ms2_mobility;

        let results = IonSearchResults {
            sequence: String::from(expect_some!(digest_slice).clone()),
            precursor_mz: ref_eg.get_monoisotopic_precursor_mz().unwrap_or(f64::NAN),
            precursor_charge: expect_some!(charge),
            precursor_mobility_query: ref_eg.mobility,
            precursor_rt_query_seconds: ref_eg.rt_seconds,
            nqueries: expect_some!(nqueries),
            is_target: expect_some!(decoy_marking).is_target(),
            main_score: expect_some!(main_score),
            delta_next: expect_some!(delta_next),
            delta_second_next: expect_some!(delta_second_next),
            delta_theo_rt,
            sq_delta_theo_rt,
            obs_rt_seconds,
            obs_mobility: expect_some!(observed_mobility),
            delta_ms1_ms2_mobility,
            sq_delta_ms1_ms2_mobility,
            npeaks: expect_some!(npeaks),
            raising_cycles: expect_some!(raising_cycles),
            falling_cycles: expect_some!(falling_cycles),

            apex_lazyerscore: expect_some!(apex_lazyerscore),
            apex_lazyerscore_vs_baseline: expect_some!(apex_lazyerscore_vs_baseline),
            // ms1_ms2_correlation: self
            //     .ms1_ms2_correlation
            //     .expect_some("ms1_ms2_correlation", "ms1_ms2_correlation")?,
            apex_norm_lazyerscore_vs_baseline: expect_some!(apex_norm_lazyerscore_vs_baseline),
            ms2_cosine_ref_similarity: expect_some!(ms2_cosine_ref_similarity),
            ms2_corr_v_gauss: expect_some!(ms2_corr_v_gauss),
            ms2_summed_transition_intensity: expect_some!(ms2_summed_transition_intensity),
            ms2_lazyerscore: expect_some!(ms2_lazyerscore),
            ms2_isotope_lazyerscore: expect_some!(ms2_isotope_lazyerscore),
            ms2_isotope_lazyerscore_ratio: expect_some!(ms2_isotope_lazyerscore_ratio),

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

            ms2_coelution_score: expect_some!(ms2_coelution_score),
            ms1_cosine_ref_similarity: expect_some!(ms1_cosine_ref_similarity),
            ms1_summed_precursor_intensity: expect_some!(ms1_summed_precursor_intensity),
            ms1_corr_v_gauss: expect_some!(ms1_corr_v_gauss),
            ms1_coelution_score: expect_some!(ms1_coelution_score),
            // ms1_coelution_score: self
            //     .ms1_coelution_score
            //     .expect_some("ms1_coelution_score", "ms1_coelution_score")?,
            ms1_mz_error_0: mz1_e0,
            ms1_mz_error_1: mz1_e1,
            ms1_mz_error_2: mz1_e2,

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
    pub delta_second_next: f32,
    pub obs_rt_seconds: f32,
    obs_mobility: f32,
    delta_theo_rt: f32,
    sq_delta_theo_rt: f32,
    delta_ms1_ms2_mobility: f32,
    // ms1_ms2_correlation: f32,
    sq_delta_ms1_ms2_mobility: f32,
    raising_cycles: u8,
    falling_cycles: u8,

    // MS2
    npeaks: u8,
    apex_lazyerscore: f32,
    apex_lazyerscore_vs_baseline: f32,
    apex_norm_lazyerscore_vs_baseline: f32,
    ms2_cosine_ref_similarity: f32,
    ms2_coelution_score: f32,
    ms2_corr_v_gauss: f32,
    ms2_summed_transition_intensity: f32,
    ms2_lazyerscore: f32,
    ms2_isotope_lazyerscore: f32,
    ms2_isotope_lazyerscore_ratio: f32,

    // MS2 - Split
    // Flattening manually bc serde(flatten)
    // is not supported by csv ...
    // https://github.com/BurntSushi/rust-csv/pull/223
    // Q: Is it supported by parquet?
    // A: As of 2025-May-20, it is not.
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
    ms1_corr_v_gauss: f32,

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

pub struct ResultParquetWriter {
    row_group_size: usize,
    writer: SerializedFileWriter<File>,
    buffer: Vec<IonSearchResults>,
}

impl ResultParquetWriter {
    pub fn new(out_path: impl AsRef<Path>, row_group_size: usize) -> Result<Self, std::io::Error> {
        let file = match File::create_new(out_path.as_ref()) {
            Ok(file) => file,
            Err(err) => {
                tracing::error!(
                    "Failed to open file {:?} with error: {}",
                    out_path.as_ref(),
                    err
                );
                return Err(err);
            }
        };
        let results: &[IonSearchResults] = &[];
        let schema = results.schema().unwrap();
        let writer = SerializedFileWriter::new(file, schema, Default::default()).unwrap();
        Ok(Self {
            buffer: Vec::with_capacity(row_group_size),
            writer,
            row_group_size,
        })
    }

    fn flush_to_file(&mut self) {
        debug!("Flushing {} results to file", self.buffer.len());
        let mut row_group = self.writer.next_row_group().unwrap();
        self.buffer
            .as_slice()
            .write_to_row_group(&mut row_group)
            .unwrap();
        row_group.close().unwrap();
        self.buffer.clear();
    }

    pub fn add(&mut self, result: IonSearchResults) {
        self.buffer.push(result);
        if self.buffer.len() >= self.row_group_size {
            self.flush_to_file();
        }
    }

    pub fn close(mut self) {
        // TODO: add some logging ...
        if !self.buffer.is_empty() {
            self.flush_to_file();
        }
        self.writer.close().unwrap();
    }
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
