use super::apex_finding::{
    ApexScore,
    CandidateContext,
    RelativeIntensities,
};
use super::offsets::MzMobilityOffsets;
use super::pipeline::SecondaryLazyScores;
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
use timsquery::TimsElutionGroup;
use tracing::debug;

#[derive(Debug, Default)]
pub struct SearchResultBuilder<'q> {
    digest_slice: SetField<&'q DigestSlice>,
    ref_eg: SetField<&'q TimsElutionGroup<IonAnnot>>,
    decoy_marking: SetField<DecoyMarking>,
    library_id: SetField<u32>,
    decoy_group_id: SetField<u32>,
    charge: SetField<u8>,
    nqueries: SetField<u8>,

    // Reference values for new metadata-based approach
    ref_precursor_mz: SetField<f64>,
    ref_rt_seconds: SetField<f32>,
    ref_mobility: SetField<f32>,

    main_score: SetField<f32>,
    delta_next: SetField<f32>,
    delta_second_next: SetField<f32>,
    rt_seconds: SetField<f32>,
    observed_mobility: SetField<f32>,
    delta_ms1_ms2_mobility: SetField<f32>,
    // ms1_ms2_correlation: SetField<f32>,
    npeaks: SetField<u8>,
    apex_lazyerscore: SetField<f32>,
    ms2_summed_transition_intensity: SetField<f32>,
    ms2_lazyerscore: SetField<f32>,
    ms2_isotope_lazyerscore: SetField<f32>,
    ms2_isotope_lazyerscore_ratio: SetField<f32>,
    lazyscore_z: SetField<f32>,
    lazyscore_vs_baseline: SetField<f32>,

    ms2_mz_errors: SetField<[f32; NUM_MS2_IONS]>,
    ms2_mobility_errors: SetField<[f32; NUM_MS2_IONS]>,

    ms1_summed_precursor_intensity: SetField<f32>,

    ms1_mz_errors: SetField<[f32; NUM_MS1_IONS]>,
    ms1_mobility_errors: SetField<[f32; NUM_MS1_IONS]>,

    // Split product & apex features
    split_product_score: SetField<f32>,
    cosine_au_score: SetField<f32>,
    scribe_au_score: SetField<f32>,
    coelution_gradient_cosine: SetField<f32>,
    coelution_gradient_scribe: SetField<f32>,
    cosine_weighted_coelution: SetField<f32>,
    cosine_gradient_consistency: SetField<f32>,
    scribe_weighted_coelution: SetField<f32>,
    scribe_gradient_consistency: SetField<f32>,
    peak_shape: SetField<f32>,
    ratio_cv: SetField<f32>,
    centered_apex: SetField<f32>,
    precursor_coelution: SetField<f32>,
    fragment_coverage: SetField<f32>,
    precursor_apex_match: SetField<f32>,
    xic_quality: SetField<f32>,
    fragment_apex_agreement: SetField<f32>,
    isotope_correlation: SetField<f32>,
    gaussian_correlation: SetField<f32>,
    per_frag_gaussian_corr: SetField<f32>,

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
    pub fn with_candidate_context(
        mut self,
        candidate_context: &'q CandidateContext<IonAnnot, DigestSlice>,
    ) -> Self {
        self.library_id = SetField::Some(candidate_context.query_values.eg.id() as u32);
        self.digest_slice = SetField::Some(&candidate_context.label);
        self.ref_eg = SetField::Some(&candidate_context.query_values.eg);
        self.nqueries = SetField::Some(candidate_context.query_values.fragments.num_ions() as u8);
        self.decoy_marking = SetField::Some(candidate_context.label.decoy);
        self.charge = SetField::Some(candidate_context.charge);
        self.decoy_group_id = SetField::Some(candidate_context.label.decoy_group);
        self
    }

    pub fn with_metadata(mut self, metadata: &'q super::apex_finding::PeptideMetadata) -> Self {
        self.library_id = SetField::Some(metadata.library_id);
        self.digest_slice = SetField::Some(&metadata.digest);
        self.decoy_marking = SetField::Some(metadata.digest.decoy);
        self.charge = SetField::Some(metadata.charge);
        self.decoy_group_id = SetField::Some(metadata.digest.decoy_group);

        // Store ref values for later use in finalize()
        self.ref_precursor_mz = SetField::Some(metadata.ref_precursor_mz);
        self.ref_rt_seconds = SetField::Some(metadata.ref_rt_seconds);
        self.ref_mobility = SetField::Some(metadata.ref_mobility_ook0);

        self
    }

    pub fn with_nqueries(mut self, nqueries: u8) -> Self {
        self.nqueries = SetField::Some(nqueries);
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

    pub fn with_apex_score(mut self, main_score: &ApexScore) -> Self {
        let ApexScore {
            score,
            retention_time_ms,
            joint_apex_cycle: _,
            split_product,
            features,
            delta_next,
            delta_second_next,
            lazyscore,
            lazyscore_vs_baseline,
            lazyscore_z,
            npeaks,
            ms2_summed_intensity,
            ms1_summed_intensity,
            raising_cycles,
            falling_cycles,
        } = *main_score;
        {
            self.main_score = SetField::Some(score);
            self.delta_next = SetField::Some(delta_next);
            self.delta_second_next = SetField::Some(delta_second_next);
            self.rt_seconds = SetField::Some(retention_time_ms as f32 / 1000.0);

            self.split_product_score = SetField::Some(split_product.base_score);
            self.cosine_au_score = SetField::Some(split_product.cosine_au);
            self.scribe_au_score = SetField::Some(split_product.scribe_au);
            self.coelution_gradient_cosine = SetField::Some(split_product.cosine_cg);
            self.coelution_gradient_scribe = SetField::Some(split_product.scribe_cg);
            self.cosine_weighted_coelution = SetField::Some(split_product.cosine_weighted_coelution);
            self.cosine_gradient_consistency = SetField::Some(split_product.cosine_gradient_consistency);
            self.scribe_weighted_coelution = SetField::Some(split_product.scribe_weighted_coelution);
            self.scribe_gradient_consistency = SetField::Some(split_product.scribe_gradient_consistency);

            self.peak_shape = SetField::Some(features.peak_shape);
            self.ratio_cv = SetField::Some(features.ratio_cv);
            self.centered_apex = SetField::Some(features.centered_apex);
            self.precursor_coelution = SetField::Some(features.precursor_coelution);
            self.fragment_coverage = SetField::Some(features.fragment_coverage);
            self.precursor_apex_match = SetField::Some(features.precursor_apex_match);
            self.xic_quality = SetField::Some(features.xic_quality);
            self.fragment_apex_agreement = SetField::Some(features.fragment_apex_agreement);
            self.isotope_correlation = SetField::Some(features.isotope_correlation);
            self.gaussian_correlation = SetField::Some(features.gaussian_correlation);
            self.per_frag_gaussian_corr = SetField::Some(features.per_frag_gaussian_corr);

            self.apex_lazyerscore = SetField::Some(lazyscore);
            self.lazyscore_z = SetField::Some(lazyscore_z);
            self.lazyscore_vs_baseline = SetField::Some(lazyscore_vs_baseline);
            self.npeaks = SetField::Some(npeaks);
            self.ms1_summed_precursor_intensity = SetField::Some(ms1_summed_intensity);
            self.ms2_summed_transition_intensity = SetField::Some(ms2_summed_intensity);
            self.raising_cycles = SetField::Some(raising_cycles);
            self.falling_cycles = SetField::Some(falling_cycles);
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

        // Use stored ref values or fallback to ref_eg if available
        let ref_mz = if self.ref_precursor_mz.is_some() {
            expect_some!(ref_precursor_mz)
        } else {
            let ref_eg = expect_some!(ref_eg);
            ref_eg.mono_precursor_mz()
        };

        let ref_rt = if self.ref_rt_seconds.is_some() {
            expect_some!(ref_rt_seconds)
        } else {
            let ref_eg = expect_some!(ref_eg);
            ref_eg.rt_seconds()
        };

        let ref_mob = if self.ref_mobility.is_some() {
            expect_some!(ref_mobility)
        } else {
            let ref_eg = expect_some!(ref_eg);
            ref_eg.mobility_ook0()
        };

        let obs_rt_seconds = expect_some!(rt_seconds);
        let delta_theo_rt = obs_rt_seconds - ref_rt;
        let sq_delta_theo_rt = delta_theo_rt * delta_theo_rt;

        let delta_ms1_ms2_mobility = expect_some!(delta_ms1_ms2_mobility);
        let sq_delta_ms1_ms2_mobility = delta_ms1_ms2_mobility * delta_ms1_ms2_mobility;

        let results = IonSearchResults {
            sequence: String::from(expect_some!(digest_slice).clone()),
            library_id: expect_some!(library_id),
            decoy_group_id: expect_some!(decoy_group_id),
            precursor_mz: ref_mz,
            precursor_charge: expect_some!(charge),
            precursor_mobility_query: ref_mob,
            precursor_rt_query_seconds: ref_rt,
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
            ms2_summed_transition_intensity: expect_some!(ms2_summed_transition_intensity),
            ms2_lazyerscore: expect_some!(ms2_lazyerscore),
            ms2_isotope_lazyerscore: expect_some!(ms2_isotope_lazyerscore),
            ms2_isotope_lazyerscore_ratio: expect_some!(ms2_isotope_lazyerscore_ratio),
            lazyscore_z: expect_some!(lazyscore_z),
            lazyscore_vs_baseline: expect_some!(lazyscore_vs_baseline),

            split_product_score: expect_some!(split_product_score),
            cosine_au_score: expect_some!(cosine_au_score),
            scribe_au_score: expect_some!(scribe_au_score),
            coelution_gradient_cosine: expect_some!(coelution_gradient_cosine),
            coelution_gradient_scribe: expect_some!(coelution_gradient_scribe),
            cosine_weighted_coelution: expect_some!(cosine_weighted_coelution),
            cosine_gradient_consistency: expect_some!(cosine_gradient_consistency),
            scribe_weighted_coelution: expect_some!(scribe_weighted_coelution),
            scribe_gradient_consistency: expect_some!(scribe_gradient_consistency),
            peak_shape: expect_some!(peak_shape),
            ratio_cv: expect_some!(ratio_cv),
            centered_apex: expect_some!(centered_apex),
            precursor_coelution: expect_some!(precursor_coelution),
            fragment_coverage: expect_some!(fragment_coverage),
            precursor_apex_match: expect_some!(precursor_apex_match),
            xic_quality: expect_some!(xic_quality),
            fragment_apex_agreement: expect_some!(fragment_apex_agreement),
            isotope_correlation: expect_some!(isotope_correlation),
            gaussian_correlation: expect_some!(gaussian_correlation),
            per_frag_gaussian_corr: expect_some!(per_frag_gaussian_corr),

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

            ms1_summed_precursor_intensity: expect_some!(ms1_summed_precursor_intensity),
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

            discriminant_score: f32::NAN,
            qvalue: f32::NAN,
            delta_group: f32::NAN,
            delta_group_ratio: f32::NAN,
            recalibrated_query_rt: ref_rt,
            calibrated_sq_delta_theo_rt: sq_delta_theo_rt,
        };

        Ok(results)
    }
}

/// Contains the results as they will be serialized to
/// parquet.
///
/// Napkin math ... as of Sept 23/2025 this struct is 265 bytes
/// Eyeballing a human proteome without mods is 1.8M peptides
/// So ... 500MB-ish / proteome in memory.
#[derive(Debug, Clone, Serialize, ParquetRecordWriter)]
pub struct IonSearchResults {
    pub sequence: String,
    pub library_id: u32,
    pub decoy_group_id: u32,
    pub precursor_mz: f64,
    pub precursor_charge: u8,
    pub precursor_mobility_query: f32,
    pub precursor_rt_query_seconds: f32,
    pub recalibrated_query_rt: f32,
    pub nqueries: u8,
    pub is_target: bool,

    // Combined
    pub main_score: f32,
    pub delta_next: f32,
    pub delta_second_next: f32,
    pub obs_rt_seconds: f32,
    pub obs_mobility: f32,
    pub delta_theo_rt: f32,
    pub sq_delta_theo_rt: f32,
    pub calibrated_sq_delta_theo_rt: f32,
    pub delta_ms1_ms2_mobility: f32,
    // ms1_ms2_correlation: f32,
    pub sq_delta_ms1_ms2_mobility: f32,
    pub raising_cycles: u8,
    pub falling_cycles: u8,
    pub delta_group: f32,
    pub delta_group_ratio: f32,

    // MS2
    pub npeaks: u8,
    pub apex_lazyerscore: f32,
    pub ms2_summed_transition_intensity: f32,
    pub ms2_lazyerscore: f32,
    pub ms2_isotope_lazyerscore: f32,
    pub ms2_isotope_lazyerscore_ratio: f32,
    pub lazyscore_z: f32,
    pub lazyscore_vs_baseline: f32,

    // Split product & apex features
    pub split_product_score: f32,
    pub cosine_au_score: f32,
    pub scribe_au_score: f32,
    pub coelution_gradient_cosine: f32,
    pub coelution_gradient_scribe: f32,
    pub cosine_weighted_coelution: f32,
    pub cosine_gradient_consistency: f32,
    pub scribe_weighted_coelution: f32,
    pub scribe_gradient_consistency: f32,
    pub peak_shape: f32,
    pub ratio_cv: f32,
    pub centered_apex: f32,
    pub precursor_coelution: f32,
    pub fragment_coverage: f32,
    pub precursor_apex_match: f32,
    pub xic_quality: f32,
    pub fragment_apex_agreement: f32,
    pub isotope_correlation: f32,
    pub gaussian_correlation: f32,
    pub per_frag_gaussian_corr: f32,

    // MS2 - Split
    pub ms2_mz_error_0: f32,
    pub ms2_mz_error_1: f32,
    pub ms2_mz_error_2: f32,
    pub ms2_mz_error_3: f32,
    pub ms2_mz_error_4: f32,
    pub ms2_mz_error_5: f32,
    pub ms2_mz_error_6: f32,
    pub ms2_mobility_error_0: f32,
    pub ms2_mobility_error_1: f32,
    pub ms2_mobility_error_2: f32,
    pub ms2_mobility_error_3: f32,
    pub ms2_mobility_error_4: f32,
    pub ms2_mobility_error_5: f32,
    pub ms2_mobility_error_6: f32,

    // MS1
    pub ms1_summed_precursor_intensity: f32,

    // MS1 Split
    pub ms1_mz_error_0: f32,
    pub ms1_mz_error_1: f32,
    pub ms1_mz_error_2: f32,
    pub ms1_mobility_error_0: f32,
    pub ms1_mobility_error_1: f32,
    pub ms1_mobility_error_2: f32,

    // Relative Intensities
    pub ms1_inten_ratio_0: f32,
    pub ms1_inten_ratio_1: f32,
    pub ms1_inten_ratio_2: f32,

    pub ms2_inten_ratio_0: f32,
    pub ms2_inten_ratio_1: f32,
    pub ms2_inten_ratio_2: f32,
    pub ms2_inten_ratio_3: f32,
    pub ms2_inten_ratio_4: f32,
    pub ms2_inten_ratio_5: f32,
    pub ms2_inten_ratio_6: f32,

    pub discriminant_score: f32,
    pub qvalue: f32,
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
