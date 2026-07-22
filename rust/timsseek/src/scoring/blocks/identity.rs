//! Identity family — peptide + precursor metadata. Hand-written because it is
//! irreducibly mixed-dtype (`Peptide`, `bool`, `f64`, `u32`, `u8`, `f32`); you
//! never add a *score* here.

use std::sync::Arc;

use crate::models::DecoyMarking;
use crate::models::sequence::Peptide;
use crate::scoring::apex_finding::PeptideMetadata;
use crate::scoring::blocks::{
    ColSink,
    FeatSink,
    NameSink,
    ScoreBlock,
};

#[derive(Debug, Clone, serde::Serialize)]
pub struct Identity {
    pub peptide: Peptide,
    pub library_id: u32,
    pub decoy_group_id: u32,
    pub precursor_mz: f64,
    pub precursor_charge: u8,
    pub precursor_mobility: f32,
    pub is_target: bool,
}

impl Identity {
    pub fn compute(metadata: &PeptideMetadata) -> Self {
        Self {
            peptide: metadata.digest.clone(),
            library_id: metadata.library_id,
            decoy_group_id: metadata.digest.decoy_group,
            precursor_mz: metadata.ref_precursor_mz,
            precursor_charge: metadata.charge,
            precursor_mobility: metadata.ref_mobility_ook0,
            is_target: metadata.digest.decoy.is_target(),
        }
    }

    /// The observed mobility is a sentinel on non-scoreable axes; drop the
    /// precursor mobility to NaN.
    pub fn neutralize_mobility(&mut self) {
        self.precursor_mobility = f32::NAN;
    }

    pub fn sample() -> Self {
        Self {
            peptide: Peptide {
                raw: Arc::from("PEPTIDEK"),
                parsed: None,
                decoy: DecoyMarking::Target,
                decoy_group: 0,
            },
            library_id: 1,
            decoy_group_id: 0,
            precursor_mz: 500.0,
            precursor_charge: 2,
            precursor_mobility: 0.9,
            is_target: true,
        }
    }
}

impl ScoreBlock for Identity {
    fn columns(&self, o: &mut ColSink) {
        o.str("sequence", self.peptide.as_str());
        o.u32("library_id", self.library_id);
        o.u32("decoy_group_id", self.decoy_group_id);
        o.f64("precursor_mz", self.precursor_mz);
        o.u8("precursor_charge", self.precursor_charge);
        o.f32("precursor_mobility", self.precursor_mobility);
        o.bool("is_target", self.is_target);
    }

    fn features(&self, o: &mut FeatSink) {
        // library_id / decoy_group_id / is_target are Parquet-only (not features).
        o.push((self.precursor_mz / 5.0).round());
        o.push(self.precursor_charge as f64);
        o.push(self.precursor_mobility as f64);
    }

    fn feature_names(o: &mut NameSink) {
        o.push("precursor_mz_round5");
        o.push("precursor_charge");
        o.push("precursor_mobility");
    }
}
