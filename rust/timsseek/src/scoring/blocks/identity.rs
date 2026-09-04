//! Identity family -- peptide + precursor metadata. Hand-written because it is
//! irreducibly mixed-dtype (`Peptide`, `bool`, `f64`, `u32`, `u8`, `f32`); you
//! never add a *score* here.

use std::sync::Arc;
use timsquery::models::{
    GroupCode,
    OwnedSourceId,
    RowIdx,
};

use crate::models::DecoyMarking;
use crate::models::sequence::Peptide;
use crate::scoring::apex_finding::PeptideMetadata;
use crate::scoring::blocks::{
    ColSink,
    NameSink,
    SchemaSink,
    ScoreBlock,
};

#[derive(Debug, Clone, serde::Serialize)]
pub struct Identity {
    pub peptide: Peptide,
    /// The arena row this result came from, and the group it competes in. Both
    /// opaque and unserializable: the row orders the q-value tie-break without
    /// a caller-supplied id reaching the sort, and the group is only ever
    /// compared. The ids a reader wants are resolved from the arena at the
    /// writer -- see [`crate::scoring::parquet_writer`].
    #[serde(skip)]
    pub row: RowIdx,
    #[serde(skip)]
    pub group: GroupCode,
    #[serde(skip)]
    pub source_id: OwnedSourceId,
    pub precursor_mz: f64,
    pub precursor_charge: u8,
    pub precursor_mobility: f32,
    pub is_target: bool,
}

impl Identity {
    /// Hand-written twin of what `#[derive(ScoreBlock)]` emits: no identity
    /// field is a linear-lane feature.
    pub const LINEAR_LEN: usize = 0;
    /// `precursor_mz_round5`, `precursor_charge`, `precursor_mobility` -- the
    /// context features. `is_target` is Parquet-only, as are the ids the writer
    /// resolves from `row`.
    pub const NONLINEAR_LEN: usize = 3;

    pub fn compute(metadata: &PeptideMetadata) -> Self {
        Self {
            peptide: metadata.digest.clone(),
            row: metadata.handles.row,
            group: metadata.handles.group,
            source_id: metadata.source_id.clone(),
            precursor_mz: metadata.ref_precursor_mz,
            precursor_charge: metadata.charge,
            precursor_mobility: metadata.ref_mobility_ook0,
            is_target: metadata.digest.decoy.is_target(),
        }
    }

    /// Whether two results compete: same decoy group, same charge.
    ///
    /// Defined in terms of [`Self::competition_key`], so sorting by that key
    /// always leaves competing results adjacent -- the predicate and the sort
    /// order cannot drift apart.
    pub fn competes_with(&self, other: &Self) -> bool {
        self.competition_key() == other.competition_key()
    }

    pub fn competition_key(&self) -> (GroupCode, u8) {
        (self.group, self.precursor_charge)
    }

    /// The observed mobility is a sentinel on non-scoreable axes; drop the
    /// precursor mobility to NaN.
    pub fn neutralize_mobility(&mut self) {
        self.precursor_mobility = f32::NAN;
    }

    pub fn linear_feature_array(&self) -> [f64; Self::LINEAR_LEN] {
        []
    }

    /// Values for [`Identity::nonlinear_feature_names`], same order.
    pub fn nonlinear_feature_array(&self) -> [f64; Self::NONLINEAR_LEN] {
        [
            (self.precursor_mz / 5.0).round(),
            self.precursor_charge as f64,
            self.precursor_mobility as f64,
        ]
    }

    pub fn sample_default() -> Self {
        Self {
            peptide: Peptide {
                raw: Arc::from("PEPTIDEK"),
                decoy: DecoyMarking::Target,
                sequence_features: false,
            },
            row: RowIdx::default(),
            group: GroupCode::default(),
            source_id: OwnedSourceId::placeholder(),
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
        o.f64("precursor_mz", self.precursor_mz);
        o.u8("precursor_charge", self.precursor_charge);
        o.f32("precursor_mobility", self.precursor_mobility);
        o.bool("is_target", self.is_target);
    }

    fn column_schema(o: &mut SchemaSink) {
        o.str("sequence");
        o.f64("precursor_mz");
        o.u8("precursor_charge");
        o.f32("precursor_mobility");
        o.bool("is_target");
    }

    fn nonlinear_feature_names(o: &mut NameSink) {
        o.push("precursor_mz_round5");
        o.push("precursor_charge");
        o.push("precursor_mobility");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The two halves of the nonlinear lane -- the per-record value array and
    /// the set-level name walk -- must agree on count, order, and which value
    /// sits under which name. Hand-written here (the derive can't check this
    /// block), so this is the only place the pairing is asserted.
    #[test]
    fn nonlinear_lane_names_match_values() {
        let identity = Identity::sample_default();
        let values = identity.nonlinear_feature_array();

        let mut names = NameSink::new();
        Identity::nonlinear_feature_names(&mut names);
        let names = names.into_names();
        assert_eq!(names.len(), Identity::NONLINEAR_LEN);

        let expected = [
            ("precursor_mz_round5", (identity.precursor_mz / 5.0).round()),
            ("precursor_charge", identity.precursor_charge as f64),
            ("precursor_mobility", identity.precursor_mobility as f64),
        ];
        assert_eq!(expected.len(), names.len());
        for (j, (name, value)) in expected.iter().enumerate() {
            assert_eq!(&*names[j], *name);
            assert_eq!(values[j], *value);
        }
    }

    #[test]
    fn column_schema_matches_columns() {
        let identity = Identity::sample_default();

        let mut cols = ColSink::new();
        identity.columns(&mut cols);
        let (col_fields, _) = cols.finish();

        let mut schema = SchemaSink::new();
        Identity::column_schema(&mut schema);
        let schema_fields = schema.into_fields();

        assert_eq!(col_fields.len(), schema_fields.len());
        for (a, b) in col_fields.iter().zip(schema_fields.iter()) {
            assert_eq!(a.name(), b.name());
            assert_eq!(a.data_type(), b.data_type());
            assert_eq!(a.is_nullable(), b.is_nullable());
        }
    }
}
