use serde::Serialize;
// The derive and the trait share a name but live in different namespaces
// (macro vs type), so both are in scope here unambiguously.
use timsseek_macros::ScoreBlock;

use super::apex_finding::{
    ApexBlocks,
    PeptideMetadata,
    RelativeIntensityCollector,
};
use super::blocks::apex_evidence::ApexEvidence;
use super::blocks::apex_features::ApexFeatures;
use super::blocks::counts::{
    ApexCounts,
    FinalizeCounts,
};
use super::blocks::identity::Identity;
use super::blocks::intensities::Intensities;
use super::blocks::ion_errors::IonErrors;
use super::blocks::lazy::{
    ApexLazyScores,
    SecondaryLazyScores,
};
use super::blocks::mobility::Mobility;
use super::blocks::primary::PrimaryScores;
use super::blocks::relative_intensities::RelativeIntensities;
use super::blocks::result_meta::ResultMeta;
use super::blocks::rt::Rt;
use super::blocks::{
    SchemaSink,
    ScoreBlock,
};
use super::offsets::MzMobilityOffsets;
use super::pipeline::SecondaryLazyScoresRaw;
use crate::models::sequence::Peptide;

/// Inputs for the finalize-stage assembly of [`ScoringFields`]. Constructing
/// this struct IS the completeness guarantee: a score that needs data not yet
/// present forces one new field here plus one line at the single construction
/// site (`pipeline::finalize_results`).
pub struct FinalizeInputs<'a> {
    pub metadata: &'a PeptideMetadata,
    pub offsets: &'a MzMobilityOffsets,
    pub rel_inten: RelativeIntensityCollector,
    pub secondary_lazy: SecondaryLazyScoresRaw,
    pub nqueries: u8,
    /// Apex-stage blocks, moved in at zero cost.
    pub apex: ApexBlocks,
}

/// Shared scoring fields produced by Phase 3, as a composition of typed
/// blocks. Each block owns its parquet/ML projections (`columns`, the lane
/// methods) in one file under [`super::blocks`]; the finalize-stage blocks own
/// their `compute` there too, while the apex-stage blocks are built in
/// [`super::apex_finding`].
///
/// The `#[block]` fields make this a DELEGATING `#[derive(ScoreBlock)]`: the
/// same derive that projects a leaf block's `#[feat(...)]` scalars walks this
/// list to emit the composed struct's Parquet columns/schema, its two lane
/// value arrays and their name walks, and `sample_default`. Field order is
/// load-bearing -- parquet column order and the positional ML value vector both
/// follow it -- so folding every walk out of the one declaration order is what
/// makes them *impossible* to desync.
///
/// The lane widths compose the same way the values do: `LINEAR_LEN` is the sum
/// of the blocks' own inherent `LINEAR_LEN` consts, so the composed lane array
/// is `[f64; N]` with `N` known at compile time and no runtime length
/// bookkeeping anywhere below it.
///
/// [`Self::compute`] (per-block dep signatures vary) and
/// [`Self::neutralize_mobility`] (only the mobility-derived blocks
/// participate) stay hand-written below: adding a block is a two-place edit
/// (this list + `compute`), not five.
#[derive(Debug, Clone, Serialize, ScoreBlock)]
pub struct ScoringFields {
    #[block]
    pub identity: Identity,
    #[block]
    pub rt: Rt,
    #[block]
    pub mobility: Mobility,
    #[block]
    pub primary: PrimaryScores,
    #[block]
    pub evidence: ApexEvidence,
    #[block]
    pub features: ApexFeatures,
    #[block]
    pub apex_lazy: ApexLazyScores,
    #[block]
    pub secondary_lazy: SecondaryLazyScores,
    #[block]
    pub counts: ApexCounts,
    #[block]
    pub finalize_counts: FinalizeCounts,
    #[block]
    pub intensities: Intensities,
    #[block]
    pub ion_errors: IonErrors,
    #[block]
    pub rel_intensities: RelativeIntensities,
}

impl ScoringFields {
    /// Assemble every finalize-stage block, moving the apex blocks in.
    pub fn compute(inp: FinalizeInputs) -> Self {
        let obs_rt_seconds = inp.apex.retention_time_ms as f32 / 1000.0;
        Self {
            identity: Identity::compute(inp.metadata),
            rt: Rt::compute(inp.metadata, obs_rt_seconds),
            mobility: Mobility::compute(inp.offsets),
            primary: inp.apex.primary,
            evidence: inp.apex.evidence,
            features: inp.apex.features,
            apex_lazy: inp.apex.apex_lazy,
            secondary_lazy: SecondaryLazyScores::from(inp.secondary_lazy),
            counts: inp.apex.counts,
            finalize_counts: FinalizeCounts {
                n_scored_fragments: inp.nqueries,
            },
            intensities: inp.apex.intensities,
            ion_errors: IonErrors::compute(inp.offsets),
            rel_intensities: RelativeIntensities::from_collector(&inp.rel_inten),
        }
    }

    /// Zero out (as NaN) every mobility-derived field. Used when the run's
    /// mobility axis is not a searchable TIMS 1/K0 (mzML/FAIMS). Each block
    /// NaNs its own mobility-derived fields (including derived squares), so a
    /// future move of a field cannot desync from its source.
    pub fn neutralize_mobility(&mut self) {
        self.identity.neutralize_mobility();
        self.mobility.neutralize();
        self.ion_errors.neutralize();
    }

    /// Baseline test fixture with every field populated. Callers (including
    /// other crates' tests) tweak the identity/score fields they care about.
    pub fn sample(peptide: Peptide) -> Self {
        let mut s = Self::sample_default();
        s.identity.peptide = peptide;
        s
    }
}

/// Phase 3 output. All scoring fields guaranteed populated.
#[derive(Debug, Clone, Serialize)]
pub struct ScoredCandidate {
    pub scoring: ScoringFields,
}

/// After target-decoy competition.
#[derive(Debug, Clone, Serialize)]
pub struct CompetedCandidate {
    pub scoring: ScoringFields,
    pub delta_group_ln1p_diff: f32,
    pub delta_group_ln1p_ratio: f32,
    /// Scratch field for CrossValidatedScorer (written during rescore)
    pub(crate) discriminant_score: f32,
    /// Scratch field for q-value assignment
    pub(crate) qvalue: f32,
}

impl CompetedCandidate {
    /// The post-model meta block (used for the ML delta-group features).
    pub(crate) fn result_meta(&self) -> ResultMeta {
        ResultMeta {
            delta_group_ln1p_diff: self.delta_group_ln1p_diff,
            delta_group_ln1p_ratio: self.delta_group_ln1p_ratio,
            discriminant_score: self.discriminant_score,
            qvalue: self.qvalue,
        }
    }
}

/// After rescoring. Written to Parquet.
#[derive(Debug, Clone, Serialize)]
pub struct FinalResult {
    pub scoring: ScoringFields,
    pub delta_group_ln1p_diff: f32,
    pub delta_group_ln1p_ratio: f32,
    pub discriminant_score: f32,
    pub qvalue: f32,
}

impl FinalResult {
    /// Schema/test fixture with zeroed meta fields.
    #[cfg(test)]
    pub fn sample() -> Self {
        Self {
            scoring: ScoringFields::sample_default(),
            delta_group_ln1p_diff: 0.0,
            delta_group_ln1p_ratio: 0.0,
            discriminant_score: 0.0,
            qvalue: 0.0,
        }
    }

    /// The post-model meta block (used for the Parquet meta columns).
    pub(crate) fn result_meta(&self) -> ResultMeta {
        ResultMeta {
            delta_group_ln1p_diff: self.delta_group_ln1p_diff,
            delta_group_ln1p_ratio: self.delta_group_ln1p_ratio,
            discriminant_score: self.discriminant_score,
            qvalue: self.qvalue,
        }
    }

    /// Value-free Parquet schema: scoring blocks (composition order), the
    /// post-model meta block, then the ids -- the same three blocks `emit_row`
    /// writes, in the same order.
    ///
    /// The ID columns have their own block because the writer resolves them
    /// from the arena, independently of scoring metadata. See `parquet_writer::Ids`.
    pub fn column_schema(o: &mut SchemaSink) {
        <ScoringFields as ScoreBlock>::column_schema(o);
        <ResultMeta as ScoreBlock>::column_schema(o);
        <crate::scoring::parquet_writer::Ids<'_> as ScoreBlock>::column_schema(o);
    }
}

// ---------------------------------------------------------------------------
// Stage conversions
// ---------------------------------------------------------------------------

impl ScoredCandidate {
    /// Convert into a `CompetedCandidate` with the given log-space delta values.
    ///
    /// Items that are alone in their group (no competitor) should pass
    /// `f32::NAN` for both deltas.
    pub fn into_competed(
        self,
        delta_group_ln1p_diff: f32,
        delta_group_ln1p_ratio: f32,
    ) -> CompetedCandidate {
        CompetedCandidate {
            scoring: self.scoring,
            delta_group_ln1p_diff,
            delta_group_ln1p_ratio,
            discriminant_score: f32::NAN,
            qvalue: f32::NAN,
        }
    }
}

impl CompetedCandidate {
    /// Promote to a `FinalResult` (all fields frozen).
    pub fn into_final(self) -> FinalResult {
        FinalResult {
            scoring: self.scoring,
            delta_group_ln1p_diff: self.delta_group_ln1p_diff,
            delta_group_ln1p_ratio: self.delta_group_ln1p_ratio,
            discriminant_score: self.discriminant_score,
            qvalue: self.qvalue,
        }
    }
}
