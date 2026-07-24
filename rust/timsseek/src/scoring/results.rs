use serde::Serialize;

use super::apex_finding::{
    ApexBlocks,
    PeptideMetadata,
    RelativeIntensityCollector,
};
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
use super::blocks::split_product::SplitProduct;
use super::blocks::{
    ColSink,
    FeatSink,
    NameSink,
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

/// Compose [`ScoringFields`] from an ordered block list, deriving the struct
/// and the four purely-mechanical walks (`push_columns`, `push_features`,
/// `push_feature_names`, `sample_default`) from that one list. Order is
/// load-bearing (parquet columns and the positional ML value vector both
/// follow it), so folding all four from a single ordered source is what makes
/// their order *impossible* to desync.
///
/// `compute` (per-block dep signatures vary) and `neutralize_mobility` (only
/// the mobility-derived blocks participate) stay hand-written below: adding a
/// block is a two-place edit (this list + `compute`), not five.
macro_rules! compose_scoring_fields {
    (
        $(#[doc = $doc:literal])*
        pub struct $Name:ident {
            $( pub $fname:ident : $fty:ty ),* $(,)?
        }
    ) => {
        $(#[doc = $doc])*
        #[derive(Debug, Clone, Serialize)]
        pub struct $Name {
            $( pub $fname : $fty ),*
        }

        impl $Name {
            /// Emit every block's Parquet columns, in the composition order.
            pub fn push_columns(&self, o: &mut ColSink) {
                $( self.$fname.columns(o); )*
            }

            /// Emit every block's Parquet SCHEMA (value-free), in composition order —
            /// same order as `push_columns`.
            pub fn push_column_schema(o: &mut $crate::scoring::blocks::SchemaSink) {
                $( <$fty as $crate::scoring::blocks::ScoreBlock>::column_schema(o); )*
            }

            /// Emit every block's direct ML feature *values* (not the
            /// cross-field derived ones, nor the conditional sequence block).
            pub fn push_features(&self, o: &mut FeatSink) {
                $( self.$fname.features(o); )*
            }

            /// Emit every block's ML feature *names*, in the same order as
            /// [`Self::push_features`] (set-level; no `&self`).
            pub fn push_feature_names(o: &mut NameSink) {
                $( <$fty as ScoreBlock>::feature_names(o); )*
            }

            /// Fixture using the placeholder peptide from [`Identity::sample`].
            pub fn sample_default() -> Self {
                Self {
                    $( $fname : <$fty>::sample() ),*
                }
            }
        }
    };
}

compose_scoring_fields! {
    /// Shared scoring fields produced by Phase 3, as a composition of typed
    /// blocks. Each block owns its parquet/ML projections (`columns`,
    /// `features`) in one file under [`super::blocks`]; the finalize-stage
    /// blocks own their `compute` there too, while the apex-stage blocks are
    /// built in [`super::apex_finding`].
    pub struct ScoringFields {
        pub identity: Identity,
        pub rt: Rt,
        pub mobility: Mobility,
        pub primary: PrimaryScores,
        pub split: SplitProduct,
        pub features: ApexFeatures,
        pub apex_lazy: ApexLazyScores,
        pub secondary_lazy: SecondaryLazyScores,
        pub counts: ApexCounts,
        pub finalize_counts: FinalizeCounts,
        pub intensities: Intensities,
        pub ion_errors: IonErrors,
        pub rel_intensities: RelativeIntensities,
    }
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
            split: inp.apex.split,
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
    pub delta_group: f32,
    pub delta_group_ratio: f32,
    /// Scratch field for CrossValidatedScorer (written during rescore)
    pub(crate) discriminant_score: f32,
    /// Scratch field for q-value assignment
    pub(crate) qvalue: f32,
}

impl CompetedCandidate {
    /// The post-model meta block (used for the ML delta-group features).
    pub(crate) fn result_meta(&self) -> ResultMeta {
        ResultMeta {
            delta_group: self.delta_group,
            delta_group_ratio: self.delta_group_ratio,
            discriminant_score: self.discriminant_score,
            qvalue: self.qvalue,
        }
    }
}

/// After rescoring. Written to Parquet.
#[derive(Debug, Clone, Serialize)]
pub struct FinalResult {
    pub scoring: ScoringFields,
    pub delta_group: f32,
    pub delta_group_ratio: f32,
    pub discriminant_score: f32,
    pub qvalue: f32,
}

impl FinalResult {
    /// Schema/test fixture with zeroed meta fields.
    pub fn sample() -> Self {
        Self {
            scoring: ScoringFields::sample_default(),
            delta_group: 0.0,
            delta_group_ratio: 0.0,
            discriminant_score: 0.0,
            qvalue: 0.0,
        }
    }

    /// The post-model meta block (used for the Parquet meta columns).
    pub(crate) fn result_meta(&self) -> ResultMeta {
        ResultMeta {
            delta_group: self.delta_group,
            delta_group_ratio: self.delta_group_ratio,
            discriminant_score: self.discriminant_score,
            qvalue: self.qvalue,
        }
    }

    /// Value-free Parquet schema: scoring blocks (composition order) then the
    /// post-model meta block — mirrors `emit_row`'s column order exactly.
    pub fn column_schema(o: &mut SchemaSink) {
        ScoringFields::push_column_schema(o);
        <ResultMeta as ScoreBlock>::column_schema(o);
    }
}

// ---------------------------------------------------------------------------
// Stage conversions
// ---------------------------------------------------------------------------

impl ScoredCandidate {
    /// Convert into a `CompetedCandidate` with the given delta-group values.
    ///
    /// Items that are alone in their group (no competitor) should pass
    /// `f32::NAN` for both deltas.
    pub fn into_competed(self, delta_group: f32, delta_group_ratio: f32) -> CompetedCandidate {
        CompetedCandidate {
            scoring: self.scoring,
            delta_group,
            delta_group_ratio,
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
            delta_group: self.delta_group,
            delta_group_ratio: self.delta_group_ratio,
            discriminant_score: self.discriminant_score,
            qvalue: self.qvalue,
        }
    }
}
