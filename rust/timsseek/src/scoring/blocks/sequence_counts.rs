//! Independently dispatched residue and modification counts; ML-only output.
use timsquery::chemistry::analyte::PeptideRef;
use timsseek_macros::ScoreBlock;

use super::{
    ColSink,
    NameSink,
};
use crate::scoring::plan::Requirement;
use timsquery::chemistry::CANONICAL_AA_LETTERS;

pub struct ResidueCounts([f64; Self::NONLINEAR_LEN]);

impl ResidueCounts {
    pub const NONLINEAR_LEN: usize = CANONICAL_AA_LETTERS.len() + 1;

    pub fn compute(peptide: PeptideRef<'_>) -> Self {
        let mut counts = [0.0; 26];
        for residue in peptide.residues.bytes() {
            counts[(residue - b'A') as usize] += 1.0;
        }
        let mut values = [0.0; Self::NONLINEAR_LEN];
        values[0] = peptide.residues.len() as f64;
        for (i, &aa) in CANONICAL_AA_LETTERS.iter().enumerate() {
            values[i + 1] = counts[(aa - b'A') as usize];
        }
        Self(values)
    }

    pub fn nonlinear_feature_array(&self) -> [f64; Self::NONLINEAR_LEN] {
        self.0
    }
}

impl super::ScoreBlock for ResidueCounts {
    fn requirement() -> Option<Requirement> {
        Some(Requirement::ResidueSequence)
    }

    fn columns(&self, _: &mut ColSink) {}

    fn nonlinear_feature_names(out: &mut NameSink) {
        out.push("peptide_length");
        for aa in CANONICAL_AA_LETTERS {
            out.push(&format!("aa_count_{}", aa as char));
        }
    }
}

#[derive(ScoreBlock)]
#[score(requires(ModificationCount))]
pub struct ModificationCounts {
    #[feat(raw, linear = false)]
    peptide_n_mods: f64,
}
impl ModificationCounts {
    pub fn compute(peptide: PeptideRef<'_>) -> Self {
        Self {
            peptide_n_mods: peptide
                .modifications
                .known()
                .expect("plan promised complete modification count")
                .len() as f64,
        }
    }
}

#[cfg(test)]
const LEN: usize = ResidueCounts::NONLINEAR_LEN + ModificationCounts::NONLINEAR_LEN;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::sequence::CANONICAL_AA_LETTERS;
    use crate::scoring::blocks::{
        NameSink,
        ScoreBlock,
    };
    use timsquery::chemistry::analyte::Analyte;
    use timsquery::models::capabilities::DecoyPolicy;
    use timsquery::models::{
        Row,
        TargetCapabilities,
        TargetColumnsBuilder,
    };

    #[test]
    fn stored_counts_match_expected_residues_and_modifications() {
        let mut names = NameSink::new();
        ResidueCounts::nonlinear_feature_names(&mut names);
        let names = names.into_names();
        for (sequence, residues, n_mods) in [
            ("PEPTIDEK", "PEPTIDEK", 0),
            ("_AC(UniMod:4)M(UniMod:35)K_", "ACMK", 2),
            ("[UNIMOD:1]-PEP[+15.99]TIDE-[UNIMOD:2]", "PEPTIDE", 3),
            ("M[Oxidation]PEPTIDE", "MPEPTIDE", 1),
            ("PEPTIDE/2", "PEPTIDE", 0),
            ("AXA", "AXA", 0),
        ] {
            let analyte = Analyte::from_sequence(sequence);
            let mut builder = TargetColumnsBuilder::<timsquery::ion::IonAnnot>::with_capabilities(
                TargetCapabilities::default_diann(),
            );
            builder.push_row(Row {
                analyte: analyte.as_input(),
                ..Default::default()
            });
            let geom = builder.seal(DecoyPolicy::Never).unwrap();
            let values = {
                let peptide = geom
                    .analyte(geom.rows().next().unwrap())
                    .peptide
                    .known()
                    .unwrap();
                let mut values = ResidueCounts::compute(peptide)
                    .nonlinear_feature_array()
                    .to_vec();
                values.extend(ModificationCounts::compute(peptide).nonlinear_feature_array());
                values
            };
            assert_eq!(values[0], residues.len() as f64, "{sequence}");
            for (i, &aa) in CANONICAL_AA_LETTERS.iter().enumerate() {
                assert_eq!(
                    values[i + 1],
                    residues.bytes().filter(|&r| r == aa).count() as f64
                );
                assert_eq!(&*names[i + 1], format!("aa_count_{}", aa as char));
            }
            assert_eq!(values[LEN - 1], n_mods as f64, "{sequence}");
        }
    }
}
