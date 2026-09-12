//! Independently dispatched residue and modification counts; ML-only output.
use timsquery::chemistry::analyte::PeptideRef;
use timsseek_macros::ScoreBlock;

use timsquery::chemistry::CANONICAL_AA_LETTERS;

const N_RESIDUES: usize = CANONICAL_AA_LETTERS.len();

const AA_SLOTS: usize = (b'Z' - b'A' + 1) as usize;
const CANONICAL_SLOTS: [usize; N_RESIDUES] = {
    let mut slots = [0; N_RESIDUES];
    let mut i = 0;
    while i < N_RESIDUES {
        slots[i] = (CANONICAL_AA_LETTERS[i] - b'A') as usize;
        i += 1;
    }
    slots
};
const COUNT_NAME_BYTES: [[u8; 10]; N_RESIDUES] = {
    let mut names = [*b"aa_count_A"; N_RESIDUES];
    let mut i = 0;
    while i < N_RESIDUES {
        names[i][9] = CANONICAL_AA_LETTERS[i];
        i += 1;
    }
    names
};
const COUNT_NAMES: [&str; N_RESIDUES] = {
    let mut names = [""; N_RESIDUES];
    let mut i = 0;
    while i < N_RESIDUES {
        names[i] = match std::str::from_utf8(&COUNT_NAME_BYTES[i]) {
            Ok(name) => name,
            Err(_) => panic!("canonical residues must be ASCII"),
        };
        i += 1;
    }
    names
};

#[derive(ScoreBlock)]
#[score(requires(ResidueSequence))]
pub struct ResidueCounts {
    #[feat(raw, linear = false)]
    peptide_length: f64,
    #[feat(raw, linear = false, names = COUNT_NAMES, indices = CANONICAL_SLOTS)]
    aa_count: [f64; AA_SLOTS],
}

impl ResidueCounts {
    pub fn compute(peptide: PeptideRef<'_>) -> Self {
        let mut aa_count = [0.0; AA_SLOTS];
        for residue in peptide.residues.bytes() {
            aa_count[(residue - b'A') as usize] += 1.0;
        }
        Self {
            peptide_length: peptide.residues.len() as f64,
            aa_count,
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
    fn alphabet_storage_retains_noncanonical_slots_without_projecting_them() {
        let analyte = Analyte::from_sequence("AXAZ");
        let mut builder = TargetColumnsBuilder::<timsquery::ion::IonAnnot>::with_capabilities(
            TargetCapabilities::default_diann(),
        );
        builder.push_row(Row {
            analyte: analyte.as_input(),
            ..Default::default()
        });
        let geom = builder.seal(DecoyPolicy::Never).unwrap();
        let counts = ResidueCounts::compute(
            geom.analyte(geom.rows().next().unwrap())
                .peptide
                .known()
                .unwrap(),
        );
        assert_eq!(counts.aa_count.len(), AA_SLOTS);
        assert_eq!(counts.aa_count[(b'X' - b'A') as usize], 1.0);
        assert_eq!(counts.aa_count[(b'Z' - b'A') as usize], 1.0);
        let projected = counts.nonlinear_feature_array();
        assert_eq!(projected[0], 4.0);
        assert_eq!(projected[1..].iter().sum::<f64>(), 2.0);
    }

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
