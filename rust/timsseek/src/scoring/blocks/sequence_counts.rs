//! Independently dispatched residue and modification counts; ML-only output.
use timsquery::chemistry::analyte::PeptideRef;
use timsseek_macros::ScoreBlock;

#[allow(non_snake_case)]
#[derive(ScoreBlock)]
#[score(requires(ResidueSequence))]
pub struct ResidueCounts {
    #[feat(raw, linear = false)]
    peptide_length: f64,
    #[feat(raw, linear = false)]
    aa_count_A: f64,
    #[feat(raw, linear = false)]
    aa_count_C: f64,
    #[feat(raw, linear = false)]
    aa_count_D: f64,
    #[feat(raw, linear = false)]
    aa_count_E: f64,
    #[feat(raw, linear = false)]
    aa_count_F: f64,
    #[feat(raw, linear = false)]
    aa_count_G: f64,
    #[feat(raw, linear = false)]
    aa_count_H: f64,
    #[feat(raw, linear = false)]
    aa_count_I: f64,
    #[feat(raw, linear = false)]
    aa_count_K: f64,
    #[feat(raw, linear = false)]
    aa_count_L: f64,
    #[feat(raw, linear = false)]
    aa_count_M: f64,
    #[feat(raw, linear = false)]
    aa_count_N: f64,
    #[feat(raw, linear = false)]
    aa_count_P: f64,
    #[feat(raw, linear = false)]
    aa_count_Q: f64,
    #[feat(raw, linear = false)]
    aa_count_R: f64,
    #[feat(raw, linear = false)]
    aa_count_S: f64,
    #[feat(raw, linear = false)]
    aa_count_T: f64,
    #[feat(raw, linear = false)]
    aa_count_V: f64,
    #[feat(raw, linear = false)]
    aa_count_W: f64,
    #[feat(raw, linear = false)]
    aa_count_Y: f64,
}
impl ResidueCounts {
    pub fn compute(peptide: PeptideRef<'_>) -> Self {
        let mut counts = [0.0; 26];
        for residue in peptide.residues.bytes() {
            counts[(residue - b'A') as usize] += 1.0;
        }
        let count = |aa: u8| counts[(aa - b'A') as usize];
        Self {
            peptide_length: peptide.residues.len() as f64,
            aa_count_A: count(b'A'),
            aa_count_C: count(b'C'),
            aa_count_D: count(b'D'),
            aa_count_E: count(b'E'),
            aa_count_F: count(b'F'),
            aa_count_G: count(b'G'),
            aa_count_H: count(b'H'),
            aa_count_I: count(b'I'),
            aa_count_K: count(b'K'),
            aa_count_L: count(b'L'),
            aa_count_M: count(b'M'),
            aa_count_N: count(b'N'),
            aa_count_P: count(b'P'),
            aa_count_Q: count(b'Q'),
            aa_count_R: count(b'R'),
            aa_count_S: count(b'S'),
            aa_count_T: count(b'T'),
            aa_count_V: count(b'V'),
            aa_count_W: count(b'W'),
            aa_count_Y: count(b'Y'),
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
