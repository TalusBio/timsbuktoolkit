use rustyms::fragment::FragmentType;
use rustyms::model::Location;
use rustyms::spectrum::MassMode;
use rustyms::system::f64::MassOverCharge;
use rustyms::system::mass_over_charge::mz;
use rustyms::system::{
    Charge,
    e,
};
use rustyms::{
    Fragment,
    LinearPeptide,
    Model,
};

use crate::{
    IonAnnot,
    IonParsingError,
};

#[derive(Debug)]
pub struct FragmentMassBuilder {
    pub model: Model,
    pub max_charge: Charge,
}

impl Default for FragmentMassBuilder {
    fn default() -> Self {
        let by_ions = Model {
            a: (Location::None, Vec::new()),
            b: (Location::SkipNC(2, 2), vec![]),
            c: (Location::None, Vec::new()),
            d: (Location::None, Vec::new()),
            v: (Location::None, Vec::new()),
            w: (Location::None, Vec::new()),
            x: (Location::None, Vec::new()),
            y: (Location::SkipNC(2, 2), vec![]),
            z: (Location::None, Vec::new()),
            precursor: vec![],
            // TODO: Fix this hard-coded value
            ppm: MassOverCharge::new::<mz>(20.0),
            glycan_fragmentation: None,
        };
        let max_charge: Charge = Charge::new::<e>(2.0);
        Self {
            model: by_ions,
            max_charge,
        }
    }
}

impl FragmentMassBuilder {
    pub fn fragment_mzs_from_linear_peptide(
        &self,
        peptide: &LinearPeptide,
    ) -> Result<Vec<(IonAnnot, f64, f32)>, IonParsingError> {
        // NOTE: I have to add this retain bc it generates precursor ions even if they are not
        // defined.
        // TODO: return a different error ... this one is very loaded.
        let ions: Vec<Fragment> = peptide
            .generate_theoretical_fragments(self.max_charge, &self.model)
            .into_iter()
            .collect();

        // Does this generate ions above the charge of the precursor?
        ions.into_iter()
            .map(|x| {
                let intensity = match x.ion {
                    FragmentType::Y(_) => 1.0,
                    FragmentType::B(_) => 0.5,
                    _ => 0.01,
                };
                Ok((
                    IonAnnot::from_fragment(x.ion.clone(), x.charge.value as i8, 0)?,
                    // IonAnnot::new(x.ion.clone(), x.charge.abs().value as u8, 0)?,
                    x.mz(MassMode::Monoisotopic).value,
                    intensity,
                ))
            })
            .collect()
    }
}
