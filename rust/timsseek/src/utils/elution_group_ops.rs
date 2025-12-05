use timsquery::TimsElutionGroup;

use crate::IonAnnot;

const NEUTRON_MASS: f64 = 1.00866491588;

pub fn isotope_offset_fragments(
    item: &TimsElutionGroup<IonAnnot>,
    offset: i8,
) -> TimsElutionGroup<IonAnnot> {
    let mut out = item.clone();
    for (k, v) in out.iter_fragments_refs_mut() {
        let new_ions = k.try_with_offset_neutrons(offset).expect(
            "Isotope offset overflow - this should never happen with realistic isotope offsets",
        );
        let mz_offset = (NEUTRON_MASS / k.get_charge() as f64) * offset as f64;
        *v += mz_offset;
        *k = new_ions;
    }

    out
}
