use timsquery::ElutionGroup;

use crate::IonAnnot;

const NEUTRON_MASS: f64 = 1.00866491588;

pub fn isotope_offset_fragments(
    item: &ElutionGroup<IonAnnot>,
    offset: i8,
) -> ElutionGroup<IonAnnot> {
    let mut fragments = Vec::with_capacity(item.fragments.len());
    for (ion, mz) in item.iter_fragments() {
        let new_ions = ion.with_offset_neutrons(offset);
        let mz_offset = (NEUTRON_MASS / ion.get_charge() as f64) * offset as f64;
        let new_mz = mz + mz_offset;
        fragments.push((new_ions, new_mz));
    }

    ElutionGroup {
        id: item.id,
        mobility: item.mobility,
        rt_seconds: item.rt_seconds,
        precursors: item.precursors.clone(),
        fragments: fragments.into(),
    }
}
