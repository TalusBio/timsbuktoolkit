use timsquery::TimsElutionGroup;
use timsquery::utils::constants::NEUTRON_MASS;

use crate::IonAnnot;

/// Buffer-override variant of `isotope_offset_fragments`: reset `dst` from `src`
/// (reusing Vec/TinyVec capacity), then apply a neutron-mass m/z shift and label
/// rewrite to every fragment in place. Zero alloc after warm-up.
///
/// `dst` and `src` must refer to distinct TimsElutionGroup values — typical use
/// is `apply_isotope_offset_fragments_into(&mut worker.isotope_scratch_eg, &item.query, 1)`.
pub fn apply_isotope_offset_fragments_into(
    dst: &mut TimsElutionGroup<IonAnnot>,
    src: &TimsElutionGroup<IonAnnot>,
    offset: i8,
) {
    dst.reset_from(src);
    for (k, v) in dst.iter_fragments_refs_mut() {
        let new_ions = k.try_with_offset_neutrons(offset).expect(
            "Isotope offset overflow - this should never happen with realistic isotope offsets",
        );
        let mz_offset = (NEUTRON_MASS / k.get_charge() as f64) * offset as f64;
        *v += mz_offset;
        *k = new_ions;
    }
}
