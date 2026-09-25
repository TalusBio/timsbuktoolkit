use timsquery::SpectralCollector;
use timsquery::utils::constants::C13_C12_MASS_DIFF;

use crate::IonAnnot;

/// Shift fragment labels and m/z in the collector's reusable buffers.
pub fn shift_fragment_isotopes(dst: &mut SpectralCollector<IonAnnot, f32>, offset: i8) {
    for (k, v) in dst
        .fragment_labels
        .iter_mut()
        .zip(dst.fragment_mzs.iter_mut())
    {
        // The error names the representable range, so a library carrying an
        // isotope offset too close to the ceiling says what the limit is rather
        // than asserting the situation is impossible.
        let new_ions = k.try_with_offset_neutrons(offset).unwrap_or_else(|e| {
            panic!("fragment {k} cannot take a {offset:+} isotope offset: {e}")
        });
        let mz_offset = (C13_C12_MASS_DIFF / k.get_charge() as f64) * offset as f64;
        *v += mz_offset;
        *k = new_ions;
    }
}
