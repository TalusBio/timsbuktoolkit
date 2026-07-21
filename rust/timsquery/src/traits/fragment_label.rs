use crate::IonAnnot;
use crate::traits::KeyLike;
use std::sync::Arc;

/// Capability to apply a decoy m/z shift to a fragment label's m/z.
///
/// Implemented by *every* label type an arena can hold so that a later
/// label-generic flyweight compiles. Labels without ion chemistry apply the
/// identity shift.
pub trait DecoyShift {
    /// Decoy m/z shift. Identity for labels without ion chemistry.
    fn decoy_shift_mz(&self, mz: f64, shift: f64) -> f64;
}

/// Marker for labels that carry ion chemistry (series / ordinal / charge).
///
/// Reserved as the bound for future annotation-dependent scores (e.g. a
/// "b-ion count" feature). Implemented only by ion-annotated labels
/// (`IonAnnot`); arbitrary string labels are `DecoyShift` but not
/// `FragmentLabel`. NOTE: decoy m/z shifting is gated on `DecoyShift`, not this
/// trait — a decoy variant computes its shift through `DecoyShift::decoy_shift_mz`,
/// which every label implements.
// TODO(fragment-features): add series() accessor once a standalone series
// enum exists (today series and ordinal are fused in `IonSeriesOrdinal`), and
// wire this bound into the score that consumes it.
pub trait FragmentLabel: KeyLike + DecoyShift {
    fn try_get_ordinal(&self) -> Option<u8>;
    fn get_charge(&self) -> i8;
}

impl DecoyShift for IonAnnot {
    fn decoy_shift_mz(&self, mz: f64, shift: f64) -> f64 {
        // Verbatim `create_mass_shifted_decoy` rule: shift ordinal>2 (and
        // ordinal-less) fragments by `shift / charge`; leave ordinal<=2 alone.
        if IonAnnot::try_get_ordinal(self).is_none_or(|o| o > 2) {
            mz + shift / IonAnnot::get_charge(self) as f64
        } else {
            mz
        }
    }
}

impl FragmentLabel for IonAnnot {
    fn try_get_ordinal(&self) -> Option<u8> {
        IonAnnot::try_get_ordinal(self)
    }

    fn get_charge(&self) -> i8 {
        IonAnnot::get_charge(self)
    }
}

/// Arbitrary string labels carry no ion chemistry: the fragment decoy m/z shift
/// is identity. This impl exists so the generic flyweight `Query<_, Arc<str>>`
/// compiles; the readers only ever build string arenas with `None`/`Passthrough`
/// decoy strategies, so an ordinal-based fragment shift is never requested here.
impl DecoyShift for Arc<str> {
    fn decoy_shift_mz(&self, mz: f64, _shift: f64) -> f64 {
        mz
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ionannot_shifts_high_ordinal_only() {
        let y1 = IonAnnot::try_from("y1").unwrap();
        let y8 = IonAnnot::try_from("y8").unwrap();
        assert!((y1.decoy_shift_mz(200.0, 14.0) - 200.0).abs() < 1e-9); // ord 1 -> no shift
        assert!((y8.decoy_shift_mz(896.5, 14.0) - (896.5 + 14.0)).abs() < 1e-9); // ord 8, z1
    }

    #[test]
    fn string_shift_is_identity() {
        let s: Arc<str> = "whatever".into();
        assert!((s.decoy_shift_mz(500.0, 14.0) - 500.0).abs() < 1e-12);
    }
}
