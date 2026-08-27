use crate::models::target::Target;
use crate::traits::KeyLike;

/// Read-only geometry contract shared by the materialized `Target`
/// and the columnar flyweight `Query<L>`. Method names mirror the aggregator
/// collectors' existing calls so they relax to `&impl QueryGeom` unchanged.
pub trait QueryGeom {
    type Label: KeyLike;

    /// What the source file called this row, when it said. `None` for formats
    /// that carry no id.
    fn source_id(&self) -> Option<crate::models::LibraryId> {
        None
    }

    /// The id a result carries; always present, since callers key results by
    /// it. Not defaulted: a default would need a position accessor on this
    /// trait, and the position belongs to whoever stores the row.
    fn output_id(&self) -> u64;
    fn mono_precursor_mz(&self) -> f64;
    fn precursor_charge(&self) -> u8;
    fn rt_seconds(&self) -> f32;
    fn mobility_ook0(&self) -> f32;
    fn precursor_mz_limits(&self) -> (f64, f64);
    fn precursor_count(&self) -> usize;
    fn fragment_count(&self) -> usize;
    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)>;
    /// Fragment (label, m/z). m/z is BY VALUE: a decoy's is a computed shift
    /// that cannot be borrowed; targets copy their stored f64 (Copy, free).
    fn iter_fragments_refs(&self) -> impl Iterator<Item = (&Self::Label, f64)>;
}

impl<T: KeyLike> QueryGeom for Target<T> {
    type Label = T;

    fn source_id(&self) -> Option<crate::models::LibraryId> {
        Some(crate::models::LibraryId::new(self.id()))
    }

    fn output_id(&self) -> u64 {
        self.id()
    }

    fn mono_precursor_mz(&self) -> f64 {
        self.mono_precursor_mz()
    }

    fn precursor_charge(&self) -> u8 {
        self.precursor_charge()
    }

    fn rt_seconds(&self) -> f32 {
        self.rt_seconds()
    }

    fn mobility_ook0(&self) -> f32 {
        self.mobility_ook0()
    }

    fn precursor_mz_limits(&self) -> (f64, f64) {
        self.precursor_mz_limits()
    }

    fn precursor_count(&self) -> usize {
        self.precursor_count()
    }

    fn fragment_count(&self) -> usize {
        self.fragment_count()
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> {
        self.iter_precursors()
    }

    fn iter_fragments_refs(&self) -> impl Iterator<Item = (&Self::Label, f64)> {
        // Inherent method resolves here (inherent-first), returns (&T, &f64);
        // map to owned m/z. No recursion.
        Target::iter_fragments_refs(self).map(|(k, v)| (k, *v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::target::Target;

    fn via_trait<G: QueryGeom>(g: &G) -> u64 {
        g.output_id()
    }

    #[test]
    fn target_output_id_is_its_source_id() {
        let eg: Target<crate::IonAnnot> = Target::builder()
            .id(7)
            .mobility_ook0(0.75)
            .rt_seconds(1.0)
            .fragment_labels([crate::IonAnnot::try_from("y3").unwrap()].as_slice().into())
            .fragment_mzs(vec![100.0])
            .precursor_labels(tinyvec::tiny_vec!(0))
            .precursor(500.0, 2)
            .try_build()
            .unwrap();
        assert_eq!(via_trait(&eg), 7);
    }

    #[test]
    fn reset_from_accepts_generic_geom() {
        fn reset_via_trait<G: QueryGeom<Label = crate::IonAnnot>>(
            dst: &mut Target<crate::IonAnnot>,
            src: &G,
        ) {
            dst.reset_from(src);
        }
        let src: Target<crate::IonAnnot> = Target::builder()
            .id(3)
            .mobility_ook0(0.5)
            .rt_seconds(2.0)
            .fragment_labels([crate::IonAnnot::try_from("y2").unwrap()].as_slice().into())
            .fragment_mzs(vec![250.0])
            .precursor_labels(tinyvec::tiny_vec!(0))
            .precursor(400.0, 2)
            .try_build()
            .unwrap();
        let mut dst = src.clone();
        reset_via_trait(&mut dst, &src);
        assert_eq!(dst.id(), 3);
    }
}
