use crate::models::elution_group::TimsElutionGroup;
use crate::traits::KeyLike;

/// Read-only geometry contract shared by the materialized `TimsElutionGroup`
/// and the columnar flyweight `Query<L>`. Method names mirror the aggregator
/// collectors' existing calls so they relax to `&impl QueryGeom` unchanged.
pub trait QueryGeom {
    type Label: KeyLike;

    /// Positional library id (`library_id`). For `TimsElutionGroup` this is the
    /// assigned `id` cast to u32; for the arena flyweight it is the target index.
    fn id(&self) -> u32;
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

impl<T: KeyLike> QueryGeom for TimsElutionGroup<T> {
    type Label = T;

    fn id(&self) -> u32 {
        // Fail loud rather than silently truncate: the library_id column and the
        // q-value determinism anchor are u32; an id beyond u32::MAX is a bug we
        // want surfaced, not wrapped.
        u32::try_from(self.id()).expect("TimsElutionGroup id exceeds u32::MAX")
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
        TimsElutionGroup::iter_fragments_refs(self).map(|(k, v)| (k, *v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::elution_group::TimsElutionGroup;

    fn geom_id<G: QueryGeom>(g: &G) -> u32 {
        g.id()
    }

    #[test]
    fn blanket_impl_exposes_id_as_u32() {
        let eg: TimsElutionGroup<crate::IonAnnot> = TimsElutionGroup::builder()
            .id(7)
            .mobility_ook0(0.75)
            .rt_seconds(1.0)
            .fragment_labels([crate::IonAnnot::try_from("y3").unwrap()].as_slice().into())
            .fragment_mzs(vec![100.0])
            .precursor_labels(tinyvec::tiny_vec!(0))
            .precursor(500.0, 2)
            .try_build()
            .unwrap();
        assert_eq!(geom_id(&eg), 7u32);
    }

    #[test]
    fn reset_from_accepts_generic_geom() {
        fn reset_via_trait<G: QueryGeom<Label = crate::IonAnnot>>(
            dst: &mut TimsElutionGroup<crate::IonAnnot>,
            src: &G,
        ) {
            dst.reset_from(src);
        }
        let src: TimsElutionGroup<crate::IonAnnot> = TimsElutionGroup::builder()
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
