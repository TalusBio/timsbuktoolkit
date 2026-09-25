//! Library coordinates never become acquisition coordinates implicitly.
use crate::models::tolerance::RtTolerance;
use crate::{
    ObservedRTSeconds,
    Target,
    TupleRange,
};
use serde::Serialize;

/// The caller must establish acquisition provenance before centering a query.
///
/// ```compile_fail
/// use timsquery::{LibraryRT, RtSelection};
/// let source_seconds = LibraryRT(120.0f32);
/// let selection = RtSelection::Centered(source_seconds);
/// ```
#[derive(Debug, Clone, Copy)]
pub enum RtSelection {
    FullRun,
    Centered(ObservedRTSeconds<f32>),
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct ResolvedRt {
    center: Option<ObservedRTSeconds<f32>>,
    bounds: TupleRange<ObservedRTSeconds<f32>>,
}

impl ResolvedRt {
    pub fn resolve(
        selection: RtSelection,
        tolerance: &RtTolerance,
        acquisition: TupleRange<ObservedRTSeconds<f32>>,
    ) -> Result<Self, crate::DataProcessingError> {
        let invalid = || crate::DataProcessingError::InvalidRtQuery;
        let (start, end) = (acquisition.start().0, acquisition.end().0);
        if !start.is_finite() || !end.is_finite() || start < 0.0 || end < start {
            return Err(invalid());
        }
        let RtSelection::Centered(center) = selection else {
            return Ok(Self {
                center: None,
                bounds: acquisition,
            });
        };
        if !center.0.is_finite() || center.0 < 0.0 {
            return Err(invalid());
        }
        // Preserve the configured tolerance's minute-domain arithmetic before
        // converting the resolved bounds to acquisition seconds.
        let center_minutes = center.0 / 60.0;
        let (left, right) = match tolerance {
            RtTolerance::Unrestricted => {
                return Ok(Self {
                    center: Some(center),
                    bounds: acquisition,
                });
            }
            RtTolerance::Minutes((lo, hi)) => (*lo, *hi),
            RtTolerance::Pct((lo, hi)) => {
                (center_minutes * lo / 100.0, center_minutes * hi / 100.0)
            }
        };
        if !left.is_finite() || !right.is_finite() || left < 0.0 || right < 0.0 {
            return Err(invalid());
        }
        let lo = ((center_minutes - left) * 60.0).max(start);
        let hi = ((center_minutes + right) * 60.0).min(end);
        if lo > hi {
            return Err(invalid());
        }
        Ok(Self {
            center: Some(center),
            bounds: TupleRange::try_new(ObservedRTSeconds(lo), ObservedRTSeconds(hi))
                .map_err(|_| invalid())?,
        })
    }

    pub fn from_mapping(
        selection: RtSelection,
        tolerance: &RtTolerance,
        mapping: &timscentroid::rt_mapping::CycleToRTMapping<
            timscentroid::rt_mapping::MS1CycleIndex,
        >,
    ) -> Result<Self, crate::DataProcessingError> {
        let (lo, hi) = mapping.range_milis();
        Self::resolve(
            selection,
            tolerance,
            TupleRange::try_new(
                ObservedRTSeconds(lo as f32 / 1000.0),
                ObservedRTSeconds(hi as f32 / 1000.0),
            )
            .map_err(|_| crate::DataProcessingError::InvalidRtQuery)?,
        )
    }

    /// Reporting/prediction metadata; filtering uses `bounds` only.
    pub fn center(&self) -> Option<ObservedRTSeconds<f32>> {
        self.center
    }

    pub fn bounds(&self) -> TupleRange<ObservedRTSeconds<f32>> {
        self.bounds
    }

    pub fn range_millis(&self) -> TupleRange<u32> {
        self.bounds
            .map_elems(|s| (s.0 * 1000.0) as u32)
            .expect("validated RT bounds")
    }
}

/// Borrows source geometry only during collector construction or reset.
/// Does not implement `Target`: source RT remains a library coordinate.
pub struct ExtractionQuery<'q, T: Target + ?Sized> {
    source: &'q T,
    rt: ResolvedRt,
    mobility_center: f32,
}

impl<'q, T: Target + ?Sized> ExtractionQuery<'q, T> {
    pub fn new(source: &'q T, rt: ResolvedRt) -> Self {
        Self {
            source,
            rt,
            mobility_center: source.mobility_ook0(),
        }
    }

    pub fn source(&self) -> &'q T {
        self.source
    }

    pub fn rt(&self) -> ResolvedRt {
        self.rt
    }

    pub fn mobility_center(&self) -> f32 {
        self.mobility_center
    }

    pub fn with_mobility(mut self, mobility: f32) -> Self {
        self.mobility_center = mobility;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn acquisition() -> TupleRange<ObservedRTSeconds<f32>> {
        TupleRange::try_new(ObservedRTSeconds(10.0), ObservedRTSeconds(200.0)).unwrap()
    }

    #[test]
    fn resolves_and_clips_once_on_acquisition_axis() {
        let rt = ResolvedRt::resolve(
            RtSelection::Centered(ObservedRTSeconds(20.0)),
            &RtTolerance::Minutes((0.5, 0.5)),
            acquisition(),
        )
        .unwrap();
        assert_eq!(rt.center(), Some(ObservedRTSeconds(20.0)));
        assert_eq!(
            rt.range_millis(),
            TupleRange::try_new(10_000, 50_000).unwrap()
        );
        let full = ResolvedRt::resolve(
            RtSelection::FullRun,
            &RtTolerance::Minutes((0.01, 0.01)),
            acquisition(),
        )
        .unwrap();
        assert_eq!(full.center(), None);
        assert_eq!(full.bounds(), acquisition());
    }

    #[test]
    fn invalid_center_or_disjoint_window_never_becomes_full_run() {
        for center in [f32::NAN, f32::INFINITY, -1.0] {
            assert!(
                ResolvedRt::resolve(
                    RtSelection::Centered(ObservedRTSeconds(center)),
                    &RtTolerance::Unrestricted,
                    acquisition()
                )
                .is_err()
            );
        }
        assert!(
            ResolvedRt::resolve(
                RtSelection::Centered(ObservedRTSeconds(300.0)),
                &RtTolerance::Minutes((0.1, 0.1)),
                acquisition()
            )
            .is_err()
        );
    }

    #[test]
    fn rejects_invalid_tolerances_and_handles_percentages() {
        for tol in [
            RtTolerance::Minutes((-1.0, 1.0)),
            RtTolerance::Pct((f32::NAN, 1.0)),
        ] {
            assert!(
                ResolvedRt::resolve(
                    RtSelection::Centered(ObservedRTSeconds(100.0)),
                    &tol,
                    acquisition()
                )
                .is_err()
            );
        }
        let rt = ResolvedRt::resolve(
            RtSelection::Centered(ObservedRTSeconds(100.0)),
            &RtTolerance::Pct((10.0, 20.0)),
            acquisition(),
        )
        .unwrap();
        assert_eq!(
            rt.range_millis(),
            TupleRange::try_new(90_000, 120_000).unwrap()
        );
    }

    #[test]
    fn collector_owns_no_source_borrow_and_reset_changes_both_center_and_bounds() {
        use crate::{
            OwnedTarget,
            SpectralCollector,
        };
        let mut collector: SpectralCollector<usize, f32> = {
            let source = OwnedTarget::empty_like();
            let rt = ResolvedRt::resolve(
                RtSelection::FullRun,
                &RtTolerance::Unrestricted,
                acquisition(),
            )
            .unwrap();
            SpectralCollector::new(&ExtractionQuery::new(&source, rt))
        };
        let source = OwnedTarget::empty_like();
        let rt = ResolvedRt::resolve(
            RtSelection::Centered(ObservedRTSeconds(100.0)),
            &RtTolerance::Minutes((0.5, 0.5)),
            acquisition(),
        )
        .unwrap();
        collector.reset_with(&ExtractionQuery::new(&source, rt));
        assert_eq!(collector.rt.center(), Some(ObservedRTSeconds(100.0)));
        assert_eq!(
            collector.rt.range_millis(),
            TupleRange::try_new(70_000, 129_999).unwrap()
        );
    }
}
