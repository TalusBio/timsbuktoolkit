//! Stable index group IDs with an m/z lookup and exact mobility projections.

use super::shadow::bounds;
use crate::errors::TimsSeekError;
use smallvec::SmallVec;
use timscentroid::IndexedTimstofPeaks;
use timscentroid::geometry::{
    QuadrupoleIsolationScheme,
    RingShape,
};
use timscentroid::indexing::IndexedPeakGroup;
use timscentroid::rt_mapping::WindowCycleIndex;

const MAX_SHAPES_PER_GROUP: usize = 8;
const MAX_SHAPE_VERTICES: usize = 16;
pub(super) const MAX_GROUPS_PER_ROW: usize = 4;
const MAX_INTERVALS_PER_GROUP: usize = MAX_SHAPES_PER_GROUP * MAX_SHAPE_VERTICES / 2;

pub(super) type Span = (f64, f64);
pub(super) type Intervals = SmallVec<[Span; MAX_INTERVALS_PER_GROUP]>;
pub(super) type Groups = SmallVec<[u16; MAX_GROUPS_PER_ROW]>;

struct MzEntry {
    id: u16,
    lower: f64,
    upper: f64,
    prefix_upper: f64,
}

pub(super) struct WindowMap<'a> {
    groups: &'a [(
        QuadrupoleIsolationScheme,
        IndexedPeakGroup<WindowCycleIndex>,
    )],
    by_mz: Vec<MzEntry>,
    pub(super) mobility: bool,
}

impl<'a> WindowMap<'a> {
    pub(super) fn new(index: &'a IndexedTimstofPeaks) -> Result<Self, TimsSeekError> {
        let groups = index.ms2_window_groups();
        if groups.len() > usize::from(u16::MAX) {
            return Err(bounds("too many MS2 groups"));
        }
        let mobility = index.mobility_kind().is_filterable();
        let mut by_mz = Vec::with_capacity(groups.len());
        for (id, (geometry, _)) in groups.iter().enumerate() {
            let shapes = geometry.shapes();
            if shapes.len() > MAX_SHAPES_PER_GROUP {
                return Err(bounds("too many isolation shapes in one group"));
            }
            for shape in shapes {
                if let RingShape::Polygon(poly) = shape {
                    let vertices = std::iter::once(poly.exterior())
                        .chain(poly.interiors())
                        .map(|ring| ring.0.len().saturating_sub(1))
                        .sum::<usize>();
                    if vertices > MAX_SHAPE_VERTICES {
                        return Err(bounds("too many vertices in one isolation shape"));
                    }
                }
            }
            let mz = if mobility {
                geometry.fragmented_range()
            } else if let [RingShape::Aabb { mz, .. }] = shapes {
                Some(*mz)
            } else {
                None
            };
            if let Some((lower, upper)) = mz {
                if !(lower.is_finite() && upper.is_finite() && lower <= upper) {
                    return Err(bounds("invalid MS2 group m/z bounds"));
                }
                by_mz.push(MzEntry {
                    id: u16::try_from(id).expect("checked group count"),
                    lower,
                    upper,
                    prefix_upper: upper,
                });
            }
        }
        by_mz.sort_unstable_by(|left, right| {
            left.lower
                .total_cmp(&right.lower)
                .then(left.id.cmp(&right.id))
        });
        let mut prefix_upper = f64::NEG_INFINITY;
        for window in &mut by_mz {
            prefix_upper = prefix_upper.max(window.upper);
            window.prefix_upper = prefix_upper;
        }
        Ok(Self {
            groups,
            by_mz,
            mobility,
        })
    }

    pub(super) fn is_empty(&self) -> bool {
        self.by_mz.is_empty()
    }

    pub(super) fn group_count(&self) -> usize {
        self.groups.len()
    }

    pub(super) fn peaks(&self, group: u16) -> &IndexedPeakGroup<WindowCycleIndex> {
        &self.groups[usize::from(group)].1
    }

    fn overlapping(&self, (lo, hi): Span) -> impl Iterator<Item = &MzEntry> {
        let start = self.by_mz.partition_point(|entry| entry.prefix_upper < lo);
        let end = self.by_mz.partition_point(|entry| entry.lower <= hi);
        self.by_mz[start..end]
            .iter()
            .filter(move |entry| entry.upper >= lo)
    }

    pub(super) fn static_group(&self, mono_mz: f64) -> Option<u16> {
        assert!(!self.mobility);
        self.overlapping((mono_mz, mono_mz))
            .min_by(|left, right| {
                let distance =
                    |entry: &MzEntry| (mono_mz - entry.lower.midpoint(entry.upper)).abs();
                distance(left)
                    .total_cmp(&distance(right))
                    .then(left.lower.total_cmp(&right.lower))
                    .then(left.upper.total_cmp(&right.upper))
                    .then(left.id.cmp(&right.id))
            })
            .map(|entry| entry.id)
    }

    pub(super) fn mobility_groups(
        &self,
        mz: Span,
        mobility: Span,
    ) -> Result<Option<(u16, Groups)>, TimsSeekError> {
        assert!(self.mobility);
        let mut groups = Groups::new();
        let mut best = None;
        let mut best_width = -1.0;
        let mut spans = Intervals::new();
        for entry in self.overlapping(mz) {
            self.project(entry.id, mz, mobility, &mut spans);
            if spans.is_empty() {
                continue;
            }
            if groups.len() == MAX_GROUPS_PER_ROW {
                return Err(bounds("too many MS2 groups for one shadow row"));
            }
            groups.push(entry.id);
            let width: f64 = spans.iter().map(|&(lo, hi)| hi - lo).sum();
            if width > best_width || (width == best_width && best.is_some_and(|id| entry.id < id)) {
                best = Some(entry.id);
                best_width = width;
            }
        }
        assert_eq!(best.is_none(), groups.is_empty());
        Ok(best.map(|id| (id, groups)))
    }

    pub(super) fn project(&self, group: u16, mz: Span, mobility: Span, spans: &mut Intervals) {
        assert!(self.mobility);
        project(&self.groups[usize::from(group)].0, mz, mobility, spans);
    }
}

pub(super) fn project(
    geometry: &QuadrupoleIsolationScheme,
    mz: Span,
    mobility: Span,
    spans: &mut Intervals,
) {
    assert!(mz.0 <= mz.1);
    assert!(mobility.0 <= mobility.1);
    spans.clear();
    for shape in geometry.shapes() {
        match shape {
            RingShape::Aabb { mz: shape_mz, im } => {
                if intersection(mz, *shape_mz).is_some() {
                    push_span(spans, *im, mobility);
                }
            }
            RingShape::Trapezoid { im, mz_lo, mz_hi } => {
                // Clip the ramp parameter to left <= query upper and right >= query lower.
                if let Some(t) = clip_upper((0.0, 1.0), *mz_lo, mz.1)
                    .and_then(|t| clip_upper(t, (-mz_hi.0, -mz_hi.1), -mz.0))
                {
                    push_span(spans, (lerp(*im, t.0), lerp(*im, t.1)), mobility);
                }
            }
            RingShape::Polygon(poly) => {
                let edges = std::iter::once(poly.exterior())
                    .chain(poly.interiors())
                    .flat_map(|ring| ring.0.windows(2))
                    .map(|edge| ((edge[0].x, edge[0].y), (edge[1].x, edge[1].y)));
                project_polygon(edges, mz, mobility, spans);
            }
        }
    }
    let len = merge_spans(spans);
    spans.truncate(len);
    assert!(
        spans
            .iter()
            .all(|&(lo, hi)| mobility.0 <= lo && lo <= hi && hi <= mobility.1)
    );
    assert!(spans.windows(2).all(|pair| pair[0].1 < pair[1].0));
}

fn intersection(left: Span, right: Span) -> Option<Span> {
    let span = (left.0.max(right.0), left.1.min(right.1));
    (span.0 <= span.1).then_some(span)
}

fn clip_upper(t: Span, x: Span, upper: f64) -> Option<Span> {
    match (x.0 <= upper, x.1 <= upper) {
        (true, true) => Some(t),
        (false, false) => None,
        (true, false) => intersection(t, (0.0, (upper - x.0) / (x.1 - x.0))),
        (false, true) => intersection(t, ((upper - x.0) / (x.1 - x.0), 1.0)),
    }
}

fn lerp(values: Span, t: f64) -> f64 {
    if t == 1.0 {
        values.1
    } else {
        values.0 + (values.1 - values.0) * t
    }
}

fn push_span(spans: &mut Intervals, span: Span, mobility: Span) {
    if let Some(span) = intersection(span, mobility) {
        assert!(spans.len() < MAX_INTERVALS_PER_GROUP);
        spans.push(span);
    }
}

fn project_polygon(
    edges: impl Iterator<Item = (Span, Span)>,
    mz: Span,
    mobility: Span,
    spans: &mut Intervals,
) {
    // Each intersecting horizontal scanline either hits an edge in the strip
    // or crosses its entire width. Edge projections cover the first case.
    // Even-odd crossings on the left strip boundary cover the second, including holes.
    const MAX_EDGE_SPANS: usize = MAX_SHAPE_VERTICES + MAX_SHAPE_VERTICES / 2;
    let mut edge_spans = SmallVec::<[Span; MAX_EDGE_SPANS]>::new();
    let mut crossings = SmallVec::<[f64; MAX_SHAPE_VERTICES]>::new();
    for ((x0, y0), (x1, y1)) in edges {
        let x = (x0, x1);
        let y = (y0, y1);
        if let Some(t) =
            clip_upper((0.0, 1.0), x, mz.1).and_then(|t| clip_upper(t, (-x0, -x1), -mz.0))
        {
            let a = lerp(y, t.0);
            let b = lerp(y, t.1);
            assert!(edge_spans.len() < MAX_EDGE_SPANS);
            edge_spans.push((a.min(b), a.max(b)));
        }
        if (x0 <= mz.0) != (x1 <= mz.0) {
            assert!(crossings.len() < MAX_SHAPE_VERTICES);
            crossings.push(lerp(y, (mz.0 - x0) / (x1 - x0)));
        }
    }
    crossings.sort_unstable_by(f64::total_cmp);
    let pairs = crossings.chunks_exact(2);
    assert!(pairs.remainder().is_empty());
    for pair in pairs {
        assert!(edge_spans.len() < MAX_EDGE_SPANS);
        edge_spans.push((pair[0], pair[1]));
    }
    let len = merge_spans(&mut edge_spans);
    assert!(len <= MAX_SHAPE_VERTICES / 2);
    for span in &edge_spans[..len] {
        push_span(spans, *span, mobility);
    }
}

fn merge_spans(spans: &mut [Span]) -> usize {
    spans.sort_unstable_by(|left, right| left.0.total_cmp(&right.0));
    let mut written = 0;
    for index in 0..spans.len() {
        let (lo, hi) = spans[index];
        if written > 0 && lo <= spans[written - 1].1 {
            spans[written - 1].1 = hi.max(spans[written - 1].1);
        } else {
            spans[written] = (lo, hi);
            written += 1;
        }
    }
    assert!(written <= spans.len());
    written
}

pub(super) fn intersects(left: &[Span], right: &[Span]) -> bool {
    let (mut left, mut right) = (left, right);
    for _ in 0..left.len() + right.len() {
        let (Some(&a), Some(&b)) = (left.first(), right.first()) else {
            return false;
        };
        if intersection(a, b).is_some() {
            return true;
        }
        if a.1 < b.0 {
            left = &left[1..];
        } else {
            right = &right[1..];
        }
    }
    false
}

pub(super) fn contains(spans: &[Span], mobility: f64) -> bool {
    let index = spans.partition_point(|span| span.1 < mobility);
    spans.get(index).is_some_and(|span| span.0 <= mobility)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polygon_projection_keeps_holes_that_span_the_query_strip() {
        let exterior = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)];
        let hole = [(1.0, 1.0), (3.0, 1.0), (3.0, 2.0), (1.0, 2.0), (1.0, 1.0)];
        let edges = exterior
            .windows(2)
            .chain(hole.windows(2))
            .map(|edge| (edge[0], edge[1]));
        let mut spans = Intervals::new();
        project_polygon(edges, (1.5, 2.5), (0.5, 3.5), &mut spans);
        assert_eq!(spans.as_slice(), &[(0.5, 1.0), (2.0, 3.5)]);
        assert!(!contains(&spans, 1.5));
    }
}
