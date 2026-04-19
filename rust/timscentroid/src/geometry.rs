use geo::algorithm::line_intersection::{
    LineIntersection,
    line_intersection,
};
use geo::{
    Coord,
    Intersects,
    MultiPolygon,
    Polygon,
    Rect,
    SimplifyVwPreserve,
    coord,
};
use serde::{
    Deserialize,
    Serialize,
};

/// Utilities to convert quad windows to geometries.
use timsrust::QuadrupoleSettings;

/// Holds the information of the quadrupole isolation windows.
/// Internally stored as a list of pre-classified `RingShape` entries so
/// the hot-path `intersects` / `intersects_ranges` calls stay alloc-free
/// on the common rectangle and ramped-trapezoid cases and only fall back
/// to the general polygon path when a ring is something exotic.
///
/// # On-disk format
///
/// Serialization emits a bare JSON array of `RingShape` (the v2 format).
/// Deserialization tolerates both v2 (array) and the legacy v1
/// (`{"inner": <MultiPolygon>}`) shape: v1 caches load with a
/// `tracing::warn!` and a one-time reclassification pass.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(into = "Vec<RingShape>", from = "QuadWire")]
pub struct QuadrupoleIsolationScheme {
    ring_shapes: Vec<RingShape>,
}

/// Ring classification. Cheap-to-intersect shapes (AABB, Trapezoid) are
/// handled inline with pure scalar math; anything else falls through to
/// the `Polygon` fallback which replays the legacy `geo` path.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum RingShape {
    /// Axis-aligned rectangle. Matches static DIA isolation windows.
    Aabb { mz: (f64, f64), im: (f64, f64) },
    /// Quadrilateral with two horizontal edges (at `im.0` and `im.1`)
    /// and two slanted sides. Matches diaPASEF ramped windows.
    /// - `mz_lo = (mz at im.0, mz at im.1)` for the left slanted side.
    /// - `mz_hi = (mz at im.0, mz at im.1)` for the right slanted side.
    Trapezoid {
        im: (f64, f64),
        mz_lo: (f64, f64),
        mz_hi: (f64, f64),
    },
    /// Generic polygon fallback. Not expected on current instruments —
    /// loading one fires a single-shot `tracing::warn!` asking the user
    /// to report the ring shape so we can extend the classifier.
    Polygon(Polygon<f64>),
}

// ---------------------------------------------------------------------------
// Serde wire: v2-then-v1 fallback
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum QuadWire {
    /// v2 on-disk shape: bare JSON array of `RingShape`.
    V2(Vec<RingShape>),
    /// v1 legacy shape: `{"inner": <MultiPolygon>}`. Kept only for
    /// backwards-compatibility with caches written before the fast-path
    /// refactor. Upgraded on the fly with a `tracing::warn!`.
    V1(QuadrupoleIsolationSchemeV1),
}

#[derive(Debug, Clone, Deserialize)]
struct QuadrupoleIsolationSchemeV1 {
    inner: MultiPolygon<f64>,
}

impl From<QuadrupoleIsolationScheme> for Vec<RingShape> {
    fn from(q: QuadrupoleIsolationScheme) -> Self {
        q.ring_shapes
    }
}

impl From<QuadWire> for QuadrupoleIsolationScheme {
    fn from(w: QuadWire) -> Self {
        match w {
            QuadWire::V2(ring_shapes) => Self { ring_shapes },
            QuadWire::V1(v1) => {
                tracing::warn!(
                    "Loaded v1 QuadrupoleIsolationScheme (raw MultiPolygon). \
                     Upgrading to v2 in memory by classifying rings. Rebuild \
                     the cache to skip this at every load."
                );
                Self {
                    ring_shapes: classify_rings(&v1.inner),
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Ring classification
// ---------------------------------------------------------------------------

fn classify_rings(mp: &MultiPolygon<f64>) -> Vec<RingShape> {
    mp.0.iter().map(classify_ring).collect()
}

/// One-shot warning when a ring lands on the generic polygon fallback.
/// Fires at most once per process; subsequent fallbacks are silent to
/// avoid spamming the log. Includes the offending ring's coords so users
/// have something concrete to paste into an issue.
static POLYGON_FALLBACK_WARNED: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);

fn warn_polygon_fallback_once(poly: &Polygon<f64>) {
    use std::sync::atomic::Ordering;
    if POLYGON_FALLBACK_WARNED
        .compare_exchange(false, true, Ordering::Relaxed, Ordering::Relaxed)
        .is_ok()
    {
        let coords: Vec<(f64, f64)> = poly.exterior().coords().map(|c| (c.x, c.y)).collect();
        let interiors_n = poly.interiors().len();
        tracing::warn!(
            exterior = ?coords,
            interiors = interiors_n,
            "QuadrupoleIsolationScheme ring did not match AABB or Trapezoid fast-paths \
             and fell back to full polygon intersect. This hurts hot-loop performance. \
             Please report the ring at https://github.com/TalusBio/timsbuktoolkit/issues \
             so we can add a fast-path variant."
        );
    }
}

fn classify_ring(poly: &Polygon<f64>) -> RingShape {
    if !poly.interiors().is_empty() {
        warn_polygon_fallback_once(poly);
        return RingShape::Polygon(poly.clone());
    }
    let coords: Vec<Coord<f64>> = poly.exterior().coords().copied().collect();
    // Expect a closed ring: 5 coords with coords[0] == coords[4].
    if coords.len() != 5 || coords[0] != coords[4] {
        warn_polygon_fallback_once(poly);
        return RingShape::Polygon(poly.clone());
    }
    let pts = &coords[..4];

    // AABB: exactly 2 distinct x, 2 distinct y, 2 vertices at each extreme.
    let xs = [pts[0].x, pts[1].x, pts[2].x, pts[3].x];
    let ys = [pts[0].y, pts[1].y, pts[2].y, pts[3].y];
    let (x_min, x_max) = minmax(&xs);
    let (y_min, y_max) = minmax(&ys);
    let x_pair_only = xs.iter().all(|&x| x == x_min || x == x_max);
    let y_pair_only = ys.iter().all(|&y| y == y_min || y == y_max);
    if x_pair_only && y_pair_only && x_min != x_max && y_min != y_max {
        let n_at_xmin = xs.iter().filter(|&&x| x == x_min).count();
        let n_at_ymin = ys.iter().filter(|&&y| y == y_min).count();
        if n_at_xmin == 2 && n_at_ymin == 2 {
            return RingShape::Aabb {
                mz: (x_min, x_max),
                im: (y_min, y_max),
            };
        }
    }

    // Trapezoid: exactly 2 distinct y values (horizontal top + bottom),
    // 2 vertices at each.
    let ys_unique = unique_sorted(&ys);
    if ys_unique.len() == 2 {
        let im_lo = ys_unique[0];
        let im_hi = ys_unique[1];
        let mut mz_at_lo: Vec<f64> = pts.iter().filter(|p| p.y == im_lo).map(|p| p.x).collect();
        let mut mz_at_hi: Vec<f64> = pts.iter().filter(|p| p.y == im_hi).map(|p| p.x).collect();
        if mz_at_lo.len() == 2 && mz_at_hi.len() == 2 {
            mz_at_lo.sort_by(|a, b| a.partial_cmp(b).unwrap());
            mz_at_hi.sort_by(|a, b| a.partial_cmp(b).unwrap());
            return RingShape::Trapezoid {
                im: (im_lo, im_hi),
                mz_lo: (mz_at_lo[0], mz_at_hi[0]),
                mz_hi: (mz_at_lo[1], mz_at_hi[1]),
            };
        }
    }

    warn_polygon_fallback_once(poly);
    RingShape::Polygon(poly.clone())
}

fn minmax(xs: &[f64]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &x in xs {
        if x < lo {
            lo = x;
        }
        if x > hi {
            hi = x;
        }
    }
    (lo, hi)
}

fn unique_sorted(xs: &[f64]) -> Vec<f64> {
    let mut v: Vec<f64> = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v.dedup();
    v
}

// ---------------------------------------------------------------------------
// Intersection primitives
// ---------------------------------------------------------------------------

#[inline]
fn overlap(a: (f64, f64), b: (f64, f64)) -> bool {
    !(a.1 < b.0 || b.1 < a.0)
}

/// Trapezoid–rect overlap. Clips the trapezoid's im range to the query im
/// range, then checks the mz band at the clipped endpoints. Conservative:
/// returns true iff the query rect overlaps the axis-aligned bounding band
/// at either clipped endpoint. Pure scalar, no heap.
fn trapezoid_rect_overlap(
    trap_im: (f64, f64),
    mz_lo: (f64, f64),
    mz_hi: (f64, f64),
    rect_im: (f64, f64),
    rect_mz: (f64, f64),
) -> bool {
    let im0 = trap_im.0.max(rect_im.0);
    let im1 = trap_im.1.min(rect_im.1);
    if im0 > im1 {
        return false;
    }
    let span = trap_im.1 - trap_im.0;
    if span <= 0.0 {
        return false;
    }
    let t0 = (im0 - trap_im.0) / span;
    let t1 = (im1 - trap_im.0) / span;
    let lo0 = mz_lo.0 + (mz_lo.1 - mz_lo.0) * t0;
    let lo1 = mz_lo.0 + (mz_lo.1 - mz_lo.0) * t1;
    let hi0 = mz_hi.0 + (mz_hi.1 - mz_hi.0) * t0;
    let hi1 = mz_hi.0 + (mz_hi.1 - mz_hi.0) * t1;
    let trap_lo = lo0.min(lo1);
    let trap_hi = hi0.max(hi1);
    !(trap_hi < rect_mz.0 || rect_mz.1 < trap_lo)
}

struct MinMax {
    inner: Option<(f64, f64)>,
}

impl MinMax {
    fn update(&mut self, v: f64) {
        if let Some((min, max)) = self.inner {
            if v < min {
                self.inner = Some((v, max));
            }
            if v > max {
                self.inner = Some((min, v));
            }
        } else {
            self.inner = Some((v, v));
        }
    }
}

impl QuadrupoleIsolationScheme {
    /// Returns true if the given mz and im ranges intersect any of the
    /// isolation rings.
    pub fn intersects(&self, mz_range: (f64, f64), im_range: (f64, f64)) -> bool {
        for shape in &self.ring_shapes {
            let hit = match shape {
                RingShape::Aabb { mz, im } => overlap(*mz, mz_range) && overlap(*im, im_range),
                RingShape::Trapezoid { im, mz_lo, mz_hi } => {
                    trapezoid_rect_overlap(*im, *mz_lo, *mz_hi, im_range, mz_range)
                }
                RingShape::Polygon(poly) => {
                    let rect = Rect::new(
                        coord! { x: mz_range.0, y: im_range.0 },
                        coord! { x: mz_range.1, y: im_range.1 },
                    );
                    poly.intersects(&rect)
                }
            };
            if hit {
                return true;
            }
        }
        false
    }

    /// Returns the im sub-range where the query rectangle overlaps any
    /// ring. Conservative over-approximation for the AABB and Trapezoid
    /// fast paths (returns the clamped im span; narrower than strictly
    /// necessary only in degenerate trapezoid cases, never wider).
    pub fn intersects_ranges(
        &self,
        mz_range: (f64, f64),
        im_range: (f64, f64),
    ) -> Option<(f64, f64)> {
        let mut out = MinMax { inner: None };
        for shape in &self.ring_shapes {
            match shape {
                RingShape::Aabb { mz, im } => {
                    if overlap(*mz, mz_range) && overlap(*im, im_range) {
                        out.update(im.0.max(im_range.0));
                        out.update(im.1.min(im_range.1));
                    }
                }
                RingShape::Trapezoid { im, mz_lo, mz_hi } => {
                    if trapezoid_rect_overlap(*im, *mz_lo, *mz_hi, im_range, mz_range) {
                        out.update(im.0.max(im_range.0));
                        out.update(im.1.min(im_range.1));
                    }
                }
                RingShape::Polygon(poly) => {
                    let rect = Rect::new(
                        coord! { x: mz_range.0, y: im_range.0 },
                        coord! { x: mz_range.1, y: im_range.1 },
                    );
                    if !poly.intersects(&rect) {
                        continue;
                    }
                    for box_line in rect.to_lines() {
                        if poly.intersects(&box_line.start_point()) {
                            out.update(box_line.start.y);
                        }
                        if poly.intersects(&box_line.end_point()) {
                            out.update(box_line.end.y);
                        }
                        for poly_line in poly.exterior().lines() {
                            match line_intersection(poly_line, box_line) {
                                None => continue,
                                Some(LineIntersection::Collinear { intersection }) => {
                                    out.update(intersection.start.y);
                                    out.update(intersection.end.y);
                                }
                                Some(LineIntersection::SinglePoint {
                                    intersection,
                                    is_proper: _,
                                }) => {
                                    out.update(intersection.y);
                                }
                            }
                        }
                    }
                }
            }
        }
        out.inner
    }

    /// Returns the overall mz span covered by all rings. Used once at
    /// index-load time by the scorer to clip precursor queries to the
    /// fragmented m/z range.
    pub fn fragmented_range(&self) -> Option<(f64, f64)> {
        let mut minmax = MinMax { inner: None };
        for shape in &self.ring_shapes {
            match shape {
                RingShape::Aabb { mz, .. } => {
                    minmax.update(mz.0);
                    minmax.update(mz.1);
                }
                RingShape::Trapezoid { mz_lo, mz_hi, .. } => {
                    minmax.update(mz_lo.0);
                    minmax.update(mz_lo.1);
                    minmax.update(mz_hi.0);
                    minmax.update(mz_hi.1);
                }
                RingShape::Polygon(poly) => {
                    for c in poly.exterior().coords() {
                        minmax.update(c.x);
                    }
                }
            }
        }
        minmax.inner
    }

    pub(crate) fn from_quad(quad: &QuadrupoleSettings, ims_converter: impl Fn(f64) -> f64) -> Self {
        let xxyys = quad_to_xxyy(quad, ims_converter);
        let mp = xxyys_to_geometry(&xxyys);
        Self {
            ring_shapes: classify_rings(&mp),
        }
    }

    /// Creates a QuadrupoleIsolationScheme from an iterator of
    /// (mz_start, mz_end, im_start, im_end) tuples. Each tuple represents
    /// a quadrupole window.
    pub fn from_xxyy<T: Iterator<Item = (f64, f64, f64, f64)>>(iter: T) -> Self {
        let xxyys: Vec<(f64, f64, f64, f64)> = iter.collect();
        let mp = xxyys_to_geometry(&xxyys);
        Self {
            ring_shapes: classify_rings(&mp),
        }
    }
}

fn connect_edges(left: &[(f64, f64)], right: &[(f64, f64)]) -> Polygon<f64> {
    let mut all_points = left.to_vec();
    all_points.extend(right.iter().rev());
    Polygon::new(all_points.into(), vec![])
}

fn quad_to_xxyy(
    quad: &QuadrupoleSettings,
    ims_converter: impl Fn(f64) -> f64,
) -> Vec<(f64, f64, f64, f64)> {
    assert!(quad.scan_starts.len() == quad.scan_ends.len());
    assert!(quad.scan_starts.len() == quad.isolation_mz.len());
    assert!(quad.scan_starts.len() == quad.isolation_width.len());
    assert!(quad.scan_starts.len() == quad.collision_energy.len());
    let mut result = Vec::new();
    for i in 0..quad.scan_starts.len() {
        let scan_start = quad.scan_starts[i];
        let scan_end = quad.scan_ends[i];
        let im_start = ims_converter(scan_start as f64);
        let im_end = ims_converter(scan_end as f64);
        let mz_center = quad.isolation_mz[i];
        let mz_width = quad.isolation_width[i];
        let mz_start = mz_center - mz_width / 2.0;
        let mz_end = mz_center + mz_width / 2.0;
        result.push((mz_start, mz_end, im_start, im_end));
    }
    result
}

// We let geo-rs do the heavy lifting of geometry operations at build
// time so we just make possibly a lot of boxes and then let it simplify
// the geometry for us.
fn xxyys_to_geometry(xxyys: &[(f64, f64, f64, f64)]) -> MultiPolygon<f64> {
    if xxyys.is_empty() {
        return MultiPolygon(vec![]);
    }
    let mut polygons = Vec::new();
    let mut curr_left_edge = Vec::new();
    let mut curr_right_edge = Vec::new();

    let mut last_mobility_end: Option<f64> = None;
    let mut last_quad_range: Option<(f64, f64)> = None;

    for (mz_start, mz_end, im_start, im_end) in xxyys.iter().copied() {
        let mut flush_polygon = false;

        if let Some((last_mz_start, last_mz_end)) = last_quad_range {
            // If the current mz range does not overlap with the last one,
            // we need to flush the current polygon.
            if mz_start > last_mz_end || mz_end < last_mz_start {
                flush_polygon = true;
            }
        }

        if let Some(last_end) = last_mobility_end {
            // If the current im_start is less than the last_im_end,
            // we have an overlap in the IMS dimension.
            // We also need the quads to overlap in the mz dimension
            if im_start != last_end {
                flush_polygon = true;
            }
        }

        if flush_polygon {
            // We have a gap, close the previous polygon
            if !curr_left_edge.is_empty() && !curr_right_edge.is_empty() {
                // Close the polygon by connecting the last right edge to the first left edge
                polygons.push(connect_edges(&curr_left_edge, &curr_right_edge));
            }
            curr_left_edge.clear();
            curr_right_edge.clear();
        }

        curr_left_edge.push((mz_start, im_start));
        curr_left_edge.push((mz_start, im_end));

        curr_right_edge.push((mz_end, im_start));
        curr_right_edge.push((mz_end, im_end));
        last_mobility_end = Some(im_end);
        last_quad_range = Some((mz_start, mz_end));
    }

    // Close the last polygon if needed
    if !curr_left_edge.is_empty() && !curr_right_edge.is_empty() {
        polygons.push(connect_edges(&curr_left_edge, &curr_right_edge));
    }

    // Here we could merge overlapping polygons if needed.
    // For simplicity, we return the list of polygons as is.
    let multy_poly = geo::MultiPolygon(polygons);
    multy_poly.simplify_vw_preserve(0.01)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_ims_converter(x: f64) -> f64 {
        // They go in inverted order from 0 == 1.5
        // to 708 == 0.5
        1.5 - (x / 708.0)
    }

    #[test]
    fn test_quad_to_geometry_diapasef() {
        let dia_pasef_window = QuadrupoleSettings {
            index: 0,
            // Note scan starts and ends are always in increasing order
            // AND never contain overlaps.
            // (if window 1 is 1-100, window 2 cannot start before 101)
            scan_starts: vec![0, 220],
            scan_ends: vec![200, 300],
            isolation_mz: vec![900.0, 500.0],
            isolation_width: vec![40.0, 40.0],
            collision_energy: vec![20.0, 22.0],
        };
        let quad = QuadrupoleIsolationScheme::from_quad(&dia_pasef_window, dummy_ims_converter);
        assert_eq!(quad.ring_shapes.len(), 2);
    }

    #[test]
    fn test_quad_to_geometry_diagonalpasef_classifies_as_trapezoid() {
        let dia_pasef_window = QuadrupoleSettings {
            index: 0,
            scan_starts: vec![0, 1, 2, 3, 4, 5],
            scan_ends: vec![1, 2, 3, 4, 5, 6],
            // We start isolating the 90-110mz window ar 1.5 ims
            // end in the 85-105mz window a bit under that
            isolation_mz: vec![100., 99., 98., 97., 96., 95.],
            isolation_width: vec![20., 20., 20., 20., 20., 20.],
            collision_energy: vec![20.0, 22.0, 24.0, 26.0, 28.0, 30.0],
        };
        let quad = QuadrupoleIsolationScheme::from_quad(&dia_pasef_window, dummy_ims_converter);
        assert_eq!(quad.ring_shapes.len(), 1);
        let max_ims = dummy_ims_converter(0.0);
        let min_ims = dummy_ims_converter(6.0);
        match &quad.ring_shapes[0] {
            RingShape::Trapezoid { im, mz_lo, mz_hi } => {
                assert!((im.0 - min_ims).abs() < 1e-9);
                assert!((im.1 - max_ims).abs() < 1e-9);
                assert!((mz_lo.0 - 85.0).abs() < 1e-3);
                assert!((mz_lo.1 - 90.0).abs() < 1e-3);
                assert!((mz_hi.0 - 105.0).abs() < 1e-3);
                assert!((mz_hi.1 - 110.0).abs() < 1e-3);
            }
            other => panic!("expected Trapezoid, got {:?}", other),
        }
    }

    #[test]
    fn test_intresection_diapasef() {
        let dia_pasef_window = QuadrupoleSettings {
            index: 0,
            scan_starts: vec![0, 220],
            scan_ends: vec![200, 300],
            isolation_mz: vec![900.0, 500.0],
            isolation_width: vec![40.0, 40.0],
            collision_energy: vec![20.0, 22.0],
        };
        let quad_polys =
            QuadrupoleIsolationScheme::from_quad(&dia_pasef_window, dummy_ims_converter);
        // This range should intersect with the first window
        // mz: 880-920, im: 1.0-1.4...ish
        assert!(quad_polys.intersects((880.0, 920.0), (1.0, 1.4)));
    }

    #[test]
    fn test_v1_json_upgrade_classifies_as_aabb() {
        // One of the MS2 window groups from the Hela bench.
        let v1_json = r#"
        {
          "inner": [
            {
              "exterior": [
                {"x": 975.0, "y": 1.3710553037445032},
                {"x": 975.0, "y": 1.1806094858619247},
                {"x": 1000.0, "y": 1.1806094858619247},
                {"x": 1000.0, "y": 1.3710553037445032},
                {"x": 975.0, "y": 1.3710553037445032}
              ],
              "interiors": []
            },
            {
              "exterior": [
                {"x": 775.0, "y": 1.1806094858619247},
                {"x": 775.0, "y": 0.9908424397939822},
                {"x": 800.0, "y": 0.9908424397939822},
                {"x": 800.0, "y": 1.1806094858619247},
                {"x": 775.0, "y": 1.1806094858619247}
              ],
              "interiors": []
            }
          ]
        }
        "#;
        let quad: QuadrupoleIsolationScheme = serde_json::from_str(v1_json).unwrap();
        assert_eq!(quad.ring_shapes.len(), 2);
        for shape in &quad.ring_shapes {
            assert!(matches!(shape, RingShape::Aabb { .. }), "got {:?}", shape);
        }
    }

    #[test]
    fn test_v2_json_roundtrip_is_bare_array() {
        let quad = QuadrupoleIsolationScheme {
            ring_shapes: vec![
                RingShape::Aabb {
                    mz: (500.0, 520.0),
                    im: (0.5, 1.5),
                },
                RingShape::Aabb {
                    mz: (520.0, 540.0),
                    im: (0.5, 1.5),
                },
            ],
        };
        let json = serde_json::to_string(&quad).unwrap();
        assert!(json.starts_with('['), "v2 should be a bare array: {json}");
        let reparsed: QuadrupoleIsolationScheme = serde_json::from_str(&json).unwrap();
        assert_eq!(reparsed.ring_shapes, quad.ring_shapes);
    }

    #[test]
    fn test_aabb_intersects_parity_vs_geo() {
        // Same Hela fixture as the v1 upgrade test.
        let v1_json = r#"
        {
          "inner": [
            {"exterior": [
              {"x": 975.0, "y": 1.3710553037445032},
              {"x": 975.0, "y": 1.1806094858619247},
              {"x": 1000.0, "y": 1.1806094858619247},
              {"x": 1000.0, "y": 1.3710553037445032},
              {"x": 975.0, "y": 1.3710553037445032}
            ], "interiors": []}
          ]
        }
        "#;
        let quad: QuadrupoleIsolationScheme = serde_json::from_str(v1_json).unwrap();
        let v1: QuadrupoleIsolationSchemeV1 = serde_json::from_str(v1_json).unwrap();
        let mp = &v1.inner;
        let cases: &[((f64, f64), (f64, f64))] = &[
            ((980.0, 990.0), (1.2, 1.3)),
            ((970.0, 980.0), (1.18, 1.19)),
            ((850.0, 900.0), (1.2, 1.3)),
            ((700.0, 700.0), (0.5, 0.55)),
        ];
        for (mz, im) in cases {
            let fast = quad.intersects(*mz, *im);
            let rect = Rect::new(coord! { x: mz.0, y: im.0 }, coord! { x: mz.1, y: im.1 });
            let slow = mp.intersects(&rect);
            assert_eq!(fast, slow, "mz={:?} im={:?}", mz, im);
        }
    }
}
