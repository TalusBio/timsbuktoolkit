use geo::algorithm::line_intersection::{
    LineIntersection,
    line_intersection,
};
use geo::{
    Intersects,
    MultiPolygon,
    Polygon,
    Rect,
    SimplifyVwPreserve,
    coord,
};

/// Utilities to convert quad windows to geometries.
use timsrust::QuadrupoleSettings;

/// Holds the information of the quadrupole isolation windows.
/// Allows to query intersections and mobility to mz ranges.
/// and to convert im to mz ranges.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct QuadrupoleIsolationScheme {
    inner: MultiPolygon<f64>,
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
    /// Returns true if the given mz and im ranges intersect with the isolation
    /// ranges of quad-mz and ims.
    pub fn intersects(&self, mz_range: (f64, f64), im_range: (f64, f64)) -> bool {
        let rect = Rect::new(
            coord! { x: mz_range.0, y: im_range.0 },
            coord! { x: mz_range.1, y: im_range.1 },
        );
        self.inner.intersects(&rect)
    }

    /// Returns the min-max mz that intersects with the polygons at the given mobility.
    /// If no intersection, returns None.
    pub fn im_to_tof_range(&self, im: f64) -> Option<(f64, f64)> {
        let line = geo::Line::new(coord! { x: f64::MIN, y: im }, coord! { x: f64::MAX, y: im });
        let mut min_mz_intersect = MinMax { inner: None };
        for poly in self.inner.iter() {
            for poly_line in poly.exterior().lines() {
                match line_intersection(poly_line, line) {
                    None => continue,
                    Some(LineIntersection::Collinear { intersection }) => {
                        min_mz_intersect.update(intersection.start.x);
                        min_mz_intersect.update(intersection.end.x);
                    }
                    Some(LineIntersection::SinglePoint {
                        intersection,
                        is_proper: _,
                    }) => {
                        // What does proper mean here?
                        min_mz_intersect.update(intersection.x);
                    }
                }
            }
        }

        min_mz_intersect.inner
    }

    /// Returns the min-max mobility that intersects with the polygons.
    pub fn intersects_ranges(
        &self,
        mz_range: (f64, f64),
        im_range: (f64, f64),
    ) -> Option<(f64, f64)> {
        let rect = Rect::new(
            coord! { x: mz_range.0, y: im_range.0 },
            coord! { x: mz_range.1, y: im_range.1 },
        );
        if !self.inner.intersects(&rect) {
            return None;
        }

        let mut min_im_intersect = MinMax { inner: None };

        for poly in self.inner.iter() {
            for box_line in rect.to_lines() {
                if poly.intersects(&box_line.start_point()) {
                    min_im_intersect.update(box_line.start.y);
                }
                if poly.intersects(&box_line.end_point()) {
                    min_im_intersect.update(box_line.end.y);
                }
                for poly_line in poly.exterior().lines() {
                    match line_intersection(poly_line, box_line) {
                        None => continue,
                        Some(LineIntersection::Collinear { intersection }) => {
                            min_im_intersect.update(intersection.start.y);
                            min_im_intersect.update(intersection.end.y);
                        }
                        Some(LineIntersection::SinglePoint {
                            intersection,
                            is_proper: _,
                        }) => {
                            // What does proper mean here?
                            min_im_intersect.update(intersection.y);
                        }
                    }
                }
            }
        }

        min_im_intersect.inner
    }

    pub fn fragmented_range(&self) -> Option<(f64, f64)> {
        let mut minmax = MinMax { inner: None };
        for poly in self.inner.iter() {
            for coord in poly.exterior().coords() {
                minmax.update(coord.x);
            }
        }
        minmax.inner
    }

    pub(crate) fn from_quad(quad: &QuadrupoleSettings, ims_converter: impl Fn(f64) -> f64) -> Self {
        let xxyys = quad_to_xxyy(quad, ims_converter);
        let geom = xxyys_to_geometry(&xxyys);
        Self { inner: geom }
    }

    /// Creates a QuadrupoleIsolationScheme from an iterator of (mz_start, mz_end, im_start, im_end) tuples.
    /// Each tuple represents a quadrupole window.
    pub fn from_xxyy<T: Iterator<Item = (f64, f64, f64, f64)>>(iter: T) -> Self {
        let xxyys: Vec<(f64, f64, f64, f64)> = iter.collect();
        let geom = xxyys_to_geometry(&xxyys);
        Self { inner: geom }
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

// We let geo-rs do the heavy lifting of geometry operations
// so we just make possibly a lot of boxes and then let it simplify
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
        let xxyys = quad_to_xxyy(&dia_pasef_window, dummy_ims_converter);
        let geom = xxyys_to_geometry(&xxyys);
        println!("Geometry: {:?}", geom);
        assert_eq!(geom.0.len(), 2);
    }

    #[test]
    fn test_quad_to_geometry_diagonalpasef() {
        let dia_pasef_window = QuadrupoleSettings {
            index: 0,
            // Note scan starts and ends are always in increasing order
            // AND never contain overlaps.
            // (if window 1 is 1-100, window 2 cannot start before 101)
            scan_starts: vec![0, 1, 2, 3, 4, 5],
            scan_ends: vec![1, 2, 3, 4, 5, 6],

            // We start isolating the 90-110mz window ar 1.5 ims
            // end in the 85-105mz window a bit under that
            isolation_mz: vec![100., 99., 98., 97., 96., 95.],
            isolation_width: vec![20., 20., 20., 20., 20., 20.],
            collision_energy: vec![20.0, 22.0, 24.0, 26.0, 28.0, 30.0],
        };
        let geom = quad_to_xxyy(&dia_pasef_window, dummy_ims_converter);
        let geom = xxyys_to_geometry(&geom);
        println!("Geometry: {:?}", geom);
        assert_eq!(geom.0.len(), 1);

        // Extract that one polygon
        let poly = &geom.0[0];
        let coords: Vec<(f64, f64)> = poly.exterior().coords().map(|c| (c.x, c.y)).collect();
        println!("Coords: {:?}", coords);
        let max_ims = dia_pasef_window
            .scan_starts
            .first()
            .map(|s| dummy_ims_converter(*s as f64))
            .unwrap();
        let min_ims = dia_pasef_window
            .scan_ends
            .last()
            .map(|s| dummy_ims_converter(*s as f64))
            .unwrap();

        // Since this layout is in essence a trapezoid, it gets summarized to 5 points
        // (the simplification removes intermediate points)
        // The first two edges are the left side (mz 90-110 at max ims to min ims)
        // The next two edges are the right side (mz 85-105 at min ims to max ims)
        // (note: this connects the last point of the first edge to the first point of the second edge)
        // And the last point closes the polygon back to the start
        let expect = [
            (90.0, max_ims),
            (85.0, min_ims),
            (105.0, min_ims),
            (110.0, max_ims),
            (90.0, max_ims),
        ];

        for i in 0..expect.len() {
            let (exp_mz, exp_im) = expect[i];
            let (got_mz, got_im) = coords[i];
            assert!(
                (exp_mz - got_mz).abs() < 0.001,
                "Expected mz {:.2}, got {:.2}",
                exp_mz,
                got_mz
            );
            assert!(
                (exp_im - got_im).abs() < 0.001,
                "Expected im {:.4}, got {:.4}",
                exp_im,
                got_im
            );
        }
    }

    #[test]
    fn test_intresection_diapasef() {
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
        let quad_polys =
            QuadrupoleIsolationScheme::from_quad(&dia_pasef_window, dummy_ims_converter);
        // This range should intersect with the first window
        // mz: 880-920, im: 1.0-1.4...ish
        assert!(quad_polys.intersects((880.0, 920.0), (1.0, 1.4)));
    }
}
