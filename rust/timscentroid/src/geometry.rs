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
        let geom = quad_to_geometry(quad, ims_converter);
        Self { inner: geom }
    }
}

fn connect_edges(left: &[(f64, f64)], right: &[(f64, f64)]) -> Polygon<f64> {
    let mut all_points = left.to_vec();
    all_points.extend(right.iter().rev());
    Polygon::new(all_points.into(), vec![])
}

// We let geo-rs do the heavy lifting of geometry operations
// so we just make possibly a lot of boxes and then let it simplify
// the geometry for us.
fn quad_to_geometry(
    quad: &QuadrupoleSettings,
    ims_converter: impl Fn(f64) -> f64,
) -> MultiPolygon<f64> {
    assert!(quad.scan_starts.len() == quad.scan_ends.len());
    assert!(quad.scan_starts.len() == quad.isolation_mz.len());
    assert!(quad.scan_starts.len() == quad.isolation_width.len());
    assert!(quad.scan_starts.len() == quad.collision_energy.len());
    if quad.scan_starts.is_empty() {
        return MultiPolygon(vec![]);
    }
    let mut polygons = Vec::new();
    let mut curr_left_edge = Vec::new();
    let mut curr_right_edge = Vec::new();

    let mut last_scan_end: Option<usize> = None;
    let mut last_quad_range: Option<(f64, f64)> = None;

    for i in 0..quad.scan_starts.len() {
        let scan_start = quad.scan_starts[i];
        let scan_end = quad.scan_ends[i];
        let im_start = ims_converter(scan_start as f64);
        let im_end = ims_converter(scan_end as f64);
        let mz_center = quad.isolation_mz[i];
        let mz_width = quad.isolation_width[i];
        let mz_start = mz_center - mz_width / 2.0;
        let mz_end = mz_center + mz_width / 2.0;

        let mut flush_polygon = false;

        if let Some((last_mz_start, last_mz_end)) = last_quad_range {
            // If the current mz range does not overlap with the last one,
            // we need to flush the current polygon.
            if mz_start > last_mz_end || mz_end < last_mz_start {
                flush_polygon = true;
            }
        }

        if let Some(last_end) = last_scan_end {
            // If the current im_start is less than the last_im_end,
            // we have an overlap in the IMS dimension.
            // We also need the quads to overlap in the mz dimension
            if scan_start != last_end {
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
        last_scan_end = Some(quad.scan_ends[i]);
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
        let geom = quad_to_geometry(&dia_pasef_window, dummy_ims_converter);
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
        let geom = quad_to_geometry(&dia_pasef_window, dummy_ims_converter);
        println!("Geometry: {:?}", geom);
        assert_eq!(geom.0.len(), 1);
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
