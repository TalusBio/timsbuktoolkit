use array2d::Array2D;
use crate::{
    CalibRtError,
    Point,
};
use tracing::info;

pub struct Grid {
    pub(crate) nodes: Vec<Node>,
    pub(crate) x_range: (f64, f64),
    pub(crate) y_range: (f64, f64),
    x_span: f64,
    y_span: f64,
    pub(crate) bins: usize,
    /// Ping-pong weight buffers. A = raw weights, B = blurred output.
    weights_a: Array2D<f64>,
    weights_b: Array2D<f64>,
}

impl Grid {
    /// Creates a new, empty grid with a fixed geometry.
    /// The center of each node is constant based on the grid resolution.
    pub fn new(
        bins: usize,
        x_range: (f64, f64),
        y_range: (f64, f64),
    ) -> Result<Self, CalibRtError> {
        if bins == 0 {
            return Err(CalibRtError::ZeroRange);
        };
        let x_span = x_range.1 - x_range.0;
        let y_span = y_range.1 - y_range.0;

        if x_span <= 0.0 || y_span <= 0.0 {
            return Err(CalibRtError::ZeroRange);
        }

        let mut nodes = Vec::with_capacity(bins * bins);
        for r in 0..bins {
            for c in 0..bins {
                // Add 0.5 to place node center at the midpoint of each bin
                let center_x = x_range.0 + (c as f64 + 0.5) * (x_span / bins as f64);
                let center_y = y_range.0 + (r as f64 + 0.5) * (y_span / bins as f64);
                nodes.push(Node {
                    center: Point {
                        library: center_x,
                        observed: center_y,
                        weight: 0.0,
                    },
                    suppressed: false,
                    sum_wx: 0.0,
                    sum_wy: 0.0,
                    sum_w: 0.0,
                });
            }
        }

        let weights_a = Array2D::from_flat_vector(vec![0.0; bins * bins], bins, bins)
            .expect("Grid dimensions are valid");
        let weights_b = Array2D::from_flat_vector(vec![0.0; bins * bins], bins, bins)
            .expect("Grid dimensions are valid");

        Ok(Self {
            nodes,
            x_range,
            y_range,
            x_span,
            y_span,
            bins,
            weights_a,
            weights_b,
        })
    }

    pub fn extend_points<'a, T>(&mut self, points: T) -> Result<(), CalibRtError>
    where
        T: IntoIterator<Item = &'a Point> + 'a,
    {
        points.into_iter().try_for_each(|p| self.add_point(p))
    }

    /// Adds a single point to the grid, incrementing the frequency of the corresponding cell.
    pub fn add_point(&mut self, point: &Point) -> Result<(), CalibRtError> {
        let Point { library, observed, weight } = point;

        if !weight.is_finite() || !library.is_finite() || !observed.is_finite() {
            return Err(CalibRtError::UnsupportedWeight(*weight));
        }

        let gx = (((library - self.x_range.0) / self.x_span) * self.bins as f64) as usize;
        let gy = (((observed - self.y_range.0) / self.y_span) * self.bins as f64) as usize;

        let gx = gx.min(self.bins - 1);
        let gy = gy.min(self.bins - 1);

        let index = gy * self.bins + gx;
        if let Some(node) = self.nodes.get_mut(index) {
            node.center.weight += weight;
            node.sum_wx += library * weight;
            node.sum_wy += observed * weight;
            node.sum_w += weight;
        }

        Ok(())
    }

    /// Applies nonmaximal suppression to the grid nodes.
    ///
    /// A node is marked as non-suppressed only if it is the maximum weight
    /// in BOTH its row AND its column. This ensures we keep only the most
    /// significant alignment points in each dimension.
    ///
    /// # Returns
    /// - `Ok(())` if at least one node remains non-suppressed
    /// - `Err(CalibRtError::NoPoints)` if all nodes have zero weight
    pub fn suppress_nonmax(&mut self) -> Result<(), CalibRtError> {
        // Initialize with 1.0 to handle empty grids gracefully
        let mut max_in_row = vec![1.; self.bins];
        let mut max_in_col = vec![1.; self.bins];

        for (r, mrow_elem) in max_in_row.iter_mut().enumerate() {
            for (c, mcol_elem) in max_in_col.iter_mut().enumerate() {
                let index = r * self.bins + c;
                let weight = self.nodes[index].center.weight;
                if &weight > mrow_elem {
                    *mrow_elem = weight;
                }
                if &weight > mcol_elem {
                    *mcol_elem = weight;
                }
            }
        }

        for (index, node) in self.nodes.iter_mut().enumerate() {
            let r = index / self.bins;
            let c = index % self.bins;
            node.suppressed = true;
            if node.center.weight == max_in_row[r] && node.center.weight == max_in_col[c] {
                node.suppressed = false;
            }
        }

        let mut suppressed_sum = 0.0;
        let mut non_suppressed_sum = 0.0;
        let mut num_unsuppressed = 0;
        self.nodes.iter().for_each(|x| {
            if x.suppressed {
                suppressed_sum += x.center.weight;
            } else {
                non_suppressed_sum += x.center.weight;
                num_unsuppressed += 1;
            }
        });
        info!(
            "Suppression complete. Suppressed weight sum: {}, Non-suppressed weight sum: {}, Num Unsuppressed: {}",
            suppressed_sum, non_suppressed_sum, num_unsuppressed
        );
        if non_suppressed_sum == 0.0 {
            return Err(CalibRtError::NoPoints);
        }

        // Replace bin centers with weighted centroids for non-suppressed nodes.
        // Clamp to grid range to avoid interpolation boundary issues.
        for node in self.nodes.iter_mut() {
            if !node.suppressed && node.sum_w > 0.0 {
                let cx = (node.sum_wx / node.sum_w).clamp(self.x_range.0, self.x_range.1);
                let cy = (node.sum_wy / node.sum_w).clamp(self.y_range.0, self.y_range.1);
                node.center.library = cx;
                node.center.observed = cy;
            }
        }

        Ok(())
    }

    /// Zero all node weights and suppression flags, preserving bin geometry.
    /// Restores each node center to the midpoint of its bin.
    pub fn reset(&mut self) {
        for (i, node) in self.nodes.iter_mut().enumerate() {
            let r = i / self.bins;
            let c = i % self.bins;
            node.center = Point {
                library: self.x_range.0 + (c as f64 + 0.5) * (self.x_span / self.bins as f64),
                observed: self.y_range.0 + (r as f64 + 0.5) * (self.y_span / self.bins as f64),
                weight: 0.0,
            };
            node.suppressed = false;
            node.sum_wx = 0.0;
            node.sum_wy = 0.0;
            node.sum_w = 0.0;
        }
        self.weights_a.reset_with_value(self.bins, self.bins, 0.0);
        self.weights_b.reset_with_value(self.bins, self.bins, 0.0);
    }

    /// Reset the grid with new dimensions and ranges. Reallocates if the
    /// bin count changes.
    pub fn reconfigure(
        &mut self,
        bins: usize,
        x_range: (f64, f64),
        y_range: (f64, f64),
    ) -> Result<(), CalibRtError> {
        *self = Self::new(bins, x_range, y_range)?;
        Ok(())
    }

    /// Read access to all grid cells.
    pub fn grid_cells(&self) -> &[Node] {
        &self.nodes
    }

    /// Copy node weights into buffer A for blur processing.
    pub(crate) fn sync_weights(&mut self) {
        for (i, node) in self.nodes.iter().enumerate() {
            let r = i / self.bins;
            let c = i % self.bins;
            self.weights_a.insert(r, c, node.center.weight);
        }
    }

    /// Apply 3x3 approximate gaussian blur: A -> B.
    pub(crate) fn blur_weights(&mut self) {
        const KERNEL: [f64; 9] = [
            1.0/10.0, 1.0/10.0, 1.0/10.0,
            1.0/10.0, 2.0/10.0, 1.0/10.0,
            1.0/10.0, 1.0/10.0, 1.0/10.0,
        ];
        self.weights_a.convolve_2d_into(&KERNEL, &mut self.weights_b);
    }

    /// Read access to blurred weight at grid position (row, col).
    pub(crate) fn blurred_weight(&self, row: usize, col: usize) -> f64 {
        self.weights_b.as_flat_slice()[row * self.bins + col]
    }
}

/// Represents a node (cell) in the grid.
#[derive(Debug, Clone, Copy, Default)]
pub struct Node {
    pub center: Point,
    pub suppressed: bool,
    // Internal accumulators — not exposed beyond crate
    pub(crate) sum_wx: f64,
    pub(crate) sum_wy: f64,
    pub(crate) sum_w: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point;

    #[test]
    fn test_grid_reset_preserves_allocation() {
        let mut grid = Grid::new(10, (0.0, 100.0), (0.0, 100.0)).unwrap();
        grid.add_point(&Point { library: 50.0, observed: 50.0, weight: 1.0 }).unwrap();

        let capacity_before = grid.nodes.capacity();
        assert!(grid.grid_cells().iter().any(|n| n.sum_w > 0.0));

        grid.reset();

        assert_eq!(grid.nodes.capacity(), capacity_before);
        assert!(grid.grid_cells().iter().all(|n| n.sum_w == 0.0));
        assert!(grid.grid_cells().iter().all(|n| !n.suppressed));
    }

    /// Helper function to print grid state and return non-suppressed nodes
    fn print_grid_state(grid: &Grid) -> Vec<(usize, usize, f64)> {
        println!("\nGrid state (S=suppressed, N=not suppressed):");
        for r in 0..grid.bins {
            for c in 0..grid.bins {
                let idx = r * grid.bins + c;
                let marker = if grid.nodes[idx].suppressed { "S" } else { "N" };
                print!("{:4.0}{} ", grid.nodes[idx].center.weight, marker);
            }
            println!();
        }

        let non_suppressed: Vec<_> = grid
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| !n.suppressed)
            .map(|(i, n)| {
                let r = i / grid.bins;
                let c = i % grid.bins;
                (r, c, n.center.weight)
            })
            .collect();

        println!(
            "Non-suppressed nodes (row, col, weight): {:?}",
            non_suppressed
        );
        non_suppressed
    }

    #[test]
    fn test_suppress_nonmax_simple_3x3() {
        // Create a 3x3 grid with known values
        // Grid layout (row, col):
        //   0   1   2
        // 0 [1] [2] [9]  <- max in row 0 is 9
        // 1 [4] [5] [6]  <- max in row 1 is 6
        // 2 [7] [8] [3]  <- max in row 2 is 8
        //   ^   ^   ^
        //   |   |   max in col 2 is 9
        //   |   max in col 1 is 8
        //   max in col 0 is 7

        let mut grid = Grid::new(3, (0.0, 3.0), (0.0, 3.0)).unwrap();

        let test_data = [
            (0.5, 0.5, 1.0), // bin (0,0) = 1
            (1.5, 0.5, 2.0), // bin (0,1) = 2
            (2.5, 0.5, 9.0), // bin (0,2) = 9
            (0.5, 1.5, 4.0), // bin (1,0) = 4
            (1.5, 1.5, 5.0), // bin (1,1) = 5
            (2.5, 1.5, 6.0), // bin (1,2) = 6
            (0.5, 2.5, 7.0), // bin (2,0) = 7
            (1.5, 2.5, 8.0), // bin (2,1) = 8
            (2.5, 2.5, 3.0), // bin (2,2) = 3
        ];

        for (x, y, weight) in test_data.iter() {
            grid.add_point(&Point {
                library: *x,
                observed: *y,
                weight: *weight,
            })
            .unwrap();
        }

        grid.suppress_nonmax().unwrap();

        let non_suppressed = print_grid_state(&grid);

        // Only nodes that are BOTH row max AND column max should be non-suppressed
        // - (0,2) = 9: max_in_row[0]=9 ✓, max_in_col[2]=9 ✓ → NOT suppressed
        // - (2,1) = 8: max_in_row[2]=8 ✓, max_in_col[1]=8 ✓ → NOT suppressed

        assert_eq!(
            non_suppressed.len(),
            2,
            "Expected 2 non-suppressed nodes, found {}",
            non_suppressed.len()
        );

        assert!(non_suppressed.contains(&(0, 2, 9.0)));
        assert!(non_suppressed.contains(&(2, 1, 8.0)));
    }

    #[test]
    fn test_suppress_nonmax_single_global_max() {
        // Create a 3x3 grid where one cell is the max in both its row and column
        let mut grid = Grid::new(3, (0.0, 3.0), (0.0, 3.0)).unwrap();

        let test_data = [
            (0.5, 0.5, 1.0),
            (1.5, 0.5, 2.0),
            (2.5, 0.5, 3.0),
            (0.5, 1.5, 4.0),
            (1.5, 1.5, 9.0),
            (2.5, 1.5, 6.0), // 9 is max
            (0.5, 2.5, 7.0),
            (1.5, 2.5, 8.0),
            (2.5, 2.5, 5.0),
        ];

        for (x, y, weight) in test_data.iter() {
            grid.add_point(&Point {
                library: *x,
                observed: *y,
                weight: *weight,
            })
            .unwrap();
        }

        grid.suppress_nonmax().unwrap();

        let non_suppressed = print_grid_state(&grid);

        // Only the center cell (1,1) with weight 9 should be non-suppressed
        // It is the max in both its row (row 1: 4,9,6) and column (col 1: 2,9,8)
        assert_eq!(
            non_suppressed.len(),
            1,
            "Expected 1 non-suppressed node, found {}",
            non_suppressed.len()
        );

        assert!(non_suppressed.contains(&(1, 1, 9.0)));
    }

    #[test]
    fn test_suppress_nonmax_diagonal_pattern() {
        // Create a diagonal pattern where each diagonal element is max in its row and column
        let mut grid = Grid::new(3, (0.0, 3.0), (0.0, 3.0)).unwrap();

        let test_data = [
            (0.5, 0.5, 9.0),
            (1.5, 0.5, 1.0),
            (2.5, 0.5, 1.0), // (0,0) = 9
            (0.5, 1.5, 1.0),
            (1.5, 1.5, 9.0),
            (2.5, 1.5, 1.0), // (1,1) = 9
            (0.5, 2.5, 1.0),
            (1.5, 2.5, 1.0),
            (2.5, 2.5, 9.0), // (2,2) = 9
        ];

        for (x, y, weight) in test_data.iter() {
            grid.add_point(&Point {
                library: *x,
                observed: *y,
                weight: *weight,
            })
            .unwrap();
        }

        grid.suppress_nonmax().unwrap();

        let non_suppressed = print_grid_state(&grid);

        // All 3 diagonal elements should be non-suppressed
        assert_eq!(
            non_suppressed.len(),
            3,
            "Expected 3 non-suppressed nodes (diagonal), found {}",
            non_suppressed.len()
        );

        assert!(
            non_suppressed.contains(&(0, 0, 9.0)),
            "Diagonal (0,0) should be non-suppressed"
        );
        assert!(
            non_suppressed.contains(&(1, 1, 9.0)),
            "Diagonal (1,1) should be non-suppressed"
        );
        assert!(
            non_suppressed.contains(&(2, 2, 9.0)),
            "Diagonal (2,2) should be non-suppressed"
        );
    }
}
