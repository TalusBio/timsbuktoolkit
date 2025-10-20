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
                let center_x = x_range.0 + (c as f64 + 0.5) * (x_span / bins as f64);
                let center_y = y_range.0 + (r as f64 + 0.5) * (y_span / bins as f64);
                nodes.push(Node {
                    center: Point {
                        x: center_x,
                        y: center_y,
                        weight: 0.0,
                    },
                    suppressed: false,
                });
            }
        }

        Ok(Self {
            nodes,
            x_range,
            y_range,
            x_span,
            y_span,
            bins,
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
        let Point { x, y, weight } = point;

        // If the weight is infinite or NaN, we yell ...
        if weight.is_infinite() || weight.is_nan() {
            return Err(CalibRtError::UnsupportedWeight(*weight));
        }

        let gx = (((x - self.x_range.0) / self.x_span) * self.bins as f64) as usize;
        let gy = (((y - self.y_range.0) / self.y_span) * self.bins as f64) as usize;

        let gx = gx.min(self.bins - 1);
        let gy = gy.min(self.bins - 1);

        let index = gy * self.bins + gx;
        if let Some(node) = self.nodes.get_mut(index) {
            node.center.weight += weight;
        }

        Ok(())
    }

    /// Applies nonmaximal suppression to the grid nodes.
    pub fn suppress_nonmax(&mut self) -> Result<(), CalibRtError> {
        // We start with 1s to prevent the max being empty in all ...
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
        self.nodes.iter().for_each(|x| {
            if x.suppressed {
                suppressed_sum += x.center.weight;
            } else {
                non_suppressed_sum += x.center.weight;
            }
        });
        info!(
            "Suppression complete. Suppressed weight sum: {}, Non-suppressed weight sum: {}",
            suppressed_sum, non_suppressed_sum
        );
        if non_suppressed_sum == 0.0 {
            return Err(CalibRtError::NoPoints);
        }
        Ok(())
    }
}

/// Represents a node (cell) in the grid.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Node {
    pub(crate) center: Point,
    pub(crate) suppressed: bool,
}
