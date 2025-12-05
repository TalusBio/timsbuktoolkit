use crate::Grid;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Color {
    Gray,
    Blue,
    Cyan,
    Red,
    Reset,
}

impl Color {
    fn to_ansi(self) -> &'static str {
        match self {
            Color::Gray => "\x1b[90m",
            Color::Blue => "\x1b[94m",
            Color::Cyan => "\x1b[96m",
            Color::Red => "\x1b[91m",
            Color::Reset => "\x1b[0m",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Cell {
    char: &'static str,
    color: Color,
}

impl Cell {
    const fn new(char: &'static str, color: Color) -> Self {
        Self { char, color }
    }
}

impl Default for Cell {
    fn default() -> Self {
        Self::new(" ", Color::Reset)
    }
}

struct Canvas {
    grid: Vec<Vec<Cell>>,
    width: usize,
    height: usize,
}

impl Canvas {
    fn new(width: usize, height: usize) -> Self {
        Self {
            grid: vec![vec![Cell::default(); width]; height],
            width,
            height,
        }
    }

    fn set(&mut self, x: usize, y: usize, cell: Cell) {
        if x < self.width && y < self.height {
            self.grid[y][x] = cell;
        }
    }

    fn render(&self, legend: &str) -> String {
        let mut output = String::new();
        output.push('╔');
        output.push_str(&"═".repeat(self.width));
        output.push_str("╗\n");
        for row in &self.grid {
            output.push('║');
            for &cell in row {
                let colored = if cell.color != Color::Reset {
                    format!(
                        "{}{}{}",
                        cell.color.to_ansi(),
                        cell.char,
                        Color::Reset.to_ansi()
                    )
                } else {
                    cell.char.to_string()
                };
                output.push_str(&colored);
            }
            output.push_str("║\n");
        }
        output.push('╚');
        output.push_str(&"═".repeat(self.width));
        output.push_str("╝\n");
        output.push_str(legend);
        output
    }
}

impl Grid {
    pub fn display_heatmap(&self) {
        println!("\n{}", self.format_heatmap());
    }

    pub fn display_heatmap_sized(&self, max_sizes: (usize, usize)) {
        println!("\n{}", self.format_heatmap_sized(max_sizes));
    }

    pub fn format_heatmap_sized(&self, max_sizes: (usize, usize)) -> String {
        let scale_factor1 = Self::calculate_scale_factor(self.bins, max_sizes.0);
        let scale_factor2 = Self::calculate_scale_factor(self.bins, max_sizes.1);
        self.format_heatmap_with_scale(scale_factor1, scale_factor2)
    }

    pub fn format_heatmap(&self) -> String {
        self.format_heatmap_sized((60, 30))
    }

    fn calculate_scale_factor(bins: usize, max_size: usize) -> usize {
        let mut scale = 1;
        while bins.div_ceil(scale) > max_size {
            scale *= 2;
        }
        scale
    }

    pub fn format_heatmap_with_scale(&self, x_scale: usize, y_scale: usize) -> String {
        let (canvas, legend) = self.rasterize_heatmap(x_scale, y_scale);
        canvas.render(&legend)
    }

    fn rasterize_heatmap(&self, x_scale: usize, y_scale: usize) -> (Canvas, String) {
        let max_freq = self
            .nodes
            .iter()
            .map(|n| n.center.weight)
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(1.0);

        let min_freq = self
            .nodes
            .iter()
            .map(|n| n.center.weight)
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(0.0);

        let display_cols = self.bins.div_ceil(x_scale);
        let display_rows = self.bins.div_ceil(y_scale);
        let mut canvas = Canvas::new(display_cols, display_rows);

        for display_row in 0..display_rows {
            for display_col in 0..display_cols {
                let mut max_node = None;
                let mut max_freq_in_block = 0.0;

                for row_offset in 0..y_scale {
                    let grid_row = display_row * y_scale + row_offset;
                    if grid_row >= self.bins {
                        break;
                    }
                    for col_offset in 0..x_scale {
                        let grid_col = display_col * x_scale + col_offset;
                        if grid_col >= self.bins {
                            break;
                        }
                        let idx = grid_row * self.bins + grid_col;
                        let node = &self.nodes[idx];
                        if node.center.weight > max_freq_in_block {
                            max_freq_in_block = node.center.weight;
                            max_node = Some(node);
                        }
                    }
                }

                // Determine color based on whether there's a non-suppressed node in the block
                // If max_node is None (all nodes have weight 0), treat as suppressed
                let (intensity, color) = if let Some(node) = max_node {
                    let intensity = if max_freq > 0.0 {
                        node.center.weight / max_freq
                    } else {
                        0.0
                    };
                    let color = if !node.suppressed {
                        Color::Red
                    } else {
                        Color::Gray
                    };
                    (intensity, color)
                } else {
                    // No nodes with weight > 0 in this block, treat as suppressed
                    (0.0, Color::Gray)
                };
                let block = get_block_char(intensity);
                canvas.set(display_col, display_row, Cell::new(block, color));
            }
        }

        let scale_info = if x_scale > 1 || y_scale > 1 {
            format!(" (Scale: {}x × {}y)", x_scale, y_scale)
        } else {
            String::new()
        };
        let legend = format!(
            "\n  Legend: {}█{} = suppressed, {}█{} = non-suppressed (max), Range: {}-{} hits{}\n",
            Color::Gray.to_ansi(),
            Color::Reset.to_ansi(),
            Color::Red.to_ansi(),
            Color::Reset.to_ansi(),
            min_freq,
            max_freq,
            scale_info
        );

        (canvas, legend)
    }
}

fn get_block_char(intensity: f64) -> &'static str {
    match intensity {
        i if i >= 0.875 => "█",
        i if i >= 0.750 => "▓",
        i if i >= 0.625 => "▒",
        i if i >= 0.500 => "░",
        i if i >= 0.375 => "▒",
        i if i >= 0.250 => "░",
        i if i >= 0.125 => "·",
        _ => " ",
    }
}

pub fn plot_function<F>(f: F, x_range: (f64, f64), width: usize, height: usize)
where
    F: Fn(f64) -> Result<f64, f64>,
{
    println!("\n{}", format_function_plot(f, x_range, width, height));
}

pub fn format_function_plot<F>(f: F, x_range: (f64, f64), width: usize, height: usize) -> String
where
    F: Fn(f64) -> Result<f64, f64>,
{
    let (canvas, legend) = rasterize_function_plot(f, x_range, width, height);
    canvas.render(&legend)
}

fn rasterize_function_plot<F>(
    f: F,
    x_range: (f64, f64),
    width: usize,
    height: usize,
) -> (Canvas, String)
where
    F: Fn(f64) -> Result<f64, f64>,
{
    let (x_min, x_max) = x_range;
    let x_span = x_max - x_min;

    let mut samples = Vec::with_capacity(width);
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;

    for i in 0..width {
        let x = x_min + (i as f64 / (width - 1) as f64) * x_span;
        let (is_err, y) = match f(x) {
            Ok(y) => (false, y),
            Err(y) => (true, y),
        };
        if y.is_finite() {
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }
        samples.push((x, y, is_err));
    }

    if !y_min.is_finite() || !y_max.is_finite() || y_min == y_max {
        y_min = -1.0;
        y_max = 1.0;
    }
    let y_span = y_max - y_min;

    let mut canvas = Canvas::new(width, height);

    let zero_row = if y_min <= 0.0 && y_max >= 0.0 {
        let normalized = (0.0 - y_min) / y_span;
        Some(((1.0 - normalized) * (height - 1) as f64) as usize)
    } else {
        None
    };

    let zero_col = if x_min <= 0.0 && x_max >= 0.0 {
        let normalized = (0.0 - x_min) / x_span;
        Some((normalized * (width - 1) as f64) as usize)
    } else {
        None
    };

    if let Some(row) = zero_row {
        for col in 0..width {
            canvas.set(col, row, Cell::new("─", Color::Gray));
        }
    }

    if let Some(col) = zero_col {
        for row in 0..height {
            let ch = if Some(row) == zero_row { "┼" } else { "│" };
            canvas.set(col, row, Cell::new(ch, Color::Gray));
        }
    }

    for (col, &(_x, y, err)) in samples.iter().enumerate() {
        if y.is_finite() {
            let normalized = (y - y_min) / y_span;
            let row = ((1.0 - normalized) * (height - 1) as f64) as usize;
            let row = row.min(height - 1);
            let color = if err { Color::Blue } else { Color::Cyan };
            canvas.set(col, row, Cell::new("●", color));
        }
    }

    let legend = format!(
        "\n  X: [{:.2}, {:.2}]  Y: [{:.2}, {:.2}]  Size: {}×{}\n",
        x_min, x_max, y_min, y_max, width, height
    );

    (canvas, legend)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Grid;

    #[test]
    fn test_plot_sine_wave() {
        let plot = format_function_plot(|x| Ok(x.sin()), (-3.14, 3.14), 60, 20);
        insta::assert_snapshot!(plot);
    }

    #[test]
    fn test_plot_quadratic_function() {
        let plot = format_function_plot(|x| Ok(x * x), (-5.0, 5.0), 80, 30);
        insta::assert_snapshot!(plot);
    }

    #[test]
    fn test_plot_with_errors() {
        let plot = format_function_plot(
            |x| {
                if x > 0.0 { Err(x.cos()) } else { Ok(x.sin()) }
            },
            (-3.14, 3.14),
            60,
            20,
        );
        insta::assert_snapshot!(plot);
    }

    fn create_test_grid(bins: usize) -> Grid {
        let mut grid = Grid::new(bins, (0.0, 1.0), (0.0, 1.0)).unwrap();
        for (i, node) in grid.nodes.iter_mut().enumerate() {
            node.center.weight = (i % 25) as f64;
            node.suppressed = i % 7 == 0;
        }
        grid
    }

    #[test]
    fn test_heatmap_default_size() {
        let grid = create_test_grid(100);
        let heatmap = grid.format_heatmap();
        insta::assert_snapshot!(heatmap);
    }

    #[test]
    fn test_heatmap_custom_size() {
        let grid = create_test_grid(100);
        let heatmap = grid.format_heatmap_sized((80, 40));
        insta::assert_snapshot!(heatmap);
    }

    #[test]
    fn test_heatmap_no_downscaling() {
        let grid = create_test_grid(20);
        let heatmap = grid.format_heatmap_sized((30, 30));
        insta::assert_snapshot!(heatmap);
    }

    #[test]
    fn test_heatmap_subsampling_empty_blocks() {
        // This test verifies that empty blocks (weight=0) are shown as suppressed (gray)
        // and don't inherit the suppression state of nodes[0]
        let mut grid = Grid::new(10, (0.0, 10.0), (0.0, 10.0)).unwrap();

        // Set nodes[0] as non-suppressed with some weight
        grid.nodes[0].center.weight = 5.0;
        grid.nodes[0].suppressed = false;

        // Set a few other nodes with weight, all suppressed
        grid.nodes[50].center.weight = 10.0;
        grid.nodes[50].suppressed = true;

        // All other nodes have weight 0 and are suppressed
        for i in 1..grid.nodes.len() {
            if i != 50 {
                grid.nodes[i].center.weight = 0.0;
                grid.nodes[i].suppressed = true;
            }
        }

        // Render with downsampling (10x10 grid into 5x5 display)
        let heatmap = grid.format_heatmap_sized((5, 5));

        // Count red (non-suppressed) cells in the output
        // Before the fix: many cells will be red due to inheriting nodes[0].suppressed
        // After the fix: only the cell containing nodes[0] should be red
        let red_ansi = Color::Red.to_ansi();
        let red_count = heatmap.matches(red_ansi).count();

        println!("Red cells found: {}", red_count);

        // Before fix: red_count will be much higher (20+ due to many empty blocks being red)
        // After fix: red_count should be 1-2 (only the block containing nodes[0])
        assert!(
            red_count <= 2,
            "Found {} red cells, expected ≤2. Empty blocks are incorrectly inheriting nodes[0].suppressed",
            red_count
        );

        // Snapshot should show only one red block (containing nodes[0])
        // and all empty blocks as gray
        insta::assert_snapshot!(heatmap);
    }
}
