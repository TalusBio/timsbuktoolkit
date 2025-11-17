// ANSI color codes
const COLOR_GRAY: &str = "\x1b[90m";
const COLOR_BLUE: &str = "\x1b[94m";
const COLOR_CYAN: &str = "\x1b[96m";
const COLOR_RED: &str = "\x1b[91m";
const COLOR_RESET: &str = "\x1b[0m";

use crate::Grid;

impl Grid {
    /// Displays the grid as a heatmap in the terminal.
    /// Uses grayscale blocks for frequency intensity and red borders for non-suppressed nodes.
    pub fn display_heatmap(&self) {
        println!("\n{}", self.format_heatmap());
    }

    /// Displays the grid with custom downscaling parameters.
    ///
    /// # Arguments
    /// * `max_size` - Maximum dimension (width or height) for the output
    pub fn display_heatmap_sized(&self, max_sizes: (usize, usize)) {
        println!("\n{}", self.format_heatmap_sized(max_sizes));
    }

    /// Formats the grid as a string heatmap with automatic downscaling.
    /// Downscales by powers of 2 until both dimensions fit within max_size.
    pub fn format_heatmap_sized(&self, max_sizes: (usize, usize)) -> String {
        // Calculate required downscale factor
        let scale_factor1 = Self::calculate_scale_factor(self.bins, max_sizes.0);
        let scale_factor2 = Self::calculate_scale_factor(self.bins, max_sizes.1);
        self.format_heatmap_with_scale(scale_factor1, scale_factor2)
    }

    /// Formats the grid as a string heatmap for logging or display.
    /// Uses a default Y-axis scale of 2 to account for character aspect ratio.
    pub fn format_heatmap(&self) -> String {
        self.format_heatmap_sized((60, 30))
    }

    /// Calculates the scale factor needed to fit dimensions within max_size.
    /// Returns the smallest power of 2 that brings the size under max_size.
    fn calculate_scale_factor(bins: usize, max_size: usize) -> usize {
        let mut scale = 1;
        while bins / scale > max_size {
            scale *= 2;
        }
        scale
    }

    /// Formats the grid as a string heatmap with custom downscaling on both axes.
    ///
    /// # Arguments
    /// * `x_scale` - Factor to downsample X-axis (e.g., 2 = half the columns)
    /// * `y_scale` - Factor to downsample Y-axis (e.g., 2 = half the rows)
    pub fn format_heatmap_with_scale(&self, x_scale: usize, y_scale: usize) -> String {
        let mut output = String::new();

        // Find max frequency for normalization
        let max_freq = self
            .nodes
            .iter()
            .map(|n| n.center.weight)
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(1.0);

        // Find min frequency for legend
        let min_freq = self
            .nodes
            .iter()
            .map(|n| n.center.weight)
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(0.0);

        // Calculate subsampled dimensions
        let display_cols = self.bins.div_ceil(x_scale);
        let display_rows = self.bins.div_ceil(y_scale);

        // Top border
        output.push('╔');
        output.push_str(&"═".repeat(display_cols));
        output.push_str("╗\n");

        // Grid rows (subsampled)
        for display_r in 0..display_rows {
            output.push('║');
            for display_c in 0..display_cols {
                // Max pooling over x_scale × y_scale block
                let mut max_node = None;
                let mut max_freq_in_block = 0.0;

                for y_offset in 0..y_scale {
                    let r = display_r * y_scale + y_offset;
                    if r >= self.bins {
                        break;
                    }

                    for x_offset in 0..x_scale {
                        let c = display_c * x_scale + x_offset;
                        if c >= self.bins {
                            break;
                        }

                        let idx = r * self.bins + c;
                        let node = &self.nodes[idx];

                        if node.center.weight > max_freq_in_block {
                            max_freq_in_block = node.center.weight;
                            max_node = Some(node);
                        }
                    }
                }

                // Use the max node from the block (or default if empty)
                let node = max_node.unwrap_or(&self.nodes[0]);

                // Normalize center.weight to 0-1 range
                let intensity = if max_freq > 0.0 {
                    node.center.weight / max_freq
                } else {
                    0.0
                };

                // Choose block character based on intensity
                let block = get_block_char(intensity);

                // Color: red for non-suppressed, gray for suppressed
                let colored = if !node.suppressed {
                    format!("{}{}{}", COLOR_RED, block, COLOR_RESET)
                } else {
                    format!("{}{}{}", COLOR_GRAY, block, COLOR_RESET)
                };

                output.push_str(&colored);
            }
            output.push_str("║\n");
        }

        // Bottom border
        output.push('╚');
        output.push_str(&"═".repeat(display_cols));
        output.push_str("╝\n");

        // Legend
        let scale_info = if x_scale > 1 || y_scale > 1 {
            format!(" (Scale: {}x × {}y)", x_scale, y_scale)
        } else {
            String::new()
        };

        output.push_str(&format!(
            "\n  Legend: {}█{} = suppressed, {}█{} = non-suppressed (max), Range: {}-{} hits{}\n",
            COLOR_GRAY, COLOR_RESET, COLOR_RED, COLOR_RESET, min_freq, max_freq, scale_info
        ));

        output
    }
}

/// Maps intensity (0.0 to 1.0) to Unicode block characters
fn get_block_char(intensity: f64) -> &'static str {
    match intensity {
        i if i >= 0.875 => "█", // Full block
        i if i >= 0.750 => "▓", // Dark shade
        i if i >= 0.625 => "▒", // Medium shade
        i if i >= 0.500 => "░", // Light shade
        i if i >= 0.375 => "▒", // Medium shade
        i if i >= 0.250 => "░", // Light shade
        i if i >= 0.125 => "·", // Dot
        _ => " ",               // Empty
    }
}

/// Plots a function in the terminal with customizable dimensions.
///
/// # Arguments
/// * `f` - The function to plot (takes f64, returns f64)
/// * `x_range` - Tuple of (min, max) for x-axis
/// * `width` - Number of columns for the plot
/// * `height` - Number of rows for the plot
///
/// # Example
/// ```
/// use calibrt::plotting::plot_function;
/// plot_function(|x| x.sin(), (-3.14, 3.14), 60, 20);
/// plot_function(|x| x * x, (-5.0, 5.0), 80, 30);
/// ```
pub fn plot_function<F>(f: F, x_range: (f64, f64), width: usize, height: usize)
where
    F: Fn(f64) -> Result<f64, f64>,
{
    println!("\n{}", format_function_plot(f, x_range, width, height));
}

/// Formats a function plot as a string for logging or display.
pub fn format_function_plot<F>(f: F, x_range: (f64, f64), width: usize, height: usize) -> String
where
    F: Fn(f64) -> Result<f64, f64>,
{
    let mut output = String::new();

    let (x_min, x_max) = x_range;
    let x_span = x_max - x_min;

    // Sample the function at each x position
    let mut samples = Vec::with_capacity(width);
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;

    for i in 0..width {
        let x = x_min + (i as f64 / (width - 1) as f64) * x_span;
        let (is_err, y) = match f(x) {
            Ok(y) => (false, y),
            Err(y) => (true, y),
        };

        // Track min/max for y-axis scaling
        if y.is_finite() {
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }

        samples.push((x, y, is_err));
    }

    // Handle edge cases
    if !y_min.is_finite() || !y_max.is_finite() || y_min == y_max {
        y_min = -1.0;
        y_max = 1.0;
    }

    let y_span = y_max - y_min;

    // Create the plot grid
    let mut grid = vec![vec![(' ', false); width]; height];

    // Plot the function
    for (col, &(x, y, err)) in samples.iter().enumerate() {
        if y.is_finite() {
            // Map y to row (inverted because row 0 is at top)
            let normalized = (y - y_min) / y_span;
            let row = ((1.0 - normalized) * (height - 1) as f64) as usize;
            let row = row.min(height - 1);

            grid[row][col] = ('●', err);
        }
    }

    // Add axes if they're in range
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

    // Draw axes
    if let Some(row) = zero_row {
        for col in 0..width {
            if grid[row][col].0 == ' ' {
                grid[row][col].0 = '─';
            }
        }
    }

    if let Some(col) = zero_col {
        for row in grid.iter_mut().take(height) {
            if row[col].0 == ' ' {
                row[col].0 = '│';
            } else if row[col].0 == '─' {
                row[col].0 = '┼';
            }
        }
    }

    // Top border
    output.push('╔');
    output.push_str(&"═".repeat(width));
    output.push_str("╗\n");

    // Render grid with colors
    for row in &grid {
        output.push('║');
        for &(ch, err) in row {
            let colored = match ch {
                '●' => format!(
                    "{}{}{}",
                    if err { COLOR_BLUE } else { COLOR_CYAN },
                    ch,
                    COLOR_RESET
                ),
                '─' | '│' | '┼' => format!("{}{}{}", COLOR_GRAY, ch, COLOR_RESET),
                _ => ch.to_string(),
            };
            output.push_str(&colored);
        }
        output.push_str("║\n");
    }

    // Bottom border
    output.push('╚');
    output.push_str(&"═".repeat(width));
    output.push_str("╝\n");

    // Legend with ranges
    output.push_str(&format!(
        "\n  X: [{:.2}, {:.2}]  Y: [{:.2}, {:.2}]  Size: {}×{}\n",
        x_min, x_max, y_min, y_max, width, height
    ));

    output
}
