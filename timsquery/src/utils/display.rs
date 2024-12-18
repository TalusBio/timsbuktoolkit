#[derive(Debug, Clone, Copy)]
pub struct GlimpseConfig {
    pub max_items: usize,
    pub padding: usize,
    pub new_line: bool,
}

impl Default for GlimpseConfig {
    fn default() -> Self {
        GlimpseConfig {
            max_items: 10,
            padding: 0,
            new_line: false,
        }
    }
}

pub fn glimpse_vec<T: std::fmt::Debug>(v: &[T], config: Option<GlimpseConfig>) -> String {
    let config = config.unwrap_or_default();
    let len = v.len();
    let mut out = String::new();

    let separator = if config.new_line { ",\n" } else { ", " };
    let padding = " ".repeat(config.padding);

    if len > config.max_items {
        out.push_str("[\n");
        v[..3]
            .iter()
            .map(|x| format!("{}{:?}{}", padding, x, separator))
            .for_each(|x| out.push_str(&x));
        out.push_str(&format!("{}...\n", padding));
        v[len - 3..]
            .iter()
            .map(|x| format!("{}{:?}{}", padding, x, separator))
            .for_each(|x| out.push_str(&x));
        out.push_str(&format!("{}] len = {}", padding, len));
        if config.new_line {
            out.push('\n');
        }
    } else {
        out.push_str("[\n");
        v.iter()
            .map(|x| format!("{}{:?}{}", padding, x, separator))
            .for_each(|x| out.push_str(&x));
        out.push(']');
    }

    // Remove trailing separator
    if config.new_line {
        out = out.trim_end_matches(",\n").to_string();
    } else {
        out = out.trim_end_matches(", ").to_string();
    }

    out
}
