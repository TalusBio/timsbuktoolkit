#[derive(Debug)]
pub enum CliError {
    Config {
        source: String,
    },
    ParseError {
        msg: String,
    },
    Io {
        source: String,
        path: Option<String>,
    },
}
