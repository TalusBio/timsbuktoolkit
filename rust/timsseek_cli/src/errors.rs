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
    DataReading {
        source: String,
    },
}

impl std::fmt::Display for CliError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CliError::Config { source } => write!(f, "Error interpreting the config: {}", source),
            CliError::ParseError { msg } => write!(f, "Error parsing config: {}", msg),
            CliError::Io { source, path } => {
                if let Some(path) = path {
                    write!(f, "Error reading file {}: {}", path, source)
                } else {
                    write!(f, "Error reading file: {}", source)
                }
            }
            CliError::DataReading { source } => write!(f, "Error reading data: {}", source),
        }
    }
}

impl From<timsquery::DataReadingError> for CliError {
    fn from(e: timsquery::DataReadingError) -> Self {
        CliError::DataReading {
            source: format!("{:?}", e),
        }
    }
}
