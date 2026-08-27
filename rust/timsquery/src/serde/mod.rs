pub mod chromatogram_output;
mod diann_io;
pub mod diann_speclib_io;
mod elution_group_inputs;
pub mod index_serde;
mod library_file;
mod mzspeclib_io;
mod skyline_io;
mod spectronaut_io;

pub use chromatogram_output::*;
pub use index_serde::*;
// The reader-internal types (`ElutionGroupCollection`, `FileReadingExtras`, the
// per-format `*PrecursorExtras`, `LibrarySniffError`) are deliberately NOT
// re-exported: they had no consumers outside this module, and every format
// already funnels into `LibraryArena`, which is the boundary worth supporting.
pub use library_file::{
    LibraryArena,
    LibraryReadingError,
    read_library_file,
};
