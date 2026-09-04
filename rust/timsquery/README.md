
# Timsquery

A library for querying and aggregating timsTOF mass spectrometry data with flexible tolerances.

## Status

Early development - API is fast-moving and may change without notice.

## Overview

Timsquery provides modular components for querying indexed timsTOF peak data:

- **Aggregators**: Different ways to collect and aggregate peaks
  - `ChromatogramCollector` - retention time profiles
  - `SpectralCollector` - m/z spectra
  - `MzMobilityStatsCollector` - statistical aggregations over m/z and ion mobility
  - `PointIntensityAggregator` - raw peak intensities

- **Tolerances**: Configure m/z, retention time, ion mobility, and quadrupole isolation tolerances via the `Tolerance` struct

- **Query traits**: `QueriableData` and `GenerallyQueriable` traits allow flexible querying of `IndexedTimstofPeaks` (from `timscentroid` crate)

- **Elution groups**: Group related precursor queries for efficient batch processing

## Usage

Peak indexing and centroiding is handled by the `timscentroid` crate. This crate focuses on querying already-indexed data with different aggregation strategies.

For a command-line interface, see the `timsquery_cli` crate.

## Dependencies

- `timscentroid` - peak detection and indexing
- `micromzpaf` - ion annotation utilities
- `timsrust` - raw timsTOF file reading
