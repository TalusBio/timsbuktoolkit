# /// script
# /// dependencies = [
# ///   "timsquery_pyo3",
# ///   "numpy",
# /// ]
# ///
# /// [tool.uv.sources]
# /// timsquery_pyo3 = { path = ".." }
# ///
"""
Streaming chromatogram extraction from timsTOF data.

This example demonstrates the three query modes in timsquery_pyo3:

  1. Single query           — one elution group at a time
  2. Aggregator reuse       — reuse allocations across sequential queries
  3. Streaming iterator     — iterator-in, iterator-out with chunked parallelism

Usage:
    uv run examples/streaming_example.py <path_to_experiment.d>
"""

import sys
import time

import numpy as np
import timsquery_pyo3 as tq


# ---------------------------------------------------------------------------
# Example elution groups (from the timsquery_cli templates)
# ---------------------------------------------------------------------------

EXAMPLE_ELUTION_GROUPS = [
    dict(
        id=0,
        precursor_mz=723.844601280237,
        precursor_charge=2,
        rt_seconds=302.2712,
        mobility=0.9851410984992981,
        fragment_mzs=[147.1128, 74.06004, 248.1604, 124.58387, 347.22889, 174.11808, 418.26601],
        fragment_labels=[0, 1, 2, 3, 4, 5, 6],
        precursor_labels=[0, 1, 2],
    ),
    dict(
        id=1,
        precursor_mz=723.844601280237,
        precursor_charge=1,
        rt_seconds=354.2712,
        mobility=0.9851410984992981,
        fragment_mzs=[147.1128, 74.06004, 248.1604, 124.58387, 347.22889, 174.11808, 418.26601],
        fragment_labels=[0, 1, 2, 3, 4, 5, 6],
        precursor_labels=[0],
    ),
]


def make_elution_groups(n: int):
    """Generate n elution groups by cycling through the templates with shifted RTs."""
    for i in range(n):
        template = EXAMPLE_ELUTION_GROUPS[i % len(EXAMPLE_ELUTION_GROUPS)]
        yield tq.PyElutionGroup(
            id=i,
            precursor_mz=template["precursor_mz"],
            precursor_charge=template["precursor_charge"],
            # Shift RT slightly for each group so they aren't identical
            rt_seconds=template["rt_seconds"] + (i * 0.5),
            mobility=template["mobility"],
            fragment_mzs=template["fragment_mzs"],
            fragment_labels=template["fragment_labels"],
            precursor_labels=template["precursor_labels"],
        )


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("Error: please provide a path to a .d or .d.idx file.")
        sys.exit(1)

    data_path = sys.argv[1]
    n_queries = 200

    # ------------------------------------------------------------------
    # Load the index
    # ------------------------------------------------------------------
    print(f"Loading index from: {data_path}")
    t0 = time.perf_counter()
    index = tq.PyTimsIndex(data_path)
    print(f"  loaded in {time.perf_counter() - t0:.2f}s  ({index})")

    # ------------------------------------------------------------------
    # Set up tolerances — narrow search window
    # ------------------------------------------------------------------
    tolerance = tq.PyTolerance(
        mz=tq.PyMzTolerance.ppm(10.0, 10.0),
        rt=tq.PyRtTolerance.minutes(0.5, 0.5),
        mobility=tq.PyMobilityTolerance.pct(5.0, 5.0),
        quad=tq.PyQuadTolerance.absolute(1.05, 1.05),
    )

    # ------------------------------------------------------------------
    # Mode 1: Single queries
    # ------------------------------------------------------------------
    print(f"\n--- Mode 1: Single queries ({n_queries} queries) ---")
    egs = list(make_elution_groups(n_queries))

    t0 = time.perf_counter()
    for eg in egs:
        result = index.query_chromatogram(eg, tolerance)
    dt = time.perf_counter() - t0
    print(f"  {dt:.3f}s total, {dt / n_queries * 1000:.2f}ms per query")
    print(f"  last result: {result}")

    # ------------------------------------------------------------------
    # Mode 2: Aggregator reuse (query_chromatogram_into)
    # ------------------------------------------------------------------
    print(f"\n--- Mode 2: Aggregator reuse ({n_queries} queries) ---")

    t0 = time.perf_counter()
    result = index.query_chromatogram(egs[0], tolerance)
    for eg in egs[1:]:
        index.query_chromatogram_into(result, eg, tolerance)
    dt = time.perf_counter() - t0
    print(f"  {dt:.3f}s total, {dt / n_queries * 1000:.2f}ms per query")
    print(f"  last result: {result}")

    # ------------------------------------------------------------------
    # Mode 3: Streaming iterator (query_chromatograms_iter)
    # ------------------------------------------------------------------
    print(f"\n--- Mode 3: Streaming iterator ({n_queries} queries, chunk_size=64) ---")

    t0 = time.perf_counter()
    total_signal = 0.0
    count = 0
    for arrays in index.query_chromatograms_iter(
        make_elution_groups(n_queries), tolerance, chunk_size=64
    ):
        total_signal += arrays.fragment_intensities.sum()
        count += 1
    dt = time.perf_counter() - t0
    print(f"  {dt:.3f}s total, {dt / n_queries * 1000:.2f}ms per query")
    print(f"  yielded {count} results, total fragment signal: {total_signal:.1f}")

    # ------------------------------------------------------------------
    # Inspect one result in detail
    # ------------------------------------------------------------------
    print("\n--- Inspecting a single result ---")
    eg = egs[0]
    result = index.query_chromatogram(eg, tolerance)
    print(f"  Elution group: {eg}")
    print(f"  Result: {result}")
    print(f"  Precursor shape: {result.precursor_intensities.shape}")
    print(f"  Fragment shape:  {result.fragment_intensities.shape}")
    print(f"  Precursor labels: {result.precursor_labels}")
    print(f"  Fragment labels:  {result.fragment_labels}")
    print(f"  RT range (ms): {result.rt_range_ms}")
    print(f"  Num cycles: {result.num_cycles}")
    print(f"  Fragment TIC per ion: {result.fragment_intensities.sum(axis=1)}")


if __name__ == "__main__":
    main()
