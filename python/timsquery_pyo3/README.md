# timsquery_pyo3

Python bindings for [timsquery](../../rust/timsquery/), a Rust library for querying timsTOF mass spectrometry data.

## Installation

Requires a Rust toolchain and Python >= 3.9.

```bash
cd python/timsquery_pyo3
pip install maturin
maturin develop --release
```

## Quick start

```python
import timsquery_pyo3 as tq

# 1. Load an index from a .d directory (builds + caches on first run)
#    or from a pre-built .d.idx cache
index = tq.PyTimsIndex("path/to/experiment.d")

# 2. Set up tolerances (or just use defaults: 20ppm, 5min RT, 3% mobility, 0.1Da quad)
tolerance = tq.PyTolerance.default()

# Override a single dimension:
tolerance = tq.PyTolerance.default().with_mz(tq.PyMzTolerance.ppm(10.0, 10.0))

# Or build from scratch:
tolerance = tq.PyTolerance(
    mz=tq.PyMzTolerance.ppm(15.0, 15.0),
    rt=tq.PyRtTolerance.minutes(3.0, 3.0),
    mobility=tq.PyMobilityTolerance.pct(5.0, 5.0),
    quad=tq.PyQuadTolerance.absolute(0.1, 0.1),
)

# 3. Define an elution group (one precursor + its fragments)
eg = tq.PyElutionGroup(
    id=1,
    precursor_mz=500.0,
    precursor_charge=2,
    rt_seconds=300.0,
    mobility=0.85,
    fragment_mzs=[600.1, 700.2, 800.3],
    fragment_labels=[0, 1, 2],
    precursor_labels=[0, 1, -1],  # M0, M+1, M-1 isotopes (optional, defaults to [0])
)

# 4. Query
result = index.query_chromatogram(eg, tolerance)

# 5. Access results as numpy arrays
result.fragment_intensities   # shape (n_fragments, n_cycles), dtype float32
result.precursor_intensities  # shape (n_precursors, n_cycles), dtype float32

result.fragment_labels        # [(label, mz), ...] — row order matches the array
result.precursor_labels       # [(isotope_offset, mz), ...]
result.rt_range_ms            # (start_ms, end_ms)
result.num_cycles             # number of RT points
result.id                     # elution group id
```

## Batch queries (parallel)

```python
elution_groups = [eg1, eg2, eg3, ...]  # list of PyElutionGroup

# Shared tolerance (applied to all queries)
results = index.query_chromatograms_batch(elution_groups, tolerance)

# Per-query tolerance (list must match length)
results = index.query_chromatograms_batch(elution_groups, [tol1, tol2, tol3])
```

## Aggregator reuse

Avoid repeated allocations by reusing a `ChromatogramResult` across queries.
The internal `Vec<f32>` capacity grows to the largest elution group and stays there.

```python
result = index.query_chromatogram(eg1, tolerance)   # first query — allocates

index.query_chromatogram_into(result, eg2, tolerance)  # reuses allocation
index.query_chromatogram_into(result, eg3, tolerance)  # same allocation
```

## Streaming queries (iterator in, iterator out)

For large-scale workloads, stream elution groups from any Python iterator.
Internally uses chunked rayon parallelism and reuses collector allocations
across chunks — after the first chunk, allocations settle and only `memcpy`
into numpy remains.

```python
# Any iterable works — generator, list, map, etc.
eg_iter = (make_eg(row) for row in dataframe.itertuples())

# Shared tolerance
for arrays in index.query_chromatograms_iter(eg_iter, tolerance, chunk_size=256):
    arrays.id                      # int
    arrays.precursor_intensities   # np.ndarray (n_prec, n_cycles), float32
    arrays.fragment_intensities    # np.ndarray (n_frag, n_cycles), float32
    arrays.precursor_labels        # list[(isotope_offset, mz)]
    arrays.fragment_labels         # list[(label, mz)]
    arrays.rt_range_ms             # (start_ms, end_ms)
    arrays.num_cycles              # int

# Per-query tolerance (iterable consumed in lockstep with elution groups)
tol_iter = (make_tol(row) for row in dataframe.itertuples())
for arrays in index.query_chromatograms_iter(eg_iter, tol_iter, chunk_size=256):
    ...
```

`ChromatogramArrays` is a lightweight frozen object that owns its numpy arrays.
The iterator's internal collector pool is never exposed — it just keeps reusing
the same Rust-side buffers across chunks.

## Spectral queries

### Summed intensity per ion

```python
result = index.query_spectrum(eg, tolerance)
result.precursor_intensities   # list[float] — one total intensity per precursor
result.fragment_intensities    # list[float] — one total intensity per fragment
result.precursor_labels        # list[(isotope_offset, mz)]
result.fragment_labels         # list[(label, mz)]
result.id                      # int
```

### Intensity-weighted mean m/z and mobility

```python
result = index.query_mz_mobility(eg, tolerance)
result.precursor_stats   # list[(weight, mean_mz, mean_mobility)]
result.fragment_stats    # list[(weight, mean_mz, mean_mobility)]
result.precursor_labels  # list[(isotope_offset, mz)]
result.fragment_labels   # list[(label, mz)]
result.id                # int
```

Each stats tuple contains:
- `weight` — total accumulated intensity
- `mean_mz` — intensity-weighted mean m/z (NaN if no peaks found)
- `mean_mobility` — intensity-weighted mean ion mobility in 1/K0 (NaN if no peaks found)

## Tolerance reference

Each dimension has its own type with `@staticmethod` constructors:

| Type | Constructors |
|---|---|
| `PyMzTolerance` | `.ppm(low, high)`, `.absolute(low, high)` |
| `PyRtTolerance` | `.minutes(low, high)`, `.pct(low, high)`, `.unrestricted()` |
| `PyMobilityTolerance` | `.absolute(low, high)`, `.pct(low, high)`, `.unrestricted()` |
| `PyQuadTolerance` | `.absolute(low, high)` |

Tolerances are symmetric ranges expressed as positive values.
A tolerance of `(5.0, 5.0)` on a value of `100.0` gives the range `[95.0, 105.0]`.

`PyTolerance` supports a builder pattern for overriding individual dimensions:

```python
tol = tq.PyTolerance.default()              # start from defaults
tol = tol.with_mz(tq.PyMzTolerance.ppm(10.0, 10.0))
tol = tol.with_rt(tq.PyRtTolerance.unrestricted())
tol = tol.with_mobility(tq.PyMobilityTolerance.absolute(0.05, 0.05))
tol = tol.with_quad(tq.PyQuadTolerance.absolute(0.2, 0.2))
```

## Lazy vs eager loading

```python
# Eager (default): loads entire index into memory — faster queries
index = tq.PyTimsIndex("experiment.d")

# Lazy: loads from cached .idx on demand — faster startup, lower memory
index = tq.PyTimsIndex("experiment.d.idx", prefer_lazy=True)

index.is_lazy  # bool
```

## RT / cycle mapping

```python
index.num_cycles          # total MS1 cycles in the acquisition
index.rt_range_ms         # (start_ms, end_ms)
index.rt_values_ms        # list[int] — RT in ms for every cycle index

# Convert between seconds and cycle indices
idx = index.rt_seconds_to_cycle_index(300.0)   # nearest cycle index
rt  = index.cycle_index_to_rt_ms(idx)          # back to ms (raises IndexError if OOB)
```

Useful for building an RT axis aligned with chromatogram arrays:

```python
import numpy as np
result = index.query_chromatogram(eg, tolerance)
rt_axis = np.array(index.rt_values_ms, dtype=np.float32) / 1000.0  # seconds
# rt_axis and result.fragment_intensities share the cycle dimension
```

## Current limitations

- **Fragment keys are `usize` only.** The Rust library is generic over key types
  (e.g. `IonAnnot`), but this binding fixes `T = usize` for simplicity.
- **Intensity values are `f32` only.**
- **PointIntensityAggregator not yet exposed.**

## Roadmap

- [x] **Aggregator reuse** — `query_chromatogram_into` reuses a `ChromatogramCollector`
      allocation across queries, avoiding repeated allocation.
- [x] **SpectralCollector** — `query_spectrum` (summed f32) and `query_mz_mobility`
      (intensity-weighted mean m/z + mobility) per ion.
- [ ] **PointIntensityAggregator** — single scalar total intensity per elution group.
- [ ] **IonAnnot key type** — support `IonAnnot` fragment labels alongside `usize`,
      enabling richer annotation round-trips between Python and Rust.
- [ ] **Zero-copy array access** — return numpy views backed by Rust-owned memory
      instead of copying, for large-scale workloads.
- [x] **CycleToRTMapping exposure** — `rt_seconds_to_cycle_index`, `cycle_index_to_rt_ms`,
      `rt_values_ms`, `num_cycles`, `rt_range_ms` on `PyTimsIndex`.
- [ ] **Library file I/O** — read DIA-NN / Spectronaut libraries directly into
      lists of `PyElutionGroup`, removing boilerplate on the Python side.
- [x] **Streaming queries** — `query_chromatograms_iter` streams from any Python
      iterator with chunked rayon parallelism and internal collector reuse.
