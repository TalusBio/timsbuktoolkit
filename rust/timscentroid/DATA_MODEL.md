# Data model: raw frame to indexed peak

What a timsTOF `.d` holds, what a frame looks like in memory, what the
centroider does to it, and what the index keeps. Numbers quoted are from a
timsTOF Ultra 2 run (`250225_Desnaux_200ng_Hela_ICC_off_DIA.d`, 21 min,
diaPASEF) and will differ on other instruments and methods.

## The raw file

A `.d` directory holds two files that matter here:

- `analysis.tdf`, SQLite. `Frames` has one row per frame: `Id`, `Time`
  (seconds), `MsMsType` (0 = MS1, 9 = diaPASEF MS2), `NumPeaks`, `NumScans`,
  and the offset of its blob. Other tables carry the per-frame calibration
  coefficients, the quadrupole isolation scheme per window group, and
  `GlobalMetadata` (digitizer, acquisition m/z range, instrument).
- `analysis.tdf_bin`, one binary blob per frame.

`NumPeaks` per frame is known before any blob is read. Summing
`min(NumPeaks, max_peaks)` over frames is an upper bound on what the index
will hold, and it is tight when the cap does not bite.

## A frame in memory

```
Frame {
    meta:  FrameMeta  { index, rt_in_seconds, ms_level, window_group,
                        intensity_correction_factor, calibration }
    peaks: FramePeaks { scan_offsets: Vec<u32>,   // len = n_scans + 1
                        tof_indices:  Vec<u32>,   // len = n_peaks
                        intensities:  Vec<u32> }  // len = n_peaks
}
```

`FramePeaks` is a CSR layout over scans. The peaks of scan `s` are the entries
`scan_offsets[s]..scan_offsets[s+1]`. Within a scan, TOF indices are ascending.
There is no per-peak scan column; `iter_peaks()` expands it from the offsets.

### Each entry is a line, not a sample

The firmware has already peak-picked every scan. One entry is one ion on one
scan: its TOF index is the integer bin the firmware's centroid fell in, and its
intensity is the summed counts. Consecutive entries within a scan sit 5 or more
bins apart in 99.9% of cases, 11 or more in 95%. The instrument's own peak width
estimate is 25 ppm, roughly 4 to 8 bins, so a real peak profile is collapsed to
one index per scan. Two consequences:

- A TOF index carries up to half a bin of quantization error. Sub-bin mass
  precision comes from averaging the same ion across scans and frames.
- Along the TOF axis there is nothing to "merge" except that quantization
  jitter: the same ion on two scans landing in neighboring bins. Measured
  across an ion's whole mobility cloud, its TOF indices span 3 bins at p50 and
  5 at p90, with no correlation to scan index. It is random per-scan jitter,
  not a tilt.

### TOF index to m/z

```
t   = index * digitizer_timebase + delay
m/z = c1 * (t - c0)^2 / 1e12
```

Coefficients are per calibration id; each frame names one and carries `t1`/`t2`
for temperature compensation, so the converter is built per frame
(`Tof2MzConverter2::try_from_calibration`). m/z grows with the square of the
index, so one bin is a different width at every m/z. On the Ultra 2 the
636029 bins span 100 to 1700 m/z:

| m/z  | TOF index | bin width  |
|------|-----------|------------|
| 200  |  84 364   | 1.4 mDa, 6.9 ppm |
| 500  | 251 734   | 2.2 mDa, 4.4 ppm |
| 700  | 335 166   | 2.6 mDa, 3.7 ppm |
| 1000 | 440 357   | 3.1 mDa, 3.1 ppm |
| 1700 | 636 029   | 4.0 mDa, 2.4 ppm |

In Da the bin grows with the square root of m/z; in ppm it shrinks with it. A
tolerance in ppm is therefore a different number of bins at each m/z, and a ppm
value below half a bin rounds to zero bins. A tolerance in bins means the same
thing everywhere on the axis. See `MzTolerance`.

### Scan index to mobility

Linear and decreasing: scan 0 is the highest 1/K0. About 0.00114 1/K0 per scan
on a 709-scan ramp. Unlike TOF, the scan axis is a true profile: an ion drifts
out of the TIMS tunnel over consecutive scans and appears on 5 to 20 of them,
with intensity rising and falling. A 3% tolerance at 1/K0 = 1.0 is about 26
scans either side. See `ImTolerance`.

### MS level and window groups

A diaPASEF cycle is one MS1 frame followed by N MS2 frames. Each MS2 frame
carries a window group: a list of (scan range, isolation m/z range) pairs, so
different scans of the same frame hold fragments from different precursor
windows. The same window group recurs every cycle. MS1 frames have none.

`intensity_correction_factor` is per frame and is applied on read
(`CorrectedFramePeak::corrected_intensity`, an `f64`).

## Centroiding, per frame

`PeakCentroider::centroid_frame` takes one frame and returns clusters. Sorted by
TOF, walked in descending order of neighborhood intensity; each unclaimed peak
becomes a parent and absorbs unclaimed peaks within `mz_tol` along TOF and
`im_tol` along scans. Absorption is only ever "uphill": a peak joins a parent at
least as intense as itself. The centroid is the intensity-weighted mean TOF and
scan of the cluster, its intensity the sum.

Three things bound a frame's output:

- `max_peaks`: whole-frame cap, in descending-intensity order, so the most
  intense clusters anywhere win.
- `window_cap`: per m/z window cap, the most intense clusters in each
  `window_da` slice win. Bounds output at `windows * max_peaks` per frame
  regardless of sample.
- `early_stop_iterations`: stop after that many consecutive singleton parents,
  on the assumption the rest is noise. Zero disables.

Peaks dropped by the window cap are never absorbed by a later, weaker parent;
that would pull the weaker centroid toward a peak that outranks it.

Growth is transitive by default: when the walk reaches a peak that was already
absorbed, that peak recruits its own unclaimed neighbors into its parent's
cluster, DBSCAN-style. With a zero-bin m/z tolerance this can only chain along
scans at one TOF index, which is how a mobility cloud wider than `im_tol` still
ends up in one cluster. With a one-bin tolerance it chains along TOF through
crowded regions and pulls centroids off their ions; `transitive = false` bounds
a cluster to the parent's own box instead.

What the defaults do, then: merge one ion's scans at the same TOF index into
one centroid per index it landed in. Since an ion lands in about 3 indices
across its cloud, it becomes 2 to 3 centroids per frame. Measured on 400
confident precursors: 2.1 centroids per MS1 isotope on the apex frame, 2.7 per
fragment isotope. Merging those with `Bins(1)` or `Bins(2)` lowered both
identifications and isotope-envelope fidelity on the runs tried, so the split
is left in place until the reason is understood.

`ClusteringSummary` reports initial, aggregated, final and window-capped counts
and the stopping reason per frame; `AggregatedClusteringSummary` sums them.

## The indexed peak

```
IndexedPeak<T> { mz: f32, intensity: f32, mobility_ook0: f16, cycle_index: T }
```

16 bytes with padding. What is gone relative to the raw entry: the scan index
(now an `f16` 1/K0, three significant digits), the TOF index (now `f32` m/z,
about 0.06 ppm at m/z 1000, well under a bin), and the frame index (now a cycle
index, see below). A build holds roughly 22 bytes per stored peak at its
transient peak, since frames are collected into per-frame vectors before the
columns are assembled.

### Cycle index replaces frame index

`MS1CycleIndex(k)` is the k-th MS1 frame. `WindowCycleIndex(k)` is the k-th
occurrence of a given window group. Each peak group carries a
`CycleToRTMapping` from cycle to retention time in milliseconds, so RT queries
resolve to a cycle range first. The type parameter keeps an MS1 cycle from
being used to query an MS2 group.

## The index

```
IndexedTimstofPeaks {
    ms1_peaks:          IndexedPeakGroup<MS1CycleIndex>,
    ms2_window_groups:  Vec<(QuadrupoleIsolationScheme, IndexedPeakGroup<WindowCycleIndex>)>,
    mobility_kind:      MobilityKind,
}
IndexedPeakGroup<T> {
    cols:              PeakColumns<T>,        // SoA: mz, intensity, mobility, cycle
    bucket_mz_ranges:  Vec<TupleRange<f32>>,
    bucket_size:       usize,
    cycle_to_rt_ms:    CycleToRTMapping<T>,
}
```

One group for MS1, one per distinct quadrupole setting for MS2. A group is
sorted by m/z and cut into buckets of `bucket_size` peaks; bucket `k` occupies
column range `[k*bucket_size, (k+1)*bucket_size)` and `bucket_mz_ranges[k]` is
its m/z span. Inside a bucket, peaks are sorted by cycle. A query binary-searches
the bucket ranges for the m/z window, binary-searches each bucket for the cycle
range, then filters mobility linearly. `rebucket` rebuilds at a different
`bucket_size`; the search uses 256.

An MS2 query with a precursor m/z first selects the window groups whose
isolation scheme intersects it (and the mobility range, if given), then runs
the same group query on each.

The on-disk `.idx` is the same structure as Parquet columns plus a JSON manifest.
The centroiding config is not recorded in it.

## Sizes on the Ultra 2 run

Without an MS2 window cap (MS1 20k/frame, MS2 50k/frame, `Bins(0)` m/z, 2% /
3% mobility), which was the default before the cap:

| | MS1 | MS2 |
|---|---|---|
| frames | 1 694 | 13 548 |
| raw entries | 492 M | 1 050 M |
| stored centroids | 32 M | 559 M |
| frames hitting `max_peaks` | 91% | 63% |

At 16 bytes a peak the resident index is 9.4 GB and the build peaks at 12.4 GB.
Stored MS2 is 17 times MS1, so MS2 is the memory lever; MS1 is the raw-read
lever. With `Bins(0)`, MS2 centroiding merges nothing along TOF (see "Each
entry is a line"), so stored MS2 is raw MS2 merged only along scans.

With the default MS2 `window_cap` of 500 per 50 Da the same run stores about
100 M MS2 centroids, the build peaks at 4.5 GB, and identifications at 1% FDR
go up by about 10%, because what the cap removes is dim clutter the scorer was
not using. Capping MS1 the same way costs identifications; its 20k global cap
is left alone.
