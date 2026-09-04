# Centroiding: per-window cap, tolerance units, transitivity

Findings from one day of measurement on `250225_Desnaux_200ng_Hela_ICC_off_DIA.d`
(timsTOF Ultra 2, 21 min diaPASEF, 3.8 GB) with `hela_gt20peps.fasta`, on a
16 GB / 8 core laptop. Search = `timsseek --fasta ... --raw-inputs ...`, IDs =
targets at q <= 0.01 in `results.parquet`. Run-to-run noise on IDs at a fixed
config is about +-50 (4433 vs 4481 for the same setting on two runs).

## 1. Memory is MS2, raw reads are MS1

| | MS1 | MS2 |
|---|---|---|
| frames | 1 694 | 13 548 |
| raw entries | 492 M | 1 050 M |
| stored centroids, default config | 32 M | 559 M |
| frames hitting the global `max_peaks` | 91% | 63% |

Resident index 9.4 GB, build peak 12.4 GB, full search peak 13.4 GB and 2.3 GB
of swap growth. A 60-minute gradient of the same sample would not fit.

`sum(min(NumPeaks, max_peaks))` from the SQLite `Frames` table bounds stored
peaks before any blob is read: 613 M bound vs 590 M actual here. A pre-flight
refuse-or-warn against free memory is one query. Note that RSS is useless for
this on macOS: it stayed at 5.3 GB while the footprint hit 12.4 GB; use
`footprint -p` / `phys_footprint`, or cgroup `memory.current` on Linux.

## 2. Per-window cap on MS2: 3x less memory, ~10% more IDs

`window_cap = { max_peaks = N, window_da = W }` keeps the N most intense
centroids per W Da per frame, inside the clustering loop (dropped peaks are
never absorbed by a weaker later parent). MS1 left at its 20k global cap.

| MS2 cap | peaks / Da | peak GB | q<=0.01 | q<=0.05 | `no_observed_signal` |
|---|---|---|---|---|---|
| none (baseline) | | 13.4 | 4036 | 4847 | 29 |
| 750 / 100 | 7.5 | 4.9 | 4268 | 5137 | 1143 |
| 1000 / 100 | 10 | 4.4 | 4372 | 5173 | 522 |
| 500 / 50 | 10 | 4.5 | 4433, 4481 | 5230 | 467 |
| 250 / 25 | 10 | 5.0 | 4437 | 5226 | 451 |
| 1500 / 100 | 15 | 5.9 | 4328 | 5290 | 201 |
| 2000 / 100 | 20 | 7.1 | 4412 | 5340 | 110 |
| 1000 / 50 | 20 | 7.1 | 4480 | 5340 | 118 |
| 4000 / 200 | 20 | 7.3 | 4495 | 5252 | 129 |
| 3000 / 100 | 30 | killed at 8.3 | | | |

IDs plateau at ~4400 from 7.5 to 20 peaks/Da; memory doubles from 10 to 20.
Finer windows edge out coarser at fixed density (within noise). Capping MS1 as
well costs IDs: MS1 3000/100 + MS2 2000/100 = 4317, MS1 2000/100 + MS2
1500/100 = 4161. Global MS2 `max_peaks` = 20k gives 4128 at 7.7 GB and 10k
gives 3880: the window cap wins on both axes.

**Recommendation**: default MS2 `window_cap = { max_peaks = 500, window_da = 50 }`.
`no_observed_signal` tracks the cap's damage well and belongs in the run report.

## 3. Raw TOF axis is line spectra; one bin is 2.4 to 7 ppm

`DigitizerNumSamples` = 636029 spans 100..1700 m/z. m/z ∝ index², so one bin
is 6.9 ppm at m/z 200 and 2.4 ppm at 1700 (1.4 to 4.0 mDa). Within a scan,
consecutive raw entries sit >= 5 bins apart in 99.9% of cases: the firmware has
already picked one integer index per ion per scan. Across an ion's mobility
cloud its indices span 3 bins (p50) to 5 (p90), uncorrelated with scan: random
per-scan jitter, not a tilt.

The old tolerances, 0.5 ppm MS1 / 1 ppm MS2, convert to a TOF delta and round
to **zero bins**. They never merged across bins. `MzTolerance::Bins(0)` is now
the default and is byte-identical on every metric. Mobility is a real profile
(5 to 20 scans per ion) and the 2% / 3% tolerances cover it.

## 4. Merging across bins is worse, and transitivity is why it collapses

| config (MS2 cap 500/50 throughout) | peak GB | q<=0.01 |
|---|---|---|
| Bins(0), transitive (default) | 5.1 | 4481 |
| MS2 Bins(1) | 4.6 | 4101 |
| MS2 Bins(2) | 4.7 | 4201 |
| MS1 Bins(1) | 4.8 | 4232 |
| both Bins(1) | 4.5 | 3979 |
| both Bins(0), **transitive = false** | 4.8 | **4544** |
| both Bins(1), transitive = false | 4.6 | 4269 |
| both Bins(1), transitive = false, im 4% / 5% | 4.5 | 4455 |
| both Bins(0), transitive = false, im 4% / 5% | 5.0 | 4289 |

Growth is transitive: an absorbed peak recruits its own neighbors into its
parent's cluster. At zero bins that chains only along scans at one index. At
one bin it chains along TOF through crowded regions: with transitive Bins(1),
39% of confident precursors were still visible at their m/z on the apex frame
(99% at Bins(0)); with `transitive = false` that returns to 99%. So the ID
loss from merging is partly chaining. But even bounded, Bins(1)/Bins(2) lose
IDs and lower isotope-envelope fidelity (below). Why a weighted-mean centroid
of one ion's jitter scores worse than the split pieces is not understood;
a plausible reason is that the scorer's mobility-profile features see a
flattened profile after merging. Left as an open question; defaults unchanged
except the Bins(0) spelling. `transitive = false` at Bins(0) is +1.4% IDs and
slightly better envelopes, within noise; worth a second sample before flipping.

## 5. Isotope envelopes as a centroiding metric (no scorer in the loop)

`examples/isotope_frames.rs`: for 400 confident precursors, centroid only the
MS1 frames around the apex and the MS2 frames of the isolating window group,
read M+0..M+4 (precursor) and M+0..M+2 (every b/y fragment) at +-10 / 15 ppm,
+-2% mobility, and compare to averagine by cosine. Null boxes at +half an
isotope. 15 s per config.

| config | MS1 centroids per isotope at apex | MS1 cosine p50 | MS2 fragments at M+0 | MS2 cosine p50 |
|---|---|---|---|---|
| default | 2.1 (63% split) | 0.951 | 10997 / 17566 | 0.973 |
| MS2 cap 500/50 | 2.1 | 0.951 | 9064 | 0.976 |
| Bins(1) transitive | 1.0 but only 39% visible | 0.976 (n=156) | 7735 | 0.958 |
| Bins(1) non-transitive | 2.1 | 0.928 | 10926 | 0.950 |
| Bins(0) non-transitive | 3.1 | 0.961 | 10988 | 0.973 |
| Bins(2) non-transitive, MS1 im 3% | 1.2 but 81% visible | 0.744 | | 0.953 |

The cap removes dim fragments (null rate 10% to 3%) and raises fragment
cosine. Every merge variant lowers it. The remaining splits at Bins(0) are 2 to
3 centroids ~10 ppm and ~1.2% mobility apart on the apex frame: the same ion's
jitter, one centroid per index it landed in.

## 6. Single-window fixture, and what it found (2026-09-04)

`isotope_frames --rt-center C --rt-window W` keeps only seeds whose apex is
within W seconds of C and centroids just their frames: one cycle is 5 frames
in 33 ms, a 30 s window is ~280 frames and ~150 seeds in 0.7 s. Seed pool is
every 1% ID from one full search (`shitshit/seeds_all.csv`). 26 configs on two
windows (450 s, 660 s) ran in 40 s (`shitshit/probe_sweep.sh`).

Mobility tolerance is the knob the fixture moves, and tighter is better on both
levels, monotonically, even though tighter means more centroids per isotope:

| MS1 `im_tol` | M+0 visible | centroids / isotope | MS1 cosine p50 | >= 0.95 |
|---|---|---|---|---|
| 4% | 86 / 88% | 2.3 | 0.973 / 0.960 | 60 / 53% |
| 2% (default) | 88 / 91% | 2.3 | 0.975 / 0.967 | 59 / 57% |
| 1% | 94 / 93% | 2.6 | 0.972 / 0.979 | 56 / 69% |
| 0.5% | 94 / 94% | 3.5 | 0.982 / 0.979 | 71 / 71% |

MS2 1% instead of 3%: fragment cosine 0.963 / 0.966 vs 0.955 / 0.962. Best
row: non-transitive, MS1 1%, MS2 1%: MS1 0.979 / 0.981 (64 / 74% >= 0.95),
MS2 0.962 / 0.968 (57 / 61%). MS1 `max_peaks` 60k adds noise (null 11-13%),
10k loses precursors (84% visible); early-stop off changes nothing.

Full searches on the fixture's top rows (MS2 cap 500/50 throughout):

| config | peak GB | q<=0.01 | q<=0.05 |
|---|---|---|---|
| reference (2% / 3%, transitive) | 5.1 | 4433, 4481 | 5230 |
| MS1 1%, MS2 1% | 5.1 | 4419 | 5087 |
| MS1 0.5% | 4.7 | 4514 | 5242 |
| **non-transitive, MS1 1%, MS2 1%** | 4.8 | **4624** | 5334 |

The fixture's best row is the ID best row, +3 to 4% over the reference and
+15% over the uncapped baseline. Its second-best is second. The middle row
did not separate from the reference on IDs (noise band is +-50). So the
30-second envelope fixture predicts full-run ranking well enough to drive
iteration, at 1/300 of the cost.

### Correction: the readout box was the metric

The cosine differences in the table above are an artifact of the +-2 % / 10 ppm
readout box. Re-read with a 4 % / 20 ppm box, every config, including
transitive Bins(1), has the same envelope fidelity (MS1 cosine 0.996-0.998,
83-89 % >= 0.95; MS2 0.965-0.974), and so does raw, un-centroided data
(0.997 / 0.952-0.958). Intensity is conserved by all of them. The envelope
metric therefore cannot rank centroiding configs; withdraw the "tighter
mobility gives better envelopes" claim and the earlier bins/transitivity
cosine claims.

What the tight box was measuring is **position**: whether the centroid lands
where the raw ion is, within a search-like tolerance. That is the metric that
tracks IDs:

| config | in 10 ppm / 2 % box | centroids / isotope | full-run q<=0.01 |
|---|---|---|---|
| raw (no centroiding) | 99 % | 33-35 | not run; old sweep: best |
| ref (Bins 0, 2 % / 3 %, transitive) | 88-91 % | 4.8-5.3 | 4433-4481 |
| non-transitive, 1 % / 1 % | 94 % | 7.5-8.2 | 4624 |
| MS1 0.5 % | 94 % | | 4514 |
| non-transitive Bins(1) | 94 % | 3.7-4.2 | 4269 |
| transitive Bins(1) | 56-60 % | 1.5 | 4101 |

Mechanism: an intensity-weighted mean across a +-2-3 % mobility box (17-26
scans against a 17-scan cloud) or across chained TOF bins lands off the ion,
by > 2 % of mobility or 10-20 ppm, and the search's +-15 ppm / calibrated
mobility window then misses part or all of it. "Less merge, more IDs" is a
position effect. Raw data is the position ceiling and the memory floor. The
MS2 window cap is the one change that beats raw on any metric: fragment
cosine 0.966-0.972 vs 0.952-0.958 and null 4.5 % vs 20 %, i.e. it removes
clutter without moving anything.

Candidate new defaults stand on the ID and position evidence: `transitive =
false`, MS1 `im_tol` 1 %, MS2 `im_tol` 1 %, MS2 `window_cap` 500/50. Needs a
second sample before flipping. The probe now takes `--box-im-pct`,
`--box-ms1-ppm`, `--box-ms2-ppm`; run a ranking at two widths before believing
it, and read the tight box as "position", the wide box as "intensity".

## 7. Mobility ablation, and what it exposed (2026-09-04)

`[indexing] ignore_mobility = true` sets the index's `MobilityKind` to `Absent`,
which unrestricts mobility in every query (Phase 1, calibration, both Phase 3
passes) and neutralizes the mobility score features. Same centroiding configs,
MS2 cap 500/50:

| config | with mobility | mobility ignored |
|---|---|---|
| ref | 4481 | 10443 |
| non-transitive, 1 % / 1 % | 4624 | 11067 |
| transitive Bins(1) | 4101 | 9314 |

Ranking preserved, so the centroiding effect is not the search's mobility
window fitting itself to the centroids. The absolute counts are not real.
Envelope check on the mobility-off IDs, 15 % box around library 1/K0 (which
covers the library's own error), averagine cosine >= 0.95 as the criterion:

| set | M+0 visible | cosine >= 0.95 |
|---|---|---|
| IDs in both searches (2994) | 96.7 % | 82 % |
| IDs only with mobility (1487) | 96.2 % | 74 % |
| IDs only without mobility (5534 with a library 1/K0) | 39 % | 23 % |
| same set, RT shifted +60 s (control) | 30 % | 2.8 % |

Real peptides score 74-82 % on this criterion, random positions 3 %; the new
set at 23 % is about one quarter real. The mobility-off search reports 1 % FDR
and is nearer 50 %, and the 1487 it dropped are real. The decoy model does not
describe the null once mobility features are gone. Two consequences: the FDR
estimate is fragile to configuration in a way only an entrapment-style check
catches, and this envelope test is a cheap stand-in for one. For any ID set,
`(rate - 0.03) / (0.80 - 0.03)` estimates its real fraction on this run.

### Calibration findings that stand on their own

- **Library 1/K0 is 8.6 % off at the median, 13.5 % at p90** against observed
  for confident IDs. The calibrant query uses a +-5 % box around library
  mobility (`processing.rs`, `MobilityTolerance::Pct((5.0, 5.0))`), so it
  keeps the minority of peptides whose prediction happened to land. Only 42 %
  of confident IDs are visible in a +-5 % box around their library 1/K0.
- **Phase 3 runs at +-60 ppm and +-14 %**, identical across every centroiding
  config (`calibration.json`: `residuals.mz_ppm` 58-63, `mobility_pct` 14).
  The windows are `median +- 3 * 1.4826 * MAD` of one intensity-weighted
  offset per calibrant (`weighted_ms1`), i.e. the mean of every peak inside the
  Phase 1 box at the apex. In a 300k-peak MS1 frame that mean is pulled by
  whatever else sits in +-15 ppm, so the residual MAD is 13.6 ppm and 3.2 %:
  uniform across the box, not the instrument's 2-3 ppm. Phase 3 has been
  running 4x looser than Phase 1. A per-calibrant robust offset (apex peak, or
  median over the box) would fix the derivation; see the sigma sweep below for
  what the loose windows cost.

## 8. Two converters for scan to 1/K0, and a corrected position table

The run-level `metadata.im_converter` and the per-frame calibration the index
build uses disagree by a mobility-dependent amount: 2.7 % at scan 50, 1.9 % at
scan 400, 0.1 % at scan 700. The probe used the run-level one until now, so
every mobility position in sections 6 and 7 carried that bias (it showed up as
a −2.1 % median offset of every fragment and isotope from the seed's observed
mobility, identical for raw data). With the per-frame calibration the median
offset is −0.1 % and the position table becomes:

| config (MS2 cap 500/50) | in 10 ppm / 2 % box | centroids / isotope |
|---|---|---|
| raw, no centroiding | 99.3 % | 24 |
| ref | 98.6 % | 3.9 |
| non-transitive, 1 % / 1 % | 97.9 % | 6.1 |
| 4 % / 5 % | 97.9 % | 3.9 |
| non-transitive Bins(1) | 98.6 % | 2.4 |
| transitive Bins(1) | 93.0 % | 1.1 |

So position no longer separates the mobility-tolerance settings; their ID
differences (4419-4624) are inside noise and the `transitive` / `im_tol`
recommendation is withdrawn to "no evidence either way". Transitive TOF
merging still loses position and IDs. Anything reading mobility off raw
frames must use the per-frame calibration; `PeakCentroider` uses the run-level
converter only for the relative width of its tolerance, where the scale
difference is harmless.

### Within-peptide mobility spread

Brightest centroid of each detected fragment M+0 and each visible precursor
isotope, relative to the seed's observed mobility, over 628 seeds:

| | median | MAD | 3 × 1.4826 × MAD | within 1 % | within 2 % |
|---|---|---|---|---|---|
| fragments M+0 (n=15189) | −0.12 % | 0.80 % | 3.6 % | 55 % | 66 % |
| precursor M+1..M+4 (n=2365) | −0.08 % | 0.31 % | 1.4 % | 80 % | 89 % |

Identical for ref, non-transitive 1 %, and raw data: a property of the run.
This is the quantity a mobility window centered on the observed mobility has
to cover, and it is not what the current derivation measures (that one
measures library-vs-observed, which is 8.6 % at the median). The fragment
tails (a third beyond 2 %) are interferences winning "brightest centroid",
which is also what pulls the pass-2 center when it is an intensity-weighted
mean.

## 9. Phase 3 windows: the largest effect found

`[calibration]` sigma multipliers, default centroiding + MS2 cap 500/50:

| `mz_sigma`, `mobility_sigma` | Phase 3 window | Phase 3 time | q<=0.01 | no_signal |
|---|---|---|---|---|
| 3, 3 (default) | 60 ppm, 14 % | 23 s | 4448 | 467 |
| 0.5, 3 | 10 ppm, 14 % | 9 s | **5967** | 4736 |
| 1, 1 | 20 ppm, 4.7 % | 9 s | 1809 | 5930 |
| 3, 0.5 | 60 ppm, 2.4 % | 16 s | 762 | 5227 |
| 0.5, 0.5 | 10 ppm, 2.4 % | 6 s | 1185 | 34872 |

The calibrant query itself uses a 50 ppm / 5 % box (`processing.rs`,
`query_tolerance`), which is where the contamination enters.

Validity of the +34 %: envelope cosine >= 0.95 at 4 % / 10 ppm around observed
coordinates is 92.8 % for the 4307 shared IDs, 75.4 % for the 1660 new ones,
11.5 % for the new set shifted +60 s. About 80 % of the gain is real
peptides. Tightening mobility to 2.4 % collapses because pass 2 centers on the
intensity-weighted fragment mobility, which is contaminated; a narrow window
around a bad center misses.

Reading: the derived m/z window is wider than the window the calibrants were
found in, which is not possible for a sound derivation. Shipped: clamp it to
the first-pass window and warn (`clamp_mz_window_to_first_pass`). Still to do:
fix the estimator (per-calibrant robust offset from the apex peak, not the
box mean). For mobility, derive the window from the within-peptide spread
above (3σ ≈ 3.6 % for fragments) and center pass 2 on a robust estimate
(median fragment mobility, or the brightest fragment) rather than the weighted
mean; every attempt to tighten it around the current center collapsed.

**Combined new defaults** (MS2 window cap 500/50, m/z window clamped to the
15 ppm first pass, everything else unchanged), one full run:

| | original | new defaults |
|---|---|---|
| q<=0.01 | 4036 | **5817** |
| q<=0.05 | 4847 | 6837 |
| peak footprint | 13.4 GB | 5.2 GB |
| wall | 259 s | 97 s |

## Recommendations, and what shipped

Shipped as defaults (one sample, both effects large and verified by the
envelope check):

- MS2 `window_cap = { max_peaks = 500, window_da = 50 }`. MS1 uncapped by
  window.
- Derived Phase 3 m/z window clamped to the first-pass window, with a warning
  naming the contaminated residual MAD.
- `MzTolerance::Bins(0)` as the spelling of the existing behavior.

Available, default off, for the next round: `[indexing]` per-level overrides
(`mz_tol` in bins or ppm, `im_tol` in pct or scans, `window_cap`, `max_peaks`,
`early_stop_iterations`, `transitive`), `rt_range_seconds` slice builds,
`ignore_mobility` ablation.

Withdrawn: the `transitive = false` / `im_tol` 1 % recommendation. Its ID
edge (4624 vs 4481) is inside noise and the position evidence for it was the
converter bias of section 8.

Open, in order of expected payoff:

1. Calibrant offset estimator: robust per-calibrant offset instead of the
   box mean; then the clamp stops firing.
2. Mobility window from within-peptide spread, pass-2 center from a robust
   fragment mobility. The library-vs-observed 8.6 % error stays a pass-1
   concern.
3. Library 1/K0 prediction is 8.6 % off at the median; either the model or
   the unit needs looking at (see #115 for the RT analogue).
4. The FDR estimate went to ~50 % true error under the mobility ablation while
   reporting 1 %. Entrapment CI (roadmap phase 6) is the real check; the
   envelope real-fraction estimate is the cheap proxy until then.
5. Pre-flight memory gate from `Frames.NumPeaks`; centroiding config into the
   `.idx` manifest; clustering summary into the run report.

## Tooling added (all in `rust/timscentroid/examples/`)

- `peak_density`: centroids per 100 Da window per frame, retention under caps.
- `index_build_rss`: full build under `/usr/bin/time -l`.
- `peak_width`: raw m/z and mobility width per ion, TOF drift across the cloud.
- `tof_bins`: bin width table, within-scan spacing histogram.
- `isotope_frames`: envelope metric above, frame-level, seconds per config.
- `isotope_envelopes`: same on a built index (slower; superseded).

Sweep drivers in `shitshit/` (not for commit): `budget_sweep.sh` runs a
prioritized point list under a wall budget with a `footprint -p` watchdog.

## Config surface added

`[indexing.ms1]` / `[indexing.ms2]` in the search TOML, every field optional:
`max_peaks`, `mz_tol = { Bins = n } | { Ppm = x }`, `im_tol = { Pct = x } |
{ Scans = n }`, `early_stop_iterations`, `window_cap = { max_peaks, window_da }`,
`transitive`. Only read when the index is built from raw.

## Open

- Why merging one ion's jitter into one centroid scores worse. Check what the
  apex / mobility-profile features see before and after.
- Cross-frame persistence filter (drop centroids with no neighbor at cycle
  t-1 / t+1): a different prior from both caps, and the one that removes
  noise rather than redistributing a budget. Needs a persistence-fraction
  measurement first.
- Pre-flight memory gate from `Frames.NumPeaks`.
- Record the centroiding config in the `.idx` manifest and the aggregated
  clustering summary in the run report.
