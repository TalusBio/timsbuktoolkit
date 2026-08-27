# Carafe contract

What Carafe assumes when calling `timsquery`. Break any of it → silent failure or NPE.
Refs: `util/CallTimsQuery.java`, `ai/AIGear.java`, `dia/{PSMQuery,PSMQueryResult,XICQueryResult}.java`.

> Vendored from the Carafe repo so the assumptions live next to the code that
> has to honour them.
>
> **Line numbers are deliberately omitted**: an earlier version of this file
> cited `AIGear.java` at `~6660–6912`, which is wrong the moment Carafe lands a
> commit above line 6660. Search for the symbol instead. If you re-verify this
> contract against Carafe, record the commit you checked below.
>
> Last verified against Carafe: `db81731e9efd644b942dece5827cb86eed46fad0`
>
> Two test modules pin the mechanically checkable parts against the literal
> JSON in this document. If you edit a payload here, edit it there:
> - `main.rs` (next to this file) — input (targets, tolerances)
> - `timsquery_cli`'s `carafe_output_contract` module — output (result field
>   names, `-a`/`-f` flag values, the `results.json` basename)

## CLI

```
<binary> query-index -a <aggregator> -r <ms_file> -t <params.json> -e <psm_query.json> -f ndjson -o <out_dir>
```

| Flag | Value |
|---|---|
| subcmd | `query-index` |
| `-a` | `spectrum-aggregator` \| `chromatogram-aggregator` (selects output schema) |
| `-r` | raw `.d` input |
| `-t` | tolerance JSON (§tolerances) |
| `-e` | query targets JSON (§targets) |
| `-f` | `ndjson` |
| `-o` | output **directory** (not file) |

Binary path: `bin/timsquery/{windows,macos,linux}/timsquery_cli[.exe]`.

**Fragile couplings** (most likely to break on a timsquery change):
- Aggregator names are literals — rename breaks Carafe, no fallback.
- `-o` is a dir; Carafe reads `<out_dir>/results.json` — basename must be exactly `results.json`.
- Output parsed line-by-line — must be ndjson, one object per line.
- Exit 0 = success; nonzero logged, no retry.

## Tolerances (`-t`)

```json
{
  "ms":       { "ppm":      [itol-itol_shift, itol+itol_shift] },
  "rt":       { "minutes":  [rt_win, rt_win] },
  "mobility": { "percent":  [mobility, mobility] },
  "quad":     { "absolute": [quad, quad] }
}
```
`ms` unit key is dynamic on Carafe's side (`itolu`, currently `ppm`), but timsquery only accepts `ppm`/`da` (and the `Ppm`/`Absolute` spellings) — any other unit fails to deserialize. Each value is `[low, high]`. A negative `low` is valid and means the window sits entirely above the target mass. Defaults: itol 15, mobility 3.0, quad 0.1, rt_win 0.1 (spectra) / `CParameter.rt_win` (xic).

### `ms` window derivation (`itol` / `itol_shift`)

The `ms` window is `[itol - itol_shift, itol + itol_shift]`, built from two independent inputs:

- **`itol`** = `CParameter.itol` — the configured fragment-ion tolerance (default 15). Static per run.
- **`itol_shift`** = per-run **median observed m/z error**, a data-driven calibration offset measured before the query (`AIGear.java`, search `error_shift`):
  - Carafe collects MS1 and MS2 mass errors from already-matched ions, takes the median of each (`ms1_error_shift`, `ms2_error_shift`).
  - If precursor and fragment units match (`CParameter.tolu` == `CParameter.itolu`): `itol_shift = max(ms1_error_shift, ms2_error_shift)`.
  - Else: `itol_shift = ms2_error_shift` (fragment only).
  - No matched ions → `itol_shift = 0`, so the window collapses to `[itol, itol]`.

timsquery interprets the pair as `[light_magnitude, heavy_magnitude]`: `[2, 7]` = the ion may be up to 2 ppm light **or** up to 7 ppm heavy (signed error range `[-2, +7]`). Under that convention Carafe's `[itol - itol_shift, itol + itol_shift]` is a `±itol` window recentered on the calibration offset `itol_shift`:

- signed error window = `[itol_shift - itol, itol_shift + itol]`
- light magnitude = `itol - itol_shift`, heavy magnitude = `itol + itol_shift`

Example: `itol` 15, `itol_shift` +2 → `[13, 17]` = "13 ppm light to 17 ppm heavy" = a ±15 ppm window centered on the +2 ppm systematic offset. This is the intended behavior — the offset shifts the window's center, `itol` sets its half-width.

## Targets (`-e`, `psm_query.json`)

JSON **array** of:
```json
{ "id": 0, "mobility": 0.95, "rt_seconds": 1234.5, "precursor": 650.32,
  "precursor_charge": 2, "precursor_isotopes": [0,1,2],
  "fragments": [175.1, 288.2], "fragment_labels": ["y1","y3^2"] }
```
`id` = row index, echoed back in results. `rt_seconds` = RT_min × 60. `fragments`↔`fragment_labels` positional. Labels: charge 1 → `y3`, else `y3^2`; precursor → `p`/`p^z`.

## Results (`-o`, `<out_dir>/results.json`)

**ndjson** — one object per line, one line per `id`. Parsed `JSON.parseObject(line, …)`, keyed by `id`.
(`.json` name is misleading; it is not a JSON array.)

### spectrum-aggregator → `PSMQueryResult`
```json
{ "id":0, "mobility_ook0":0.95, "rt_seconds":1234.5, "precursor_mz":650.32,
  "precursor_charge":2, "precursor_intensities":[1200,800,300], "precursor_labels":[0,1,2],
  "fragment_mzs":[175.1,288.2], "fragment_intensities":[500,0] }
```
Scalar intensity per ion (no RT axis). `precursor_intensities`/`fragment_intensities` are **1-D**, positionally paired with their m/z arrays.

### chromatogram-aggregator → `XICQueryResult`
```json
{ "id":0, "mobility_ook0":0.95, "rt_seconds":1234.5,
  "precursor_mzs":[650.32,650.82], "precursor_intensities":[[…],[…]],
  "fragment_mzs":[175.1,288.2], "fragment_labels":["y1","b2"],
  "fragment_intensities":[[…],[…]], "retention_time_results_seconds":[1230,1231] }
```
2-D matrices: `[ion][rt_point]`. Every row length == `retention_time_results_seconds.length`.

## Invariants

1. `id` echoed from input; unique + present. Dup → overwrite; missing → NPE downstream.
2. ndjson: one complete object per line, no array wrapper, no pretty-print.
3. Exact field names (fastjson, no remap). Note **singular vs plural** across modes: spectrum `precursor_mz`+`precursor_labels`(int[]); chromatogram `precursor_mzs`+`precursor_intensities`(2-D). Two distinct schemas.
4. m/z ↔ intensity arrays positionally paired (spectrum 1-D, chromatogram row-per-ion).
5. Output basename exactly `results.json` in `-o` dir. Exit 0 on success.

## Known violations

**Invariant 1, every `id` is present.** The chromatogram path drops targets that
produce no data (`timsquery_cli`'s `processing.rs`, `ExpectedNonEmptyData` →
`None`), so a requested `id` can be absent from `results.json`. Per this
document that NPEs downstream rather than reading as an empty result. Still
open.

Fixed: `id` used to be *renumbered* rather than echoed — the input id was
dropped at `push_row` and results carried the arena position instead. Carafe
was unaffected only because it defines `id` as the row index, so the two agreed
by construction; any other caller got rows silently relabelled `0..n-1`. The
caller's id is now carried in `QueryCollection::source_ids` and reported by
`QueryGeom::output_id`, which falls back to the position for formats that carry no
id (the DIA-NN/Skyline/Spectronaut readers).
