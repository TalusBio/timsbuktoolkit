# mzSpecLib test fixtures

Verbatim from [HUPO-PSI/mzSpecLib](https://github.com/HUPO-PSI/mzSpecLib)
`examples/`, Apache-2.0 (same licence as this project).

| file | why |
|---|---|
| `diann.mzSpecLib.txt` | the shape `speclib_build` will emit; all peaks annotated, every mass error exactly `0.0` |
| `spectronaut.mzSpecLib.txt` | carries `-H2O`/`-NH3` losses and a unit-tagged retention time in minutes |

Both are fully representable, so neither exercises the unknown-label or
skip paths — those are covered by unit tests over `resolve_annotation`.

Deliberately NOT vendored: the NIST and SpectraST examples. They are
dominated by internal fragments, immonium ions and unannotated (`?`) peaks,
and by consensus spectra whose observed m/z carries real calibration error.
Worth adding when the resolution policy needs end-to-end coverage; today
they would only make these two harder to read.
