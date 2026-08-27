# mzSpecLib test fixtures

Verbatim from [HUPO-PSI/mzSpecLib](https://github.com/HUPO-PSI/mzSpecLib)
`examples/`, Apache-2.0 (same licence as this project).

| file | why |
|---|---|
| `diann.mzSpecLib.txt` | the shape `speclib_build` will emit; all peaks annotated, every mass error exactly `0.0` |
| `spectronaut.mzSpecLib.txt` | carries neutral losses (`-H2O`, `-NH3`), which exercise the unrepresentable-annotation path |

Deliberately NOT vendored: the NIST and SpectraST examples. They are
dominated by internal fragments, immonium ions and unannotated (`?`) peaks,
and by consensus spectra whose observed m/z carries real calibration error.
Useful later for the resolution-policy tests, but they would make these two
fixtures harder to read for no gain.
