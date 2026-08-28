# Timsbuktoolkit

Querying and searching timsTOF mass spectrometry data. Two crates dominate the
vocabulary and they mean different things by similar words: `timsquery` extracts
signal for a set of coordinates, `timsseek` scores peptide-data matches.

## Language

### What gets searched for

**Target**:
One analyte's constraints in (m/z, retention time, ion mobility) space, plus the
fragment m/z values to extract.
_Avoid_: elution group, precursor, row

**Query**:
A target together with the tolerances applied to it against an index. Targets are
data; a query is the act of looking.

**Library**:
A collection of targets that also carries expected or observed fragment
intensities. Intensities are what make it a library -- without them it is targets.
Note this is narrower than the field's usage, where "spectral library" covers
anything mapping analytes to m/z values.
_Avoid_: spectral library (ambiguous), speclib

**Fragment label**:
The mzPAF annotation identifying a fragment ion (`y3`, `b2^2`, `y5-H2O`). Carries
ion chemistry; distinct from an opaque string label, which does not.

### How rows are named

**Arena index**:
A row's position in the in-memory columnar store. Self-incrementing, assigned on
insertion, meaningful only within one process. Feeds decoy grouping and q-value
determinism, so it is never caller-supplied.
_Avoid_: id, library_id, row id

**Source id**:
What the source file called a precursor -- the JSON target payload's `id`,
mzSpecLib's `<Spectrum=N>` key, DIA-NN's `transition_group_id`. Opaque: carried
through and echoed back, never used to address anything. Absent in some formats.
_Avoid_: id, library id

**Output id**:
The identifier a result carries, so a caller can map results onto the request
they sent. The source id where there is one, the arena index otherwise.

### Decoys

**Decoy**:
A deliberately wrong analyte scored alongside real ones to estimate the false
discovery rate. Either shipped by the library or generated as a mass shift.

**Decoy group**:
A target and its decoy variants, competing as a unit so exactly one survives.

**Variant**:
One member of a decoy group -- the target itself, or one of its mass-shifted
decoys. A stored row expands into several scored variants, so "one row" is not
"one result".

### Capabilities

**Capability**:
Something a loaded library supports, decided at load time from the file's
contents: whether sequences are available, whether fragment labels carry ion
chemistry, how isotopes are derived, how decoys are obtained.

**Graceful degradation**:
Proceeding with a capability absent rather than failing -- skipping FDR when there
are no decoys, skipping sequence-dependent scores when sequences are unavailable.
The gate is per-score, not per-run; only some scores need sequences.

### External contracts

**Contract**:
Assumptions another tool makes about this one's command-line surface or output
format, vendored here because breaking them fails silently on the other side.
See `rust/timsquery/tests/carafe_contract/`.
