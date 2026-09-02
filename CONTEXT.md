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

**Retention time**:
Two different quantities share the name, and a target's `rt_seconds` may hold
either.

- A _measured_ retention time: when the analyte left the column, a duration.
- A _normalized_ retention time, or index: where the analyte sits on a reference
  scale anchored by standard peptides. Biognosys iRT spans about -20 to 120;
  msspeculator writes a PROCAL scale anchored `TFAHTESHISK = 0`,
  `SILDYVSLVEK = 100`. It is dimensionless.

So on a normalized scale **zero and negative values are ordinary**, zero being an
anchor rather than an absence, and a conversion between time units is arithmetic
on an index. Only the source file knows which it wrote, and the arena does not
record it: absent, zero and "index that happens to be zero" are one value.

Treat `rt_seconds` as ordered, not as a duration. Calibration fits a monotone
path from it to observed time, which is why a library on either scale searches
correctly and why nothing has forced the distinction yet.
_Avoid_: iRT (says normalized without saying it is an index), RT (says neither)

### How rows are named

**Arena index**:
Where a row sits in memory. Internal to the process and never an identifier: it
addresses storage, and nothing else.
_Avoid_: id, library_id, row id

**Source id**:
What the source file called a precursor -- the JSON target payload's `id`,
mzSpecLib's `<Spectrum=N>` key, DIA-NN's `transition_group_id`. Opaque: carried
through and echoed back, never used to address anything.

Every row has one. A format that names nothing gets ids minted at load, so a
result is never keyed by where its row happened to land.
_Avoid_: id, library id, output id

### Decoys

**Decoy**:
A deliberately wrong analyte scored alongside real ones to estimate the false
discovery rate. Either shipped by the library or generated as a mass shift.

How a shipped decoy was made is the library's business, not ours, and the methods
have nothing in common:

| method | sequence | composition and mass |
| --- | --- | --- |
| pseudo-reverse | permuted | preserved |
| edge-mutate | substituted | changed |
| mass-shift ±CH2 | identical | changed, by construction |

Nothing is stable across all three, which is why **pairing can only be read, never
derived**. A reader that tries to recognise a decoy's target by sequence, mass or
composition is right for one method and wrong for the others. The only sound
source is an explicit statement by the writer.

**Decoy group**:
The set of analytes that compete, so that one result survives per group and
charge. Two ways a group arises:

- _Declared_ by the library, which is the only way a target and a decoy of that
  target end up in one group: `PEPTIDEK/2` and its reversed partner
  `PEDITPEP/2`.
- _Derived_ when the library declares nothing, where a target competes only with
  its own decoy variants.

A group holds at most one target per charge. More than one means competition
discards a real identification to keep another, so a group spanning several
precursor masses -- two modified forms of one peptide, say -- is a mis-declared
group rather than a large one. A library whose pair identifier names a peptide
rather than a precursor produces exactly that, and the loss is silent: the run
reports a lower target count and nothing else.

Sealing an arena refuses such a group rather than loading it. The arena holds the
decoy flag and the charge beside the group labels, so it is the one place that can
check this for every format at once; the alternative is every reader getting it
right on its own, which is the situation that produced the silent loss in the
first place.

Only mzSpecLib declares them. It is the one format with a way to say that one
entry is the counterpart of another (`related spectrum keys`); a DIA-NN,
Spectronaut or Skyline `transition_group_id` names the row it sits on, not a
partner. So a decoy from any other format is a row flagged as a decoy and
nothing more, and its group is derived.

**Variant**:
One member of a single target's decoy expansion -- the target itself, or one of
its mass-shifted decoys. Scoped to one row, not to a group: a declared group
holds a target and its partners, each row with its own variants. A stored row
expands into several scored variants, so "one row" is not "one result".

**Shipped and derived**:
The two independent ways a scored slot can be a decoy, and neither implies the
other. A shipped decoy is a row the library marked as one; it sits in the
variant-0 slot exactly as a target does. A derived decoy is a mass-shifted
variant of a row that is itself a target. So "is this a decoy" is a question
about the row *and* the variant, and a reader that asks only one of them calls
half the decoys targets -- which leaves the false discovery rate estimated from
nothing while reporting every candidate as passing.

The two never coexist in one arena: a library that ships decoys is expanded 1:1,
and mass-shift variants are derived only when nothing was shipped. That is what
makes the variant-0 slot of a decoy row mean "shipped". It holds because the
decoy policy is applied when rows are read, so it is a property of how an arena
was filled rather than of the arena itself.

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
