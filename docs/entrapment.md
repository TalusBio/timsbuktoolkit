# Peptide entrapment experiment

`bench/entrapment.py` makes a one-fold peptide entrapment library, searches raw
files with `timsseek`, and estimates false discovery proportion (FDP) from the
reported q-values. It uses the peptide-pair TSV format of
[FDRBench](https://github.com/Noble-Lab/FDRBench). The FDRBench JAR is optional;
`analyze --reference-jar` runs it against the same precursor table and checks
the three FDP estimates.

Build `timsseek` first. Run these commands from the repository root:

```sh
cargo build --release --bin timsseek
uv run python bench/entrapment.py generate \
  --fasta /data/parent.fasta --out /data/entrapment
uv run python bench/entrapment.py run \
  --fasta /data/entrapment/entrapment.fasta \
  --library /data/entrapment/entrapment.mzspeclib.txt.gz \
  --raw /data/sample.d --out /data/entrapment/search
uv run --group interactive python bench/entrapment.py analyze \
  --results /data/entrapment/search/sample/results.parquet \
  --pairs /data/entrapment/peptide_pairs.tsv \
  --out /data/entrapment/sample_fdp
```

Repeat `--raw` to search more files. Add `--reuse-library` on later `run`
invocations. `--config` passes a timsseek config to both library prediction and
search. The generated FASTA and pair file stay separate from the raw search
results. Analyze each sample's `results.parquet` separately.
Add `--reference-jar /path/to/fdrbench.jar` to `analyze` to compare its FDP
estimates with FDRBench.

The generator digests the parent FASTA with trypsin, zero missed cleavages,
7–35 residue peptides, and I converted to L. It keeps one copy of each unique
target peptide, shuffles interior residues for one distinct entrapment, and
keeps the first and last residue fixed. Low-complexity peptides for which it
cannot make a distinct entrapment are omitted and counted in `skipped`.
Peptides with an internal K or R are also omitted so the shuffled sequence
does not gain a new tryptic cleavage site.
`generate` prints these settings before processing. It makes unmodified peptide
pairs; `run` applies timsseek's separate library-prediction modification defaults
(fixed `C[UNIMOD:4]`, variable `M[UNIMOD:35]`, at most one variable mod) unless
`--config` overrides them.
`timsseek build-library` predicts both classes and their pseudo-reversed search
decoys. Search uses `--decoy-strategy never` to keep those decoys and
`--max-qvalue 1` to retain the full reported range.

`analyze` removes search decoys (`is_target=false`), groups modified forms by
stripped peptide and charge, and keeps the form with the lowest q-value, then
highest discriminant score. Peptides absent from the pair file are an error.
An external FDRBench pair file must use its `-I2L` option to match the
generator and the adapter's I/L normalization.
For an externally built library with extra peptides, use `--allow-unpaired` to
drop them explicitly. The output contains `fdrbench_input.tsv`,
`entrapment_fdp.tsv`, `entrapment_fdp.png`, `entrapment_fdp.pdf`, and, when
requested, `fdrbench_reference.csv`. The plots show the three FDP estimates
against reported q-value, both below 0.1 and across the full range. Each panel
retains the last cutoff in up to 1,200 q-value bins.

To plot an existing report without rerunning `analyze`:

```sh
uv run --group interactive python bench/entrapment.py plot \
  --report shitshit/hela_entrapment/fdp/entrapment_fdp.tsv
```

At each q-value cutoff, let `T` be target hits and `E` entrapment hits. The
lower-bound estimate is `E/(T+E)` and the combined estimate is `2E/(T+E)`.
The one-fold paired estimate is `(E + N(p>s>t) + 2N(p>t>s))/(T+E)`, where the
two extra counts compare each entrapment to its paired target at the same
charge. All hits with the same q-value receive the estimate at the end of that
q-value group. These are FDP estimates for the expanded search space; compare
them with the q-values reported by `timsseek` to assess calibration. Exact
score ties between a pair can make FDRBench's sort order differ from this
script's deterministic order. Reference comparison allows the resulting
paired-estimate difference, and still compares the other two estimates exactly.
