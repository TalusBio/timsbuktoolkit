"""Cut an mzSpecLib library down to the peptides that can elute in an RT slice.

Fits observed RT (seconds) against the library's normalized RT (minutes) on a
results.parquet from a full run, then keeps every spectrum whose predicted
observed RT is inside [lo - margin, hi + margin]. Decoys carry their own RT
attribute and are filtered the same way.

usage: uv run python slice_library.py LIB.mzspeclib.txt.gz RESULTS.parquet LO HI OUT.mzspeclib.txt.gz [MARGIN_S]
"""

import gzip
import re
import sys

import numpy as np
import polars as pl

lib_path, results_path, lo, hi, out_path = sys.argv[1:6]
lo, hi = float(lo), float(hi)
margin = float(sys.argv[6]) if len(sys.argv) > 6 else 90.0

RT_RE = re.compile(r"^\[\d+\]MS:1000896\|normalized retention time=([-\d.eE+]+)")
NAME_RE = re.compile(r"^MS:1003061\|library spectrum name=(.+)$")

# Pass 1: name -> library RT (minutes), to fit the mapping.
name_rt = {}
with gzip.open(lib_path, "rt") as f:
    name = None
    for line in f:
        if line.startswith("<Spectrum="):
            name = None
        elif name is None and line.startswith("MS:1003061|"):
            name = NAME_RE.match(line.rstrip("\n")).group(1)
        elif name is not None and line.startswith("[") and "MS:1000896|" in line:
            m = RT_RE.match(line)
            if m:
                name_rt[name] = float(m.group(1))
print(f"library spectra with RT: {len(name_rt)}", file=sys.stderr)

res = (
    pl.read_parquet(results_path)
    .filter(pl.col("is_target") & (pl.col("qvalue") <= 0.01))
    .select("sequence", "obs_rt_seconds")
)
pairs = [(name_rt[s], o) for s, o in zip(res["sequence"], res["obs_rt_seconds"]) if s in name_rt]
x = np.array([p[0] for p in pairs])
y = np.array([p[1] for p in pairs])
a, b = np.polyfit(x, y, 1)
resid = np.abs(y - (a * x + b))
print(
    f"fit on {len(pairs)} IDs: obs_s = {a:.3f} * lib_min + {b:.1f}; "
    f"residual p50 {np.median(resid):.1f} s, p95 {np.percentile(resid, 95):.1f} s",
    file=sys.stderr,
)

# Pass 2: stream, keep header + spectra whose predicted RT is in range.
keep_lo, keep_hi = lo - margin, hi + margin
kept = total = 0
with gzip.open(lib_path, "rt") as f, gzip.open(out_path, "wt", compresslevel=1) as out:
    block = []
    keep = True  # header block
    in_spectrum = False

    def flush():
        global kept
        if block and keep:
            out.write("".join(block))
            if in_spectrum:
                kept += 1

    for line in f:
        if line.startswith("<Spectrum="):
            flush()
            block = [line]
            keep = False
            in_spectrum = True
            total += 1
            continue
        block.append(line)
        if in_spectrum and line.startswith("[") and "MS:1000896|" in line:
            m = RT_RE.match(line)
            if m:
                pred = a * float(m.group(1)) + b
                keep = keep_lo <= pred <= keep_hi
    flush()
print(f"kept {kept} of {total} spectra for rt {lo}-{hi} s (+-{margin} s) -> {out_path}", file=sys.stderr)
