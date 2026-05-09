"""Entrapment classification + FDR walk + plot.

The classification half lives here (Task 6). The FDR walk and plot land in
Task 7. The CLI entry-point lands in Task 8.
"""

from __future__ import annotations

import enum
import re
from pathlib import Path

import ahocorasick  # ty: ignore[unresolved-import]
import polars as pl

_MOD_RE = re.compile(r"\[[^\]]*\]|\([^)]*\)|[0-9.]+")
"""Strip bracketed mods (`[U:4]`, `[42]`), parenthesised mods (`(Phospho)`),
and bare numeric mass shifts (`123.45`)."""


def strip_mods(seq: str) -> str:
    """Strip mod annotations to return a bare AA sequence (alpha chars only)."""
    return _MOD_RE.sub("", seq)


def parse_fasta(path: str | Path) -> dict[str, str]:
    """Parse a FASTA file into {accession: concatenated_sequence}.

    Accession is taken as the full header line minus the leading `>`,
    stripped of trailing whitespace. The full header is used so callers can
    later parse it however they want; we don't impose uniprot's `sp|...|`
    grammar here.
    """
    out: dict[str, str] = {}
    current_acc: str | None = None
    parts: list[str] = []
    with Path(path).open("r") as f:
        for raw_line in f:
            line = raw_line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_acc is not None:
                    out[current_acc] = "".join(parts)
                current_acc = line[1:].strip()
                parts = []
            else:
                parts.append(line)
    if current_acc is not None:
        out[current_acc] = "".join(parts)
    return out


class PeptideClass(enum.Enum):
    TARGET = "target"
    ENTRAPMENT = "entrapment"
    SHARED_DROPPED = "shared_dropped"
    UNKNOWN = "unknown"


def _build_hits(patterns: set[str], proteins: dict[str, str]) -> set[str]:
    """Return the subset of `patterns` that occurs as a substring of any value
    in `proteins`."""
    if not patterns:
        return set()
    aut = ahocorasick.Automaton()
    for pat in patterns:
        aut.add_word(pat, pat)
    aut.make_automaton()

    hits: set[str] = set()
    for seq in proteins.values():
        for _, pat in aut.iter(seq):
            hits.add(pat)
            if len(hits) == len(patterns):
                return hits
    return hits


def classify_peptides(
    results: pl.DataFrame,
    target_fasta: str | Path,
    entrapment_fasta: str | Path,
) -> pl.DataFrame:
    """Add `class` and `is_entrapment` columns to a results DataFrame.

    `results` must have a `sequence` column. Sequences are mod-stripped
    before substring matching. Shared peptides (present in both fastas) are
    marked SHARED_DROPPED -- callers exclude them from FDR.
    """
    if "sequence" not in results.columns:
        raise ValueError("results dataframe missing required 'sequence' column")

    target = parse_fasta(target_fasta)
    entrap = parse_fasta(entrapment_fasta)

    stripped = results["sequence"].map_elements(strip_mods, return_dtype=pl.Utf8)
    patterns = set(stripped.to_list())

    hits_t = _build_hits(patterns, target)
    hits_e = _build_hits(patterns, entrap)

    def _classify(s: str) -> str:
        in_t, in_e = s in hits_t, s in hits_e
        if in_t and in_e:
            return PeptideClass.SHARED_DROPPED.value
        if in_t:
            return PeptideClass.TARGET.value
        if in_e:
            return PeptideClass.ENTRAPMENT.value
        return PeptideClass.UNKNOWN.value

    classes = stripped.map_elements(_classify, return_dtype=pl.Utf8)
    is_entrap = classes == PeptideClass.ENTRAPMENT.value

    return results.with_columns(
        classes.alias("class"),
        is_entrap.alias("is_entrapment"),
    )
