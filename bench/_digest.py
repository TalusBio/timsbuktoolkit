"""Trypsin digestion + length filter helpers.

Trypsin rule: cleave C-terminal to K or R, NOT when followed by P.
Missed cleavages are concatenations of N+1 contiguous base segments.
"""

from __future__ import annotations

import re
from pathlib import Path

# Cleavage site: after K or R that is NOT followed by P.
# Use a regex split that emits peptides ending in K|R (or the final tail).
_CLEAVE = re.compile(r"(?<=[KR])(?!P)")


def parse_fasta(path: str | Path) -> dict[str, str]:
    """Parse a FASTA into {accession: concatenated_sequence} (header line minus `>`)."""
    out: dict[str, str] = {}
    cur_acc: str | None = None
    parts: list[str] = []
    with Path(path).open("r") as f:
        for raw in f:
            line = raw.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_acc is not None:
                    out[cur_acc] = "".join(parts)
                cur_acc = line[1:].strip()
                parts = []
            else:
                parts.append(line)
    if cur_acc is not None:
        out[cur_acc] = "".join(parts)
    return out


def digest_protein(sequence: str, missed_cleavages: int = 1) -> list[str]:
    """Digest one protein into peptides. Returns base segments plus all
    contiguous merges of up to (missed_cleavages+1) segments."""
    base = [s for s in _CLEAVE.split(sequence) if s]
    out: list[str] = []
    n = len(base)
    for i in range(n):
        for j in range(i + 1, min(i + 2 + missed_cleavages, n + 1)):
            out.append("".join(base[i:j]))
    return out


def digest_proteins(
    proteins: dict[str, str], missed_cleavages: int = 1
) -> set[str]:
    """Digest a {accession: sequence} dict; return the deduplicated peptide set."""
    out: set[str] = set()
    for seq in proteins.values():
        out.update(digest_protein(seq, missed_cleavages=missed_cleavages))
    return out


def length_filter(peptides: set[str], min_len: int = 7, max_len: int = 30) -> set[str]:
    """Keep peptides with length in [min_len, max_len]."""
    return {p for p in peptides if min_len <= len(p) <= max_len}
