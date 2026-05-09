"""Polymorphic --db spec parsing.

Classifies one CLI value into one of: local fasta file, local accession-list
text file, remote s3 fasta, uniprot proteome ID, or uniprot accession.
"""

from __future__ import annotations

import enum
import os
import re
from dataclasses import dataclass

# Longest suffixes first so 'foo.fasta.gz' does not short-circuit on '.fasta'.
_FASTA_EXTS = (".fasta.gz", ".fa.gz", ".fasta", ".fa")
_TXT_EXTS = (".txt",)
_UNIPROT_PROTEOME_RE = re.compile(r"^UP\d{9}$")
_UNIPROT_ACCESSION_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


class DbSpecKind(enum.Enum):
    LOCAL_FASTA = "local_fasta"
    ACCESSION_LIST_FILE = "accession_list_file"
    S3_FASTA = "s3_fasta"
    UNIPROT_PROTEOME = "uniprot_proteome"
    UNIPROT_ACCESSION = "uniprot_accession"


@dataclass(frozen=True)
class DbSpec:
    kind: DbSpecKind
    value: str  # original spec string


def classify_db_spec(spec: str) -> DbSpec:
    """Classify a single --db value. Raises ValueError for unrecognised input."""
    if spec.startswith("s3://"):
        return DbSpec(DbSpecKind.S3_FASTA, spec)

    # Local file? Path-shaped strings get checked first so a stray file named
    # `UP000005640` on disk still resolves as local.
    if os.path.exists(spec):
        lower = spec.lower()
        if any(lower.endswith(ext) for ext in _FASTA_EXTS):
            return DbSpec(DbSpecKind.LOCAL_FASTA, spec)
        if any(lower.endswith(ext) for ext in _TXT_EXTS):
            return DbSpec(DbSpecKind.ACCESSION_LIST_FILE, spec)
        # File exists but unrecognised extension — refuse rather than guess.
        raise ValueError(
            f"unrecognised --db spec: {spec!r}"
            " (file exists but extension is not .fasta/.fa/.txt)"
        )

    if _UNIPROT_PROTEOME_RE.match(spec):
        return DbSpec(DbSpecKind.UNIPROT_PROTEOME, spec)

    if _UNIPROT_ACCESSION_RE.match(spec):
        return DbSpec(DbSpecKind.UNIPROT_ACCESSION, spec)

    raise ValueError(
        f"unrecognised --db spec: {spec!r}"
        " (not s3://, not a local .fasta/.txt, not UP..., not an accession)"
    )
