"""Polymorphic --db spec parsing and resolution.

Classifies one CLI value into one of: local fasta file, local accession-list
text file, remote s3 fasta, uniprot proteome ID, or uniprot accession.
Also provides resolve_dbs() to turn a list of specs into a merged FASTA file.
"""

from __future__ import annotations

import enum
import gzip
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path

from bench._s3 import s3_download_file
from bench._uniprot import fetch_accession_batch, fetch_proteome

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
    SHUFFLED = "shuffled"


@dataclass(frozen=True)
class DbSpec:
    kind: DbSpecKind
    value: str  # original spec string


def classify_db_spec(spec: str) -> DbSpec:
    """Classify a single --db value. Raises ValueError for unrecognised input."""
    if spec == "SHUFFLED":
        return DbSpec(DbSpecKind.SHUFFLED, spec)

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


# ---------------------------------------------------------------------------
# Resolution helpers
# ---------------------------------------------------------------------------


def _read_local_fasta_text(path: str) -> str:
    p = Path(path)
    if str(p).lower().endswith(".gz"):
        with gzip.open(p, "rt") as f:
            return f.read()
    return p.read_text()


def _read_accession_list_file(path: str) -> list[str]:
    accs: list[str] = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line:
            accs.append(line)
    return accs


def resolve_dbs(specs: list[str], output_path: Path) -> None:
    """Resolve a list of --db specs into a concatenated FASTA at output_path.

    Individual UniProt accessions (UNIPROT_ACCESSION kind) are coalesced into
    a single batched fetch call (one HTTP round-trip). Accession-list files
    each produce their own batch call. Other spec kinds are fetched
    independently and appended in CLI order.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    classified = [classify_db_spec(raw) for raw in specs]

    for spec in classified:
        if spec.kind is DbSpecKind.SHUFFLED:
            raise ValueError(
                "SHUFFLED is not resolvable to a fasta"
                " — push_fixture handles it directly"
            )

    # Coalesce all bare UNIPROT_ACCESSION specs into one batch.
    bare_accessions: list[str] = [
        s.value for s in classified if s.kind is DbSpecKind.UNIPROT_ACCESSION
    ]

    fragments: list[str] = []
    pending_accessions_flushed = False

    for spec in classified:
        if spec.kind is DbSpecKind.UNIPROT_ACCESSION:
            # Flush the whole accession batch on the first encounter.
            if not pending_accessions_flushed:
                fragments.append(fetch_accession_batch(bare_accessions))
                pending_accessions_flushed = True
            # Subsequent accessions are already included in the batch above.
        elif spec.kind is DbSpecKind.LOCAL_FASTA:
            fragments.append(_read_local_fasta_text(spec.value))
        elif spec.kind is DbSpecKind.S3_FASTA:
            with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
                tmp_path = tmp.name
            try:
                s3_download_file(spec.value, tmp_path)
                fragments.append(_read_local_fasta_text(tmp_path))
            finally:
                Path(tmp_path).unlink(missing_ok=True)
        elif spec.kind is DbSpecKind.UNIPROT_PROTEOME:
            fragments.append(fetch_proteome(spec.value))
        elif spec.kind is DbSpecKind.ACCESSION_LIST_FILE:
            accs = _read_accession_list_file(spec.value)
            fragments.append(fetch_accession_batch(accs))
        else:
            raise AssertionError(f"unhandled kind {spec.kind}")

    with output_path.open("w") as f:
        for chunk in fragments:
            if not chunk.endswith("\n"):
                chunk = chunk + "\n"
            f.write(chunk)
