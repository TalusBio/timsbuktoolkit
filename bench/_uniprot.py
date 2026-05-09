"""Uniprot REST helpers.

We use the *stream* endpoint which returns full FASTA payloads without
pagination. See https://www.uniprot.org/help/api_queries.
"""

from __future__ import annotations

from typing import Iterable

import requests
from loguru import logger

_BASE = "https://rest.uniprot.org/uniprotkb/stream"
_TIMEOUT = 120
_BATCH_SIZE = 100  # accession IDs per query — keeps URL well under the ~8KB limit


def _get(params: dict[str, str]) -> str:
    logger.info("uniprot GET {} {}", _BASE, params)
    r = requests.get(_BASE, params=params, timeout=_TIMEOUT)
    r.raise_for_status()
    return r.text


def fetch_proteome(proteome_id: str, reviewed_only: bool = True) -> str:
    """Fetch a uniprot proteome (e.g. UP000005640) as FASTA text.

    Defaults to Swiss-Prot only (`reviewed:true`) to keep search spaces
    tractable; full proteome (incl. TrEMBL) is rarely what bench fixtures
    actually want and is much slower to build a speclib over. Pass
    `reviewed_only=False` for the unfiltered set.
    """
    query = f"proteome:{proteome_id}"
    if reviewed_only:
        query += " AND reviewed:true"
    return _get({"query": query, "format": "fasta"})


def fetch_accession_batch(accessions: Iterable[str]) -> str:
    """Fetch FASTA for a list of accessions. Chunks under URL-length limits."""
    accs = list(accessions)
    out: list[str] = []
    for i in range(0, len(accs), _BATCH_SIZE):
        chunk = accs[i : i + _BATCH_SIZE]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        out.append(_get({"query": query, "format": "fasta"}))
    return "".join(out)
