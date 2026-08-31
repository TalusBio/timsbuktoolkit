#!/usr/bin/env python3
"""Vendor the PSI-MS `spectrum origin type` subtree as Rust.

No stack we depend on carries the PSI-MS controlled vocabulary. mzcv knows `MS`
as a namespace, so an `MS:` CURIE parses and compares, but there is no index of
MS terms and therefore no `is_a` to walk. The one question the mzSpecLib reader
needs answered -- is this spectrum a decoy -- is a subsumption question, so the
hierarchy has to come from somewhere.

It comes from here. This resolves `is_a` against the published ontology and
emits the *answers* as sorted accession lists, rather than emitting the edges
for Rust to walk at runtime. Two consequences worth the trade:

  - Subsumption is derived from the real CV, once, at generation time. A list in
    the output is not a hand-picked set of accessions; it is whatever the
    ontology says descends from the named root.
  - The reader does a binary search over a static slice. No graph, no parse, no
    allocation.

Run `--check` in CI: one ranged HTTP GET of the first few hundred bytes, and a
comparison against the `DATA_VERSION` recorded in the generated file. That turns
a stale snapshot from an invisible gap into a failing job.

Usage:
    python3 scripts/gen_psims_origin_type.py            # regenerate
    python3 scripts/gen_psims_origin_type.py --check    # CI freshness check

Standard library only, so CI can call `python3` directly.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys
import urllib.request

OBO_URL = "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"

OUT_PATH = (
    pathlib.Path(__file__).resolve().parent.parent
    / "rust/timsquery/src/serde/psims_origin_type.rs"
)

# `MS:1003072|spectrum origin type`. Every term below is reachable from it.
ROOT = 1003072

# The lists to emit, as `(rust_name, root_accession, doc)`. Each is the
# transitive closure of `is_a` under its root, intersected with the ROOT
# subtree.
#
# `decoy spectrum` is the one the reader's FDR depends on. The two theoretical
# roots are separate because there is no single parent term meaning "these
# masses are calculated": `predicted spectrum` says the whole spectrum was, and
# `selected fragment theoretical m/z observed intensity spectrum` says it of the
# m/z alone.
SUBSETS = [
    (
        "DECOY",
        [1003192],
        "Spectra whose analyte is deliberately wrong, so an FDR can be estimated.",
    ),
    (
        "THEORETICAL_MZ",
        [1003074, 1003424],
        "Spectra whose peak m/z values are calculated rather than measured.",
    ),
]

ID_RE = re.compile(r"^id: MS:(\d+)$", re.MULTILINE)
NAME_RE = re.compile(r"^name: (.+)$", re.MULTILINE)
IS_A_RE = re.compile(r"^is_a: MS:(\d+)", re.MULTILINE)
VERSION_RE = re.compile(r"^data-version: (.+)$", re.MULTILINE)


def fetch(url: str, byte_range: str | None = None) -> str:
    request = urllib.request.Request(url)
    if byte_range is not None:
        request.add_header("Range", f"bytes={byte_range}")
    with urllib.request.urlopen(request) as response:
        return response.read().decode("utf-8", errors="replace")


def published_version() -> str:
    """The ontology's own `data-version`, from its first few hundred bytes."""
    head = fetch(OBO_URL, byte_range="0-400")
    match = VERSION_RE.search(head)
    if match is None:
        sys.exit("psi-ms.obo carries no data-version line; cannot check freshness")
    return match.group(1).strip()


def vendored_version() -> str:
    if not OUT_PATH.exists():
        sys.exit(f"{OUT_PATH} does not exist; run this script without --check")
    match = re.search(r'DATA_VERSION: &str = "([^"]+)"', OUT_PATH.read_text())
    if match is None:
        sys.exit(f"{OUT_PATH} carries no DATA_VERSION; regenerate it")
    return match.group(1)


def parse_terms(obo: str) -> tuple[dict[int, str], dict[int, list[int]]]:
    """`(name, is_a parents)` per MS accession, over every `[Term]` stanza."""
    names: dict[int, str] = {}
    parents: dict[int, list[int]] = {}
    for stanza in obo.split("\n\n"):
        if not stanza.startswith("[Term]"):
            continue
        accession = ID_RE.search(stanza)
        name = NAME_RE.search(stanza)
        if accession is None or name is None:
            continue
        key = int(accession.group(1))
        names[key] = name.group(1).strip()
        parents[key] = [int(m) for m in IS_A_RE.findall(stanza)]
    return names, parents


def descendants(roots: list[int], parents: dict[int, list[int]]) -> set[int]:
    """Every accession reaching any of `roots` through `is_a`, roots included.

    Iterated to a fixed point rather than recursed: the graph is a DAG, a term
    may have several parents, and a malformed ontology could contain a cycle.
    """
    closure = set(roots)
    changed = True
    while changed:
        changed = False
        for accession, own_parents in parents.items():
            if accession in closure:
                continue
            if any(parent in closure for parent in own_parents):
                closure.add(accession)
                changed = True
    return closure


def rust_list(accessions: list[int]) -> str:
    return "".join(f"    {a:_},\n" for a in accessions)


def render(version: str, names: dict[int, str], parents: dict[int, list[int]]) -> str:
    subtree = sorted(descendants([ROOT], parents))

    body = [
        "//! The PSI-MS `spectrum origin type` subtree, as sorted accession lists.",
        "//!",
        "//! @generated by `scripts/gen_psims_origin_type.py`. Do not edit: run that",
        "//! script instead. `--check` compares the version below against the",
        "//! published ontology, so a stale snapshot fails CI rather than going",
        "//! unnoticed.",
        "//!",
        "//! Every list is the transitive closure of `is_a` under a named root,",
        "//! resolved against the ontology at generation time. So membership is",
        "//! derived from the controlled vocabulary rather than hand-picked, and the",
        "//! reader needs only a binary search.",
        "",
        "/// The `data-version` of the ontology these lists were taken from.",
        f'pub const DATA_VERSION: &str = "{version}";',
        "",
        "/// Every accession under `MS:1003072|spectrum origin type`, sorted.",
        "///",
        "/// A `spectrum origin type` value outside this list is a term this build",
        "/// has never heard of. The reader reports it rather than guessing, because",
        "/// guessing \"target\" on an unknown decoy subtype produces an FDR that is",
        "/// wrong instead of absent.",
        "pub const ALL: &[u32] = &[",
        rust_list(subtree).rstrip("\n"),
        "];",
        "",
    ]

    for name, roots, doc in SUBSETS:
        members = sorted(descendants(roots, parents) & set(subtree))
        rendered_roots = ", ".join(
            f"`MS:{r}|{names.get(r, '?')}`" for r in roots
        )
        body += [
            f"/// {doc}",
            "///",
            f"/// Closure of {rendered_roots}.",
            f"pub const {name}: &[u32] = &[",
            rust_list(members).rstrip("\n"),
            "];",
            "",
        ]

    # A comment block naming what landed, so a reviewer can see the membership
    # without re-running the script or opening the OBO.
    body += ["// Terms in this snapshot:"]
    for accession in subtree:
        tags = [name for name, roots, _ in SUBSETS
                if accession in (descendants(roots, parents) & set(subtree))]
        suffix = f"  [{', '.join(tags)}]" if tags else ""
        body.append(f"//   MS:{accession}  {names[accession]}{suffix}")
    body.append("")

    return "\n".join(body)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare the vendored data-version against the published one",
    )
    args = parser.parse_args()

    if args.check:
        published = published_version()
        vendored = vendored_version()
        if published == vendored:
            print(f"psi-ms data-version {vendored} is current")
            return
        sys.exit(
            f"psi-ms.obo has moved to data-version {published}; "
            f"{OUT_PATH.name} was generated from {vendored}.\n"
            f"Re-run: python3 scripts/gen_psims_origin_type.py\n"
            f"Then check whether any new term changes the reader's behaviour."
        )

    obo = fetch(OBO_URL)
    version = published_version()
    names, parents = parse_terms(obo)
    if ROOT not in names:
        sys.exit(f"MS:{ROOT} is absent from psi-ms.obo; the root has been renamed")
    OUT_PATH.write_text(render(version, names, parents))
    print(f"wrote {OUT_PATH} from psi-ms data-version {version}")


if __name__ == "__main__":
    main()
