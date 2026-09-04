#!/usr/bin/env python3
"""Vendor the PSI-MS `spectrum origin type` subtree as Rust.

No stack we depend on carries the PSI-MS controlled vocabulary. mzcv knows `MS`
as a namespace, so an `MS:` CURIE parses and compares, but there is no index of
MS terms and therefore no `is_a` to walk. The one question the mzSpecLib reader
needs answered -- is this spectrum a decoy -- is a subsumption question, so the
hierarchy has to come from somewhere.

It comes from here. This resolves `is_a` against the published ontology and
emits the *answers* as one row per term, rather than emitting the edges for Rust
to walk at runtime. Two consequences worth the trade:

  - Subsumption is derived from the real CV, once, at generation time. A flag in
    the output is not a hand-picked set of accessions; it is whatever the
    ontology says descends from the named root.
  - The reader does a lookup over a static table. No graph, no parse, no
    allocation, and no ordering invariant to state and test.

Run `--check` in CI: one ranged HTTP GET of the first few hundred bytes, compared
against the `DATA_VERSION` in the generated file. Cheap enough to run on every
build, and it turns a stale snapshot from an invisible gap into a failing job.

It reports a release that left this subtree alone, which is most of them. That is
the price of not downloading 1.2 MB per CI run: regenerate, and if the table is
unchanged the only edit is the version string.

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


def render(version: str, names: dict[int, str], parents: dict[int, list[int]]) -> str:
    subtree = sorted(descendants([ROOT], parents))

    # Each closure is computed once. Doing it inside a loop over the subtree ran
    # the fixed-point walk over all ~4000 psi-ms terms per row.
    memberships = {
        name: descendants(roots, parents) & set(subtree) for name, roots, _ in SUBSETS
    }

    body = [
        "//! The PSI-MS `spectrum origin type` subtree: one row per term, with",
        "//! what the term implies.",
        "//!",
        "//! @generated by `scripts/gen_psims_origin_type.py`. Do not edit: run that",
        "//! script instead.",
        "//!",
        "//! Each flag is the transitive closure of `is_a` under a named root,",
        "//! resolved against the ontology at generation time. Membership is derived",
        "//! from the controlled vocabulary rather than hand-picked, and one row per",
        "//! term means the flags cannot disagree with each other.",
        "",
        "/// The `data-version` of the ontology this table was taken from.",
        "///",
        "/// `--check` compares this against the published ontology, so most",
        "/// releases will report a difference that leaves the table unchanged.",
        "/// Regenerating then edits only this line.",
        f'pub const DATA_VERSION: &str = "{version}";',
        "",
        "/// One `spectrum origin type` term, and what it implies for a reader.",
        "#[derive(Debug, Clone, Copy, PartialEq, Eq)]",
        "pub struct OriginType {",
        "    /// The PSI-MS accession, without the `MS:` prefix.",
        "    pub accession: u32,",
    ]
    for name, roots, doc in SUBSETS:
        rendered_roots = ", ".join(f"`MS:{r}|{names.get(r, '?')}`" for r in roots)
        body += [
            f"    /// {doc}",
            "    ///",
            f"    /// Closure of {rendered_roots}.",
            f"    pub {name.lower()}: bool,",
        ]
    body += [
        "}",
        "",
        "/// Every term under `MS:1003072|spectrum origin type`.",
        "///",
        "/// An origin type absent from this table is a term this build has never",
        "/// heard of. The reader reports it rather than guessing, because guessing",
        '/// "target" on an unknown decoy subtype produces an FDR that is wrong',
        "/// instead of absent.",
        "pub const ORIGIN_TYPES: &[OriginType] = &[",
    ]
    for accession in subtree:
        flags = "".join(
            f" {name.lower()}: {str(accession in memberships[name]).lower()},"
            for name, _, _ in SUBSETS
        )
        body.append(
            f"    // MS:{accession} {names[accession]}\n"
            f"    OriginType {{ accession: {accession:_},{flags} }},"
        )
    body += [
        "];",
        "",
        "/// The term, if this build knows it.",
        "///",
        "/// A linear scan: the table is ten rows, so an ordering invariant would",
        "/// cost more to state and test than it saves.",
        "pub fn lookup(accession: u32) -> Option<&'static OriginType> {",
        "    ORIGIN_TYPES.iter().find(|term| term.accession == accession)",
        "}",
        "",
    ]
    return "\n".join(body)


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
            f"Re-run: python3 scripts/{pathlib.Path(__file__).name}\n"
            f"Most releases leave this subtree alone, in which case the only "
            f"change is the version string. If a term did move, check whether it "
            f"changes the reader's behaviour."
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
