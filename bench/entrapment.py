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


def count_kmers(fasta_path: str | Path, k: int = 7) -> set[str]:
    """Return the set of unique k-mers across all proteins in a FASTA.

    Window slides over each protein sequence; sequences shorter than `k`
    contribute nothing.
    """
    proteins = parse_fasta(fasta_path)
    out: set[str] = set()
    for seq in proteins.values():
        if len(seq) < k:
            continue
        for i in range(len(seq) - k + 1):
            out.add(seq[i : i + k])
    return out


def kmer_normalization_factor(
    target_fasta: str | Path,
    entrapment_fasta: str | Path,
    k: int = 7,
) -> float:
    """Compute |T_unique_kmers| / |E_unique_kmers| after dropping shared k-mers.

    Used to rescale entrapment hit counts so empirical FDR is comparable
    across proteomes of unequal search-space size. Clamped at 1.0 in the
    denominator to avoid div-by-zero on tiny / empty entrapment fastas.
    """
    t = count_kmers(target_fasta, k=k)
    e = count_kmers(entrapment_fasta, k=k)
    shared = t & e
    t_only = t - shared
    e_only = e - shared
    return len(t_only) / max(1, len(e_only))


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


def compute_fdr_curve(
    classified: pl.DataFrame,
    normalization_factor: float = 1.0,
) -> pl.DataFrame:
    """Sort by qvalue, accumulate target/entrapment counts, return curve.

    Rows whose class is SHARED_DROPPED or UNKNOWN are excluded from both
    numerator and denominator.

    `normalization_factor` rescales entrapment counts to compensate for
    differences in target vs entrapment search-space size (e.g., from
    `kmer_normalization_factor`). With factor=1.0 (default), `empirical_fdr_norm`
    equals `empirical_fdr_raw`.
    """
    if "qvalue" not in classified.columns:
        raise ValueError("classified dataframe missing 'qvalue' column")

    keep = classified.filter(
        pl.col("class").is_in(
            [PeptideClass.TARGET.value, PeptideClass.ENTRAPMENT.value]
        )
    ).sort("qvalue")

    n_target = (
        (keep["class"] == PeptideClass.TARGET.value).cast(pl.UInt32).cum_sum()
    )
    n_entrap = (
        (keep["class"] == PeptideClass.ENTRAPMENT.value).cast(pl.UInt32).cum_sum()
    )
    raw_fdr = (n_entrap.cast(pl.Float64) / (n_target + n_entrap)).fill_nan(0.0)
    n_entrap_norm = n_entrap.cast(pl.Float64) * normalization_factor
    norm_fdr = (
        n_entrap_norm / (n_target.cast(pl.Float64) + n_entrap_norm)
    ).fill_nan(0.0)

    return keep.with_columns(
        n_target.alias("n_target"),
        n_entrap.alias("n_entrap"),
        n_entrap_norm.alias("n_entrap_norm"),
        raw_fdr.alias("empirical_fdr_raw"),
        norm_fdr.alias("empirical_fdr_norm"),
    )


def plot_fdr_curve(
    curve: pl.DataFrame,
    output_path: str | Path,
    title: str = "Reported q-value vs empirical entrapment FDR",
) -> None:
    """Render a FDR-vs-qvalue plot to a PNG file. Plots the kmer-normalized
    curve (primary) and the raw curve (faded) for comparison."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)
    ax.plot(
        curve["qvalue"],
        curve["empirical_fdr_norm"],
        lw=1.5,
        label="empirical (kmer-normalized)",
    )
    ax.plot(
        curve["qvalue"],
        curve["empirical_fdr_raw"],
        lw=1.0,
        ls=":",
        alpha=0.6,
        label="empirical (raw, unnormalized)",
    )
    _max_qv = curve["qvalue"].cast(pl.Float64).max()
    lim: float = max(0.05, _max_qv if isinstance(_max_qv, float) else 0.05)
    ax.plot([0, lim], [0, lim], color="grey", ls="--", lw=1.0, label="y=x")
    ax.set_xlabel("reported q-value")
    ax.set_ylabel("empirical FDR")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_title(title)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def _scalar_at_q(
    curve: pl.DataFrame, q_threshold: float, suffix: str
) -> dict[str, float | int]:
    """Read off n_target / n_entrap / empirical FDR (raw + norm) at q <= threshold."""
    sub = curve.filter(pl.col("qvalue") <= q_threshold)
    if sub.height == 0:
        return {
            f"entrap/n_target_at_{suffix}": 0,
            f"entrap/n_entrap_at_{suffix}": 0,
            f"entrap/empirical_fdr_raw_at_{suffix}": 0.0,
            f"entrap/empirical_fdr_norm_at_{suffix}": 0.0,
        }
    last = sub.row(-1, named=True)
    return {
        f"entrap/n_target_at_{suffix}": int(last["n_target"]),
        f"entrap/n_entrap_at_{suffix}": int(last["n_entrap"]),
        f"entrap/empirical_fdr_raw_at_{suffix}": float(last["empirical_fdr_raw"]),
        f"entrap/empirical_fdr_norm_at_{suffix}": float(last["empirical_fdr_norm"]),
    }


def analyse(
    results_parquet: str | Path,
    target_fasta: str | Path,
    entrapment_fasta: str | Path,
    out_parquet: str | Path,
    out_plot: str | Path,
    title: str = "Reported q-value vs empirical entrapment FDR",
    kmer_k: int = 7,
) -> dict[str, float | int]:
    """End-to-end: classify → kmer norm → FDR walk → write parquet+plot → scalars."""
    results = pl.read_parquet(results_parquet)
    classified = classify_peptides(results, target_fasta, entrapment_fasta)
    Path(out_parquet).parent.mkdir(parents=True, exist_ok=True)
    classified.write_parquet(out_parquet)

    factor = kmer_normalization_factor(target_fasta, entrapment_fasta, k=kmer_k)
    curve = compute_fdr_curve(classified, normalization_factor=factor)
    Path(out_plot).parent.mkdir(parents=True, exist_ok=True)
    plot_fdr_curve(curve, out_plot, title=title)

    scalars: dict[str, float | int] = {"entrap/normalization_factor": float(factor)}
    scalars.update(_scalar_at_q(curve, 0.01, "q01"))
    scalars.update(_scalar_at_q(curve, 0.05, "q05"))
    return scalars
