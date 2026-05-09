"""Entrapment classification + FDR estimators (Noble et al, FDRBench Table S2).

The fixture's peptide-list inputs (target_peptides, entrapment_peptides) feed
sequence-level set membership classification of result PSMs. Three FDP
estimators are emitted from compute_fdr_curve:

- empirical_fdr_lower      = n_e / (n_e + n_t)               # lower bound
- empirical_fdr_combined   = n_e × (1 + 1/r) / (n_e + n_t)   # avg upper bound
- empirical_fdr_matched    = (n_e + n_p_s_t + 2·n_p_t_s) / (n_e + n_t)
                              # matched estimator, k=1; only when pairing supplied

Where:
- r is the entrapment_ratio recorded on the fixture
- n_p_t_s = entrap hits whose paired target also scored ≥ threshold
            AND entrap_score > paired_target_score
- n_p_s_t = entrap hits whose paired target scored < threshold

The walk filters to is_target=True before accumulating counts (post-competition
target winners only; decoy wins are TDA-style FPs, not entrapment FPs).
"""

from __future__ import annotations

import enum
import re
from pathlib import Path

import polars as pl

_MOD_RE = re.compile(r"\[[^\]]*\]|\([^)]*\)|[0-9.]+")
"""Strip [...], (...), and bare numeric mass shifts from peptide sequences."""


def strip_mods(seq: str) -> str:
    return _MOD_RE.sub("", seq)


def load_peptide_set(path: str | Path) -> set[str]:
    """Load a peptide list .txt file; returns a set of bare AA sequences."""
    out: set[str] = set()
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if s:
            out.add(s)
    return out


def load_pairing(path: str | Path) -> dict[str, str]:
    """Load a pairing.tsv with header `target_peptide\\tentrap_peptide`."""
    lines = Path(path).read_text().splitlines()
    if not lines:
        return {}
    out: dict[str, str] = {}
    # Skip header
    for raw in lines[1:]:
        s = raw.strip()
        if not s:
            continue
        parts = s.split("\t")
        if len(parts) != 2:
            raise ValueError(f"bad pairing row: {raw!r}")
        out[parts[0]] = parts[1]
    return out


class PeptideClass(enum.Enum):
    TARGET = "target"
    ENTRAPMENT = "entrapment"
    SHARED_DROPPED = "shared_dropped"
    UNKNOWN = "unknown"


def classify_peptides(
    results: pl.DataFrame,
    target_peptides: set[str],
    entrapment_peptides: set[str],
) -> pl.DataFrame:
    """Add `class` and `is_entrapment` columns via set membership."""
    if "sequence" not in results.columns:
        raise ValueError("results dataframe missing 'sequence' column")

    shared = target_peptides & entrapment_peptides
    target_only = target_peptides - shared
    entrap_only = entrapment_peptides - shared

    def _classify(seq: str) -> str:
        s = strip_mods(seq)
        # Prefer the stripped form; fall back to original for sets that store
        # sequences without any mod annotations (e.g., peptide-list .txt files
        # whose entries were never annotated, so strip_mods is a no-op on them
        # but may modify the query sequence if it contains numeric characters).
        key = s if (s in target_only or s in entrap_only or s in shared) else seq
        if key in target_only:
            return PeptideClass.TARGET.value
        if key in entrap_only:
            return PeptideClass.ENTRAPMENT.value
        if key in shared:
            return PeptideClass.SHARED_DROPPED.value
        return PeptideClass.UNKNOWN.value

    classes = results["sequence"].map_elements(_classify, return_dtype=pl.Utf8)
    is_entrap = classes == PeptideClass.ENTRAPMENT.value
    return results.with_columns(
        classes.alias("class"),
        is_entrap.alias("is_entrapment"),
    )


def _walk_matched_counts(
    keep: pl.DataFrame,
    pairing: dict[str, str],
) -> tuple[list[int], list[int]]:
    """Walk rows in qvalue order, accumulate (n_p_t_s, n_p_s_t) at each row.

    n_p_t_s = entrap hits whose paired_target ALSO entered the keep set
              before this row AND entrap_score > paired_target_score.
    n_p_s_t = entrap hits whose paired_target has NOT entered the keep set
              by this row (paired target failed the threshold).

    Walks rows in the order they appear in `keep` (already sorted by qvalue).
    """
    # Per-row maps for lookup
    seq_col = keep["sequence"].to_list()
    cls_col = keep["class"].to_list()
    score_col = keep["main_score"].to_list()
    # target peptide -> (row_index, score) when discovered; otherwise absent
    discovered_target_score: dict[str, float] = {}
    pts = 0  # n_p_t_s cumulative
    pst = 0  # n_p_s_t cumulative
    pts_col: list[int] = []
    pst_col: list[int] = []
    # Build target → entrap reverse index from pairing (we look up paired target
    # of an entrap hit; pairing is target -> entrap by convention).
    entrap_to_target = {e: t for t, e in pairing.items()}
    for i in range(len(seq_col)):
        seq = seq_col[i]
        cls = cls_col[i]
        score = float(score_col[i])
        if cls == PeptideClass.TARGET.value:
            discovered_target_score[seq] = score
        elif cls == PeptideClass.ENTRAPMENT.value:
            paired_target = entrap_to_target.get(seq)
            if paired_target is None:
                # Entrap not in pairing dict → treat as if its paired target
                # was not discovered; counts as n_p_s_t per the "paired target
                # didn't reach threshold" branch.
                pst += 1
            elif paired_target in discovered_target_score:
                # Paired target already in keep set; compare scores
                tscore = discovered_target_score[paired_target]
                if score > tscore:
                    pts += 1
                # else: entrap < target; not an upper-bound contribution
            else:
                # Paired target not discovered yet → entrap above s, target below
                pst += 1
        pts_col.append(pts)
        pst_col.append(pst)
    return pts_col, pst_col


def compute_fdr_curve(
    classified: pl.DataFrame,
    ratio: float,
    pairing: dict[str, str] | None = None,
) -> pl.DataFrame:
    """Sort by qvalue (ascending), accumulate target/entrap counts, return curve.

    Filters to is_target=True before walking (decoys are TDA FPs, not entrapment FPs).
    Only rows whose class is TARGET or ENTRAPMENT contribute.

    Emits columns: n_target, n_entrap, empirical_fdr_lower, empirical_fdr_combined.
    If `pairing` provided AND `main_score` column present, also emits
    empirical_fdr_matched.
    """
    if "qvalue" not in classified.columns:
        raise ValueError("classified dataframe missing 'qvalue' column")
    if "is_target" not in classified.columns:
        raise ValueError("classified dataframe missing 'is_target' column")
    if ratio < 1.0:
        raise ValueError(f"ratio must be >= 1.0, got {ratio}")

    keep = classified.filter(
        pl.col("is_target")
        & pl.col("class").is_in(
            [PeptideClass.TARGET.value, PeptideClass.ENTRAPMENT.value]
        )
    ).sort("qvalue")

    n_target = (
        (keep["class"] == PeptideClass.TARGET.value).cast(pl.UInt32).cum_sum()
    )
    n_entrap = (
        (keep["class"] == PeptideClass.ENTRAPMENT.value).cast(pl.UInt32).cum_sum()
    )
    n_total = (n_target + n_entrap).cast(pl.Float64)
    lower = (n_entrap.cast(pl.Float64) / n_total).fill_nan(0.0)
    combined = (n_entrap.cast(pl.Float64) * (1.0 + 1.0 / ratio) / n_total).fill_nan(
        0.0
    )

    out = keep.with_columns(
        n_target.alias("n_target"),
        n_entrap.alias("n_entrap"),
        lower.alias("empirical_fdr_lower"),
        combined.alias("empirical_fdr_combined"),
    )

    if pairing is not None and "main_score" in keep.columns:
        pts_col, pst_col = _walk_matched_counts(out, pairing)
        out = out.with_columns(
            pl.Series("n_p_t_s", pts_col, dtype=pl.UInt32),
            pl.Series("n_p_s_t", pst_col, dtype=pl.UInt32),
        )
        matched_num = (
            n_entrap.cast(pl.Float64)
            + pl.Series("n_p_s_t", pst_col).cast(pl.Float64)
            + 2.0 * pl.Series("n_p_t_s", pts_col).cast(pl.Float64)
        )
        matched = (matched_num / n_total).fill_nan(0.0)
        out = out.with_columns(matched.alias("empirical_fdr_matched"))

    return out


def plot_fdr_curve(
    curve: pl.DataFrame,
    output_path: str | Path,
    title: str = "Reported q-value vs empirical entrapment FDR",
    xlim: float | None = None,
) -> None:
    """Render curves to a PNG. Plots all available estimators."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)
    if "empirical_fdr_combined" in curve.columns:
        ax.plot(
            curve["qvalue"], curve["empirical_fdr_combined"], lw=1.5,
            label="combined (avg upper bound)",
        )
    if "empirical_fdr_matched" in curve.columns:
        ax.plot(
            curve["qvalue"], curve["empirical_fdr_matched"], lw=1.5,
            label="matched (k=1, avg upper bound)",
        )
    if "empirical_fdr_lower" in curve.columns:
        ax.plot(
            curve["qvalue"], curve["empirical_fdr_lower"], lw=1.0,
            ls=":", alpha=0.7, label="lower bound",
        )

    if xlim is None:
        _max_qv = curve["qvalue"].cast(pl.Float64).max()
        lim: float = max(0.05, _max_qv if isinstance(_max_qv, float) else 0.05)
    else:
        lim = float(xlim)
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


def plot_score_histogram(
    classified: pl.DataFrame,
    output_path: str | Path,
    title: str = "main_score by class × is_target",
) -> None:
    """Render a 4-group main_score histogram (log-x, log-y count)."""
    import matplotlib
    import numpy as np

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if "main_score" not in classified.columns:
        raise ValueError("classified dataframe missing 'main_score' column")

    df = classified
    groups = {
        "target × class=target": df.filter(
            pl.col("is_target") & (pl.col("class") == "target")
        ),
        "target × class=entrap": df.filter(
            pl.col("is_target") & (pl.col("class") == "entrapment")
        ),
        "decoy × class=target": df.filter(
            ~pl.col("is_target") & (pl.col("class") == "target")
        ),
        "decoy × class=entrap": df.filter(
            ~pl.col("is_target") & (pl.col("class") == "entrapment")
        ),
    }
    all_scores = df["main_score"].to_numpy()
    all_scores = all_scores[np.isfinite(all_scores) & (all_scores > 0)]
    if len(all_scores) == 0:
        # Empty plot rather than crash
        fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
        ax.set_title(title + " (no positive scores)")
        fig.savefig(output_path)
        plt.close(fig)
        return
    bins = np.logspace(np.log10(all_scores.min()), np.log10(all_scores.max()), 80)

    fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
    colors = {
        "target × class=target": "C0",
        "target × class=entrap": "C1",
        "decoy × class=target":  "C2",
        "decoy × class=entrap":  "C3",
    }
    for label, sub in groups.items():
        s = sub["main_score"].to_numpy()
        s = s[np.isfinite(s) & (s > 0)]
        if len(s) == 0:
            continue
        ax.hist(
            s, bins=bins, alpha=0.6, label=f"{label} (n={len(s)})",
            color=colors[label], histtype="step", linewidth=1.5,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("main_score (log)")
    ax.set_ylabel("count (log)")
    ax.legend(loc="best", fontsize=9)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def _scalar_at_q(
    curve: pl.DataFrame, q_threshold: float, suffix: str
) -> dict[str, float | int]:
    sub = curve.filter(pl.col("qvalue") <= q_threshold)
    if sub.height == 0:
        out: dict[str, float | int] = {
            f"entrap/n_target_at_{suffix}": 0,
            f"entrap/n_entrap_at_{suffix}": 0,
            f"entrap/empirical_fdr_lower_at_{suffix}": 0.0,
            f"entrap/empirical_fdr_combined_at_{suffix}": 0.0,
        }
        if "empirical_fdr_matched" in curve.columns:
            out[f"entrap/empirical_fdr_matched_at_{suffix}"] = 0.0
        return out
    last = sub.row(-1, named=True)
    out = {
        f"entrap/n_target_at_{suffix}": int(last["n_target"]),
        f"entrap/n_entrap_at_{suffix}": int(last["n_entrap"]),
        f"entrap/empirical_fdr_lower_at_{suffix}": float(last["empirical_fdr_lower"]),
        f"entrap/empirical_fdr_combined_at_{suffix}": float(
            last["empirical_fdr_combined"]
        ),
    }
    if "empirical_fdr_matched" in curve.columns:
        out[f"entrap/empirical_fdr_matched_at_{suffix}"] = float(
            last["empirical_fdr_matched"]
        )
    return out


def analyse(
    results_parquet: str | Path,
    target_peptides: str | Path,
    entrapment_peptides: str | Path,
    ratio: float,
    pairing_path: str | Path | None,
    out_parquet: str | Path,
    out_fdr_plot: str | Path,
    out_hist_plot: str | Path,
    title: str = "Reported q-value vs empirical entrapment FDR",
) -> dict[str, float | int]:
    """End-to-end: classify → curve → write parquet + 2 plots → return scalars."""
    results = pl.read_parquet(results_parquet)
    target_set = load_peptide_set(target_peptides)
    entrap_set = load_peptide_set(entrapment_peptides)
    classified = classify_peptides(results, target_set, entrap_set)
    Path(out_parquet).parent.mkdir(parents=True, exist_ok=True)
    classified.write_parquet(out_parquet)

    pairing = None
    if pairing_path is not None:
        pairing = load_pairing(pairing_path)

    curve = compute_fdr_curve(classified, ratio=ratio, pairing=pairing)
    Path(out_fdr_plot).parent.mkdir(parents=True, exist_ok=True)
    plot_fdr_curve(curve, out_fdr_plot, title=title)

    Path(out_hist_plot).parent.mkdir(parents=True, exist_ok=True)
    plot_score_histogram(classified, out_hist_plot, title=f"{title} (score hist)")

    scalars: dict[str, float | int] = {"entrap/ratio": float(ratio)}
    scalars.update(_scalar_at_q(curve, 0.01, "q01"))
    scalars.update(_scalar_at_q(curve, 0.05, "q05"))
    return scalars
