# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "polars",
#     "pyarrow",
#     "matplotlib",
#     "numpy",
# ]
# ///
"""Diagnostic plots for a timsseek results directory.

Reads `<results-dir>/<run>/results.parquet`, filters at the given q-value
threshold, writes PNGs to an output directory.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl


def load_results(
    results_dir: Path, qvalue_max: float
) -> tuple[pl.DataFrame, pl.DataFrame]:
    parquets = list(results_dir.rglob("results.parquet"))
    if not parquets:
        raise FileNotFoundError(f"no results.parquet under {results_dir}")
    raw = pl.concat([pl.read_parquet(p) for p in parquets], how="vertical_relaxed")
    passing = raw.filter(pl.col("qvalue") <= qvalue_max)
    return raw, passing


def _finite(df: pl.DataFrame, col: str) -> np.ndarray:
    arr = df[col].to_numpy()
    return arr[np.isfinite(arr)]


def plot_mass_errors(df: pl.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, col, title in [
        (axes[0], "ms1_mz_error_0", "MS1 primary m/z error"),
        (axes[1], "ms2_mz_error_0", "MS2 primary m/z error"),
    ]:
        vals = _finite(df, col)
        ax.hist(vals, bins=80, color="#3366cc", alpha=0.85)
        ax.axvline(
            np.median(vals),
            color="k",
            linestyle="--",
            linewidth=1,
            label=f"median={np.median(vals):.2f}",
        )
        ax.set_xlabel("ppm")
        ax.set_title(f"{title} (n={vals.size})")
        ax.legend()
    axes[0].set_ylabel("count")
    fig.tight_layout()
    fig.savefig(out / "mass_errors.png", dpi=150)
    plt.close(fig)


def plot_mobility_error(df: pl.DataFrame, out: Path) -> None:
    vals = _finite(df, "ms1_mobility_error_0")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(vals, bins=80, color="#cc6633", alpha=0.85)
    ax.axvline(
        np.median(vals),
        color="k",
        linestyle="--",
        linewidth=1,
        label=f"median={np.median(vals):.3f}",
    )
    ax.set_xlabel("mobility error (primary MS1)")
    ax.set_ylabel("count")
    ax.set_title(f"MS1 primary mobility error (n={vals.size})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "mobility_error.png", dpi=150)
    plt.close(fig)


def plot_rt_calibration(df: pl.DataFrame, out: Path) -> None:
    lib = df["library_rt"].to_numpy()
    obs = df["obs_rt_seconds"].to_numpy()
    mask = np.isfinite(lib) & np.isfinite(obs)
    lib, obs = lib[mask], obs[mask]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(lib, obs, s=4, alpha=0.3, color="#3366cc", rasterized=True)
    ax.set_xlabel("library iRT")
    ax.set_ylabel("observed RT (s)")
    ax.set_title(f"Library RT vs observed RT (n={lib.size})")
    fig.tight_layout()
    fig.savefig(out / "rt_calibration.png", dpi=150)
    plt.close(fig)


def plot_rt_vs_mass_error(df: pl.DataFrame, out: Path) -> None:
    rt = df["obs_rt_seconds"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, col, title in [
        (axes[0], "ms1_mz_error_0", "MS1"),
        (axes[1], "ms2_mz_error_0", "MS2"),
    ]:
        err = df[col].to_numpy()
        mask = np.isfinite(rt) & np.isfinite(err)
        ax.scatter(
            rt[mask], err[mask], s=3, alpha=0.25, color="#3366cc", rasterized=True
        )
        ax.axhline(0, color="k", linewidth=0.5)
        ax.set_xlabel("observed RT (s)")
        ax.set_title(f"{title} m/z error vs RT (n={mask.sum()})")
    axes[0].set_ylabel("ppm")
    fig.tight_layout()
    fig.savefig(out / "rt_vs_mass_error.png", dpi=150)
    plt.close(fig)


def _needs_log_scale(vals: np.ndarray) -> bool:
    pos = vals[vals > 0]
    if pos.size < 10:
        return False
    lo, hi = np.quantile(pos, [0.01, 0.99])
    return lo > 0 and hi / lo > 1000


def plot_score_target_decoy(raw: pl.DataFrame, out: Path) -> None:
    targets = _finite(raw.filter(pl.col("is_target")), "main_score")
    decoys = _finite(raw.filter(~pl.col("is_target")), "main_score")
    combined = np.concatenate([targets, decoys]) if decoys.size else targets
    log_x = _needs_log_scale(combined)

    if log_x:
        targets = targets[targets > 0]
        decoys = decoys[decoys > 0]
        lo = min(targets.min(), decoys.min()) if decoys.size else targets.min()
        hi = max(targets.max(), decoys.max()) if decoys.size else targets.max()
        bins = np.geomspace(lo, hi, 80)
    else:
        lo = min(targets.min(), decoys.min()) if decoys.size else targets.min()
        hi = max(targets.max(), decoys.max()) if decoys.size else targets.max()
        bins = np.linspace(lo, hi, 80)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(
        targets,
        bins=bins,
        alpha=0.7,
        label=f"target (n={targets.size})",
        color="#3366cc",
    )
    ax.hist(
        decoys, bins=bins, alpha=0.7, label=f"decoy (n={decoys.size})", color="#cc3333"
    )
    ax.set_xlabel("main_score" + (" (log)" if log_x else ""))
    ax.set_ylabel("count")
    ax.set_title("Target vs decoy main_score (unfiltered)")
    if log_x:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "score_target_decoy.png", dpi=150)
    plt.close(fig)


def plot_qvalue_curve(raw: pl.DataFrame, out: Path) -> None:
    q = raw.filter(pl.col("is_target"))["qvalue"].to_numpy()
    q = np.sort(q[np.isfinite(q)])
    ids = np.arange(1, q.size + 1)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(q, ids, color="#3366cc")
    for thr in (0.01, 0.05):
        n = int((q <= thr).sum())
        ax.axvline(thr, color="k", linestyle="--", linewidth=0.8)
        ax.text(thr, n, f"  q<={thr}: {n}", va="bottom")
    ax.set_xlabel("q-value")
    ax.set_ylabel("cumulative target IDs")
    ax.set_xscale("log")
    ax.set_title("ID curve (targets)")
    fig.tight_layout()
    fig.savefig(out / "qvalue_curve.png", dpi=150)
    plt.close(fig)


def plot_ids_by_charge(df: pl.DataFrame, out: Path) -> None:
    counts = (
        df
        .filter(pl.col("is_target"))
        .group_by("precursor_charge")
        .len()
        .sort("precursor_charge")
    )
    charges = counts["precursor_charge"].to_list()
    n = counts["len"].to_list()
    positions = list(range(len(charges)))
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(positions, n, color="#3366cc")
    ax.set_xticks(positions, [str(c) for c in charges])
    for x, y in zip(positions, n):
        ax.text(x, y, str(y), ha="center", va="bottom")
    ax.set_xlabel("precursor charge")
    ax.set_ylabel("target IDs")
    ax.set_title("Target IDs by charge")
    fig.tight_layout()
    fig.savefig(out / "ids_by_charge.png", dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "results_dir",
        type=Path,
        help="timsseek output dir (searched recursively for results.parquet)",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        default=None,
        help="output dir for PNGs (default: <results_dir>/diagnostic_plots)",
    )
    parser.add_argument(
        "-q",
        "--qvalue",
        type=float,
        default=0.01,
        help="q-value threshold for filtered plots (default: 0.01)",
    )
    args = parser.parse_args()

    out = args.out or (args.results_dir / "diagnostic_plots")
    out.mkdir(parents=True, exist_ok=True)

    raw, passing = load_results(args.results_dir, args.qvalue)
    print(f"loaded {raw.height} rows, {passing.height} pass q<={args.qvalue}")

    plot_mass_errors(passing, out)
    plot_mobility_error(passing, out)
    plot_rt_calibration(passing, out)
    plot_rt_vs_mass_error(passing, out)
    plot_ids_by_charge(passing, out)
    plot_score_target_decoy(raw, out)
    plot_qvalue_curve(raw, out)

    print(f"wrote plots to {out}")


if __name__ == "__main__":
    main()
