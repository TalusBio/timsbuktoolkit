from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from rich.pretty import pprint


def plot_importances(importances: dict[str, list[float]]):
    # Lollipop plot showing the importance of each feature
    fig, ax = plt.subplots(figsize=(6, 8))
    rev_imps = [(k, v) for k, v in importances.items()]
    rev_imps = rev_imps[::-1]
    for feature, importance in rev_imps:
        expanded_feat = [feature] * len(importance)
        # Stem of the lollipop
        ax.hlines(expanded_feat, 0, importance)

        # Dot at the top of the lollipop
        ax.scatter(importance, expanded_feat, c="k", marker="o")

    # Rotate the labels
    ax.tick_params(axis="x", rotation=45)

    ax.set_ylabel("Importance ('gain' as defined by xgboost)")
    ax.set_xlabel("Feature")
    # square root scale the x axis
    ax.set_xscale("log")
    fig.tight_layout()
    return fig


def main_score_hist(df: pl.LazyFrame, output_dir: Path):
    pprint("Plotting main scores")
    scores_df = df.select(["is_target", "main_score", "precursor_charge"]).collect()
    target_scores = scores_df.filter(pl.col("is_target") == "true")[
        "main_score"
    ].to_numpy()
    target_charges = scores_df.filter(pl.col("is_target") == "true")["precursor_charge"]
    decoy_scores = scores_df.filter(pl.col("is_target") == "false")[
        "main_score"
    ].to_numpy()
    decoy_charges = scores_df.filter(pl.col("is_target") == "false")["precursor_charge"]

    pprint("Number of targets: {}".format(len(target_scores)))
    pprint("Number of decoys: {}".format(len(decoy_scores)))
    if (len(target_scores) + len(decoy_scores)) != len(scores_df):
        raise ValueError("Error filtering targets and decoys")

    scores = np.log1p(np.concatenate((target_scores, decoy_scores)))
    bins = np.histogram_bin_edges(scores[~np.isnan(scores)], bins=50)

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True)

    ax[0].hist(scores, bins, alpha=0.5, label="All Scores")
    ax[0].hist(np.log1p(target_scores), bins, alpha=0.5, label="Target Scores")
    ax[0].hist(np.log1p(decoy_scores), bins, alpha=0.5, label="Decoy Scores")
    ax[0].set_xlabel("Main Score (log1p)")
    ax[0].set_ylabel("Count")
    ax[0].legend(loc="upper right")
    ax[0].set_yscale("log")

    ax[1].hist(scores, bins, alpha=0.5, label="All Scores")
    ax[1].hist(np.log1p(target_scores), bins, alpha=0.5, label="Target Scores")
    ax[1].hist(np.log1p(decoy_scores), bins, alpha=0.5, label="Decoy Scores")
    ax[1].set_xlabel("Main Score (log1p)")
    ax[1].set_ylabel("Count")
    ax[1].legend(loc="upper right")

    target_file = output_dir / "plot_mainscores.png"
    pprint(f"Saving plot to {target_file}")
    plt.title("Histogram of 'main_score' scores.")
    plt.savefig(target_file)
    plt.close()

    ## Re-make the same plot but filtering for charge states
    uniq_charges = np.unique(np.concatenate((target_charges, decoy_charges)))
    pprint(uniq_charges)

    fig, ax = plt.subplots(
        nrows=len(uniq_charges), ncols=1, figsize=(10, 8), sharex=True
    )

    for i, charge in enumerate(uniq_charges):
        target_scores_loc = target_scores[target_charges == charge]
        decoy_scores_loc = decoy_scores[decoy_charges == charge]
        loc_scores = np.log1p(np.concatenate((target_scores_loc, decoy_scores_loc)))

        ax[i].hist(loc_scores, bins, alpha=0.5, label="All Scores")
        ax[i].hist(np.log1p(target_scores_loc), bins, alpha=0.5, label="Target Scores")
        ax[i].hist(np.log1p(decoy_scores_loc), bins, alpha=0.5, label="Decoy Scores")
        ax[i].set_xlabel(f"Main Score (log1p); charge={charge}")
        ax[i].set_ylabel("Count")
        ax[i].legend(loc="upper right")
        ax[i].set_yscale("log")

    target_file = output_dir / "plot_mainscores_by_charge.png"
    pprint(f"Saving plot to {target_file}")
    # plt.title("Histogram of 'main_score' scores.")
    plt.savefig(target_file)
    plt.close()


def main_score_hist(df: pl.LazyFrame, output_dir: Path):
    pprint("Plotting main scores")
    scores_df = df.select(["is_target", "main_score", "precursor_charge"]).collect()
    target_scores = scores_df.filter(pl.col("is_target") == "true")[
        "main_score"
    ].to_numpy()
    target_charges = scores_df.filter(pl.col("is_target") == "true")["precursor_charge"]
    decoy_scores = scores_df.filter(pl.col("is_target") == "false")[
        "main_score"
    ].to_numpy()
    decoy_charges = scores_df.filter(pl.col("is_target") == "false")["precursor_charge"]

    pprint("Number of targets: {}".format(len(target_scores)))
    pprint("Number of decoys: {}".format(len(decoy_scores)))
    if (len(target_scores) + len(decoy_scores)) != len(scores_df):
        raise ValueError("Error filtering targets and decoys")

    scores = np.log1p(np.concatenate((target_scores, decoy_scores)))
    bins = np.histogram_bin_edges(scores[~np.isnan(scores)], bins=50)

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True)

    ax[0].hist(scores, bins, alpha=0.5, label="All Scores")
    ax[0].hist(np.log1p(target_scores), bins, alpha=0.5, label="Target Scores")
    ax[0].hist(np.log1p(decoy_scores), bins, alpha=0.5, label="Decoy Scores")
    ax[0].set_xlabel("Main Score (log1p)")
    ax[0].set_ylabel("Count")
    ax[0].legend(loc="upper right")
    ax[0].set_yscale("log")

    ax[1].hist(scores, bins, alpha=0.5, label="All Scores")
    ax[1].hist(np.log1p(target_scores), bins, alpha=0.5, label="Target Scores")
    ax[1].hist(np.log1p(decoy_scores), bins, alpha=0.5, label="Decoy Scores")
    ax[1].set_xlabel("Main Score (log1p)")
    ax[1].set_ylabel("Count")
    ax[1].legend(loc="upper right")

    target_file = output_dir / "plot_mainscores.png"
    pprint(f"Saving plot to {target_file}")
    plt.title("Histogram of 'main_score' scores.")
    plt.savefig(target_file)
    plt.close()

    ## Re-make the same plot but filtering for charge states
    uniq_charges = np.unique(np.concatenate((target_charges, decoy_charges)))
    pprint(uniq_charges)

    fig, ax = plt.subplots(
        nrows=len(uniq_charges), ncols=1, figsize=(10, 8), sharex=True
    )

    for i, charge in enumerate(uniq_charges):
        target_scores_loc = target_scores[target_charges == charge]
        decoy_scores_loc = decoy_scores[decoy_charges == charge]
        loc_scores = np.log1p(np.concatenate((target_scores_loc, decoy_scores_loc)))

        ax[i].hist(loc_scores, bins, alpha=0.5, label="All Scores")
        ax[i].hist(np.log1p(target_scores_loc), bins, alpha=0.5, label="Target Scores")
        ax[i].hist(np.log1p(decoy_scores_loc), bins, alpha=0.5, label="Decoy Scores")
        ax[i].set_xlabel(f"Main Score (log1p); charge={charge}")
        ax[i].set_ylabel("Count")
        ax[i].legend(loc="upper right")
        ax[i].set_yscale("log")

    target_file = output_dir / "plot_mainscores_by_charge.png"
    pprint(f"Saving plot to {target_file}")
    # plt.title("Histogram of 'main_score' scores.")
    plt.savefig(target_file)
    plt.close()


def plot_scores_hist(df: pl.LazyFrame, columns: list[str], output_dir: Path):
    targets_df = df.filter(pl.col("is_target") == "true")
    decoys_df = df.filter(pl.col("is_target") == "false")

    pprint(f"Plotting scores for columns: {columns}")
    for score in columns:
        pprint(f"Plotting histogram for score: {score}")
        target_scores = targets_df[score].to_numpy()
        decoy_scores = decoys_df[score].to_numpy()

        scores = np.concatenate((target_scores, decoy_scores))
        bins = np.histogram_bin_edges(scores[~np.isnan(scores)], bins=50)

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True)

        ax[0].hist(scores, bins, alpha=0.5, label="All Scores")
        ax[0].hist(target_scores, bins, alpha=0.5, label="Target Scores")
        ax[0].hist(decoy_scores, bins, alpha=0.5, label="Decoy Scores")
        ax[0].set_xlabel(f"Score: {score}")
        ax[0].set_ylabel("Count")
        ax[0].legend(loc="upper right")
        ax[0].set_yscale("log")

        ax[1].hist(scores, bins, alpha=0.5, label="All Scores")
        ax[1].hist(target_scores, bins, alpha=0.5, label="Target Scores")
        ax[1].hist(decoy_scores, bins, alpha=0.5, label="Decoy Scores")
        ax[1].set_xlabel(f"Score: {score}")
        ax[1].set_ylabel("Count")
        ax[1].legend(loc="upper right")

        target_file = output_dir / f"score_plot_{score}.png"
        pprint(f"Saving plot to {target_file}")
        plt.title(f"Histogram of '{score}' scores.")
        plt.savefig(target_file)
        plt.close()
