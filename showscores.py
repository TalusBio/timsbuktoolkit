# /// script
# dependencies = [
#   "polars",
#   "rich",
#   "matplotlib",
#   "numpy",
#   "tqdm",
#   "mokapot @ git+https://github.com/jspaezp/mokapot.git@feat/re_add_compound_index_spec",
#   "xgboost",
#   "torch",
#   "uniplot",
# ]
# ///

import argparse
import logging
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from turtle import forward
from typing import List, Optional

import matplotlib.pyplot as plt
import mokapot
import numpy as np
import pandas as pd
import polars as pl
import torch
import torch.nn as nn
import torch.nn.functional as F
import xgboost as xgb
from mokapot.column_defs import ColumnGroups, OptionalColumns
from rich.pretty import pprint
from torch.utils.data import DataLoader, TensorDataset
from tqdm.auto import tqdm
from uniplot import histogram

# from sklearn.model_selection import GridSearchCV
# from xgboost import XGBClassifier

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results_dir",
        type=str,
        nargs="+",
        help="Path to the directories containing the results.csv files",
    )
    parser.add_argument(
        "-o", "--output_dir", type=str, default=".", help="Output directory"
    )
    return parser


def read_files(results_dirs: list[Path]) -> pl.LazyFrame:
    files = set()
    for results_dir in results_dirs:
        files.update(results_dir.glob("*.parquet"))

    files = list(files)
    pprint(f"Scanning {len(files)} files")
    data = pl.scan_parquet(files)
    pprint("Done scanning")
    return data


def lazy_abs_and_maxfill(df: pl.LazyFrame, columns: list[str]) -> pl.LazyFrame:
    exprs_first = []
    exprs_later = []
    for column in columns:
        exprs_first.append(pl.col(column).abs())
        exprs_later.append(pl.col(column).fill_nan(pl.col(column).max()))

    return df.with_columns(exprs_first).with_columns(exprs_later)


def lazy_zero_fill(df: pl.LazyFrame, columns: list[str]) -> pl.LazyFrame:
    exprs = []
    for column in columns:
        exprs.append(pl.col(column).fill_nan(0))
    return df.with_columns(exprs)


def log_cols(df: pl.LazyFrame, columns: list[str]) -> pl.LazyFrame:
    exprs = []
    for col in columns:
        exprs.append(pl.col(col).log1p().fill_nan(0))
    return df.with_columns(exprs)


def add_id(df: pl.LazyFrame) -> pl.LazyFrame:
    # Add an id column ... pretty simple incremental integer
    df = df.with_row_index(name="id")
    return df


def check_noninf(df: pl.DataFrame, columns: list[str]):
    pprint("Checking for infinite values")
    any_inf = False
    df_inf = df.select(columns).filter(pl.any_horizontal(pl.all().is_infinite()))
    pprint(df_inf)
    for col in columns:
        if df_inf[col].is_infinite().any():
            pprint(f"Column {col} has infinite values")
            any_inf = True
    if any_inf:
        raise ValueError("Data contains infinite values")


def check_nonnan(df: pl.LazyFrame, columns: list[str]):
    pprint("Checking for NaN values")
    nan_cols = []
    df_nan = df.select(columns).filter(pl.any_horizontal(pl.all().is_nan()))
    pprint(df_nan)
    for col in columns:
        if df_nan[col].is_nan().any():
            pprint(f"Column {col} has NaN values")
            nan_cols.append(col)
    if nan_cols:
        raise ValueError(f"Data contains NaN values: {nan_cols}")


def check_nonexp(df: pl.DataFrame, columns: list[str]):
    # checks that things are not exponential
    # The heuristic here is that stuff that is log transformed
    # should only span 3 orders of magnitude

    norders = {}
    failing_cols = []
    df = df.select(columns)
    for col in columns:
        try:
            mag_diff = df[col].abs().max() - df[col].abs().min()
            nmags = math.log10(mag_diff)
        except ValueError as e:
            if "math domain error" in str(e):
                raise ValueError(f"Column {col} has 0 variance values")
        norders[col] = nmags
        if nmags > 3:
            failing_cols.append(col)
    pprint(norders)
    if failing_cols:
        ranges = [
            f"{col}: min={df[col].min()}, max={df[col].max()} norders={norders[col]}"
            for col in failing_cols
        ]
        pprint(ranges)
        raise ValueError(f"Data contains exponential values: {failing_cols}")


def cast_f32(df: pl.LazyFrame, cols: list[str]) -> pl.LazyFrame:
    exprs = []
    for col in cols:
        exprs.append(pl.col(col).cast(pl.Float32))
    return df.with_columns(exprs)


def ohe_charges(df: pl.LazyFrame, charges: list[int]) -> tuple[pl.LazyFrame, list[str]]:
    exprs = []
    colnames = []
    for charge in charges:
        colname = f"charge_{charge}"
        colnames.append(colname)
        exprs.append((pl.col("precursor_charge") == charge).alias(colname))
    return df.with_columns(exprs), tuple(colnames)


def scale_columns(
    df: pl.LazyFrame, cols: list[tuple[str, float]]
) -> tuple[pl.LazyFrame, tuple[str, ...]]:
    exprs = []
    cols_out = []
    for col, factor in cols:
        exprs.append(pl.col(col).cast(pl.Float32) / factor)
        cols_out.append(col)
    return df.with_columns(exprs), tuple(cols_out)


def td_compete(df_use: pl.DataFrame) -> pl.DataFrame:
    pprint("Stripping sequences")
    stripped_seqs = (
        df_use["sequence"].str.replace_all("\/\d+", "").str.replace_all("\[.*?\]", "")
    )
    mods = [
        tuple(x)
        for x in df_use["sequence"]
        .str.replace_all("\/\d+", "")
        .str.extract_all("\[.*?\]")
        .list.sort()
        .to_list()
    ]
    df_use = df_use.with_columns(
        td_id=pl.Series(
            derive_td_pair(
                stripped_seqs.to_list(),
                mods,
                charges=df_use["precursor_charge"].to_list(),
            )
        )
    )
    df_use = df_use.with_columns(
        delta_td_score=pl.when(pl.col("main_score").count() > 1)
        .then(
            pl.col("main_score")
            - pl.col("main_score").min().over(["td_id", "precursor_charge"])
        )
        .otherwise(pl.col("main_score")),
    )

    pprint("Initial T/D competition")
    init_nrow = len(df_use)
    df_use = df_use.filter(
        pl.col("main_score")
        == pl.col("main_score").max().over(["td_id", "precursor_charge"])
    )
    pprint(
        f"T/D competition Removed {init_nrow - len(df_use)} rows (kept {len(df_use)})"
    )

    return df_use


def to_mokapot_df(df: pl.LazyFrame) -> tuple[pd.DataFrame, ColumnGroups]:
    pprint("Starting to_mokapot")
    loggable_cols = (
        # Log
        "npeaks",
        "lazyerscore",
        "main_score",
        "delta_next",
        "delta_second_next",
        "lazyerscore_vs_baseline",
        "norm_lazyerscore_vs_baseline",
        "ms1_summed_precursor_intensity",
        "ms2_summed_transition_intensity",
        # TODO: consider clamping instead of logging here.
        "sq_delta_theo_rt",
    )
    imputable_cols = (
        # Abs impute
        "ms2_mz_error_0",
        "ms2_mz_error_1",
        "ms2_mz_error_2",
        "ms2_mz_error_3",
        "ms2_mz_error_4",
        "ms2_mz_error_5",
        "ms2_mz_error_6",
        "ms2_mobility_error_0",
        "ms2_mobility_error_1",
        "ms2_mobility_error_2",
        "ms2_mobility_error_3",
        "ms2_mobility_error_4",
        "ms2_mobility_error_5",
        "ms2_mobility_error_6",
        "ms1_mz_error_0",
        "ms1_mz_error_1",
        "ms1_mz_error_2",
        "ms1_mobility_error_0",
        "ms1_mobility_error_1",
        "ms1_mobility_error_2",
        "sq_delta_ms1_ms2_mobility",
        "delta_ms1_ms2_mobility",
    )
    scaling_cols = (
        ("precursor_rt_query_seconds", 60),
        ("obs_rt_seconds", 60),
        ("delta_theo_rt", 60),
    )
    # zero_imputable_cols = ("ms1_ms2_correlation",)
    zero_imputable_cols = ()
    generated_cols = [
        "delta_td_score",
    ]

    df_use = log_cols(df, loggable_cols)
    df_use = lazy_abs_and_maxfill(df_use, imputable_cols)
    df_use = lazy_zero_fill(df_use, zero_imputable_cols)
    df_use, scaling_cols = scale_columns(df_use, scaling_cols)
    pprint("Collecting")
    df_use, ohe_cols = ohe_charges(df_use, charges=[2, 3, 4])
    generated_cols += ohe_cols
    df_use = add_id(df_use).collect(streaming=True)
    df_use = td_compete(df_use)

    feat_cols = (
        (
            "precursor_charge",
            "precursor_mz",
            "precursor_mobility_query",
            "obs_mobility",
            "ms2_cosine_ref_similarity",
            "ms2_coelution_score",
            "ms1_cosine_ref_similarity",
            "ms1_coelution_score",
            "nqueries",
            "ms1_inten_ratio_2",
            "ms2_inten_ratio_4",
            "ms2_inten_ratio_6",
            "ms1_inten_ratio_1",
            "ms2_inten_ratio_2",
            "ms2_inten_ratio_1",
            "ms2_inten_ratio_3",
            "ms1_inten_ratio_0",
            "ms2_inten_ratio_5",
            "ms2_inten_ratio_0",
        )
        + loggable_cols
        + imputable_cols
        + zero_imputable_cols
        + tuple(generated_cols)
        + scaling_cols
    )

    # This requires all columns to exist, so we do it after all preprocessing
    df_use = cast_f32(df_use, feat_cols)
    check_noninf(df_use, feat_cols)
    check_nonnan(df_use, feat_cols)
    check_nonexp(df_use, feat_cols)
    pprint("Converting to pandas")
    df_use = df_use.filter(pl.col("main_score") > 1)
    df_use = df_use.to_pandas()
    cols = ColumnGroups(
        columns=df_use.columns,
        target_column="is_target",
        peptide_column="sequence",
        # spectrum_columns=("id", "td_id"),
        spectrum_columns=("td_id",),
        feature_columns=feat_cols,
        extra_confidence_level_columns=(),
        optional_columns=OptionalColumns(
            id=None,
            filename=None,
            scan=None,
            calcmass="precursor_mz",
            # charge="precursor_charge",
            charge=None,
            expmass=None,
            rt="obs_rt_seconds",
            protein=None,
        ),
    )
    nonfeat_cols = set(cols.columns) - set(cols.feature_columns)
    pprint(f"Non-feature columns: {nonfeat_cols}")
    return df_use, cols


def derive_td_pair(
    sequences: list[str], mods: list[tuple[str, ...]], charges: list[int]
) -> list[int]:
    pprint("Deriving TD pairs")
    td_pairs = {}
    max_id = 0

    out = []
    for seq, mod, charge in tqdm(
        zip(sequences, mods, charges, strict=True), total=len(sequences)
    ):
        if (seq, mod, charge) in td_pairs:
            out.append(td_pairs[(seq, mod, charge)])
            continue

        dec_seq = seq[0] + seq[1:-1][::-1] + seq[-1]
        td_pairs[(seq, mod, charge)] = td_pairs[(dec_seq, mod, charge)] = max_id
        out.append(max_id)
        max_id += 1

    pprint(f"Found {max_id} TD pairs in {len(sequences)} sequences/charge pairs")
    return out


def to_mokapot(df: pl.LazyFrame) -> mokapot.LinearPsmDataset:
    df_use, cols = to_mokapot_df(df)
    return mokapot.LinearPsmDataset(
        df_use,
        column_groups=cols,
    )


def to_folds_xgb(
    *,
    shuffled_df: pd.DataFrame,
    cols: ColumnGroups,
    num_folds: int,
) -> list[xgb.DMatrix]:
    tmp = to_folds(shuffled_df=shuffled_df, cols=cols, num_folds=num_folds)
    out = [
        xgb.DMatrix(
            x[0],
            label=x[1],
            feature_names=cols.feature_columns,
        )
        for x in tmp
    ]
    return out


def to_folds(
    *,
    shuffled_df: pd.DataFrame,
    cols: ColumnGroups,
    num_folds: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    out = np.array_split(shuffled_df, num_folds)

    out = [
        (
            x.loc[:, cols.feature_columns].to_numpy(),
            np.where(x.loc[:, cols.target_column].to_numpy(), 1, 0),
        )
        for x in out
    ]
    return out


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


def xgboost_stuff(df: pl.LazyFrame, output_dir: Path):
    df_use, cols = to_mokapot_df(df)
    pprint("Shuffling")
    df_use = df_use.sample(frac=1).reset_index(drop=True, inplace=False)
    folds = to_folds_xgb(shuffled_df=df_use, num_folds=5, cols=cols)
    fold_model = KFoldModel.from_folds(folds)
    fold_model.train()
    fold_model.score()
    importances = fold_model.get_importances()
    pprint(importances)
    fig = plot_importances(importances)
    outfile = output_dir / "importances.png"
    fig.savefig(outfile)
    pprint(f"Wrote {outfile}")
    plt.close()

    ctargs = fold_model.concat_targets()
    cscores = fold_model.concat_scores()
    df_use["rescore_score"] = cscores
    df_use["qvalue"] = mokapot.qvalues.qvalues_from_scores(cscores, ctargs == 1)
    df_use = df_use.sort_values("rescore_score", ascending=False)
    outfile = output_dir / "rescored_values.parquet"
    df_use.to_parquet(outfile, index=False)
    pprint(f"Wrote {outfile}")
    order = np.argsort(-cscores)
    ctargs = ctargs[order]
    cscores = cscores[order]
    cumtargs = np.cumsum(ctargs)

    target_preds = cscores[ctargs == 1]
    decoy_preds = cscores[ctargs == 0]

    histogram(target_preds, title="Target scores")
    histogram(decoy_preds, title="Decoy scores")

    qvals = mokapot.qvalues.qvalues_from_scores(cscores, ctargs == 1)
    for ct in [0.01, 0.05, 0.1, 0.5, 1.0]:
        num_at_thresh = int(np.sum(ctargs[qvals < ct]))
        ssc = cscores[qvals < ct]
        if len(ssc) == 0:
            pprint(f"No scores at {ct}")
            continue
        score_at_thresh = np.min(ssc)
        pprint(f"Score at {ct}: {score_at_thresh}")
        pprint(f"Number of targets at {ct}: {num_at_thresh}")

    report = Report(
        targets_at_1=np.sum(ctargs[qvals < 0.01]).item(),
        targets_at_5=np.sum(ctargs[qvals < 0.05]).item(),
        targets_at_10=np.sum(ctargs[qvals < 0.1]).item(),
    )
    pprint(report)
    outfile = output_dir / "report.toml"
    report.save_to_toml(outfile)
    pprint(f"Wrote {outfile}")

    # plt.plot(qvals, cumtargs, label="All Scores")

    mask = qvals < 0.1
    plt.plot(qvals[mask], cumtargs[mask], label="QValues < 0.1")
    plt.legend(loc="upper right")
    plt.xlim(0, 0.1)
    plt.title("Cumulative number of accepted peptides.")
    outfile = output_dir / "plot_qvalues.png"
    plt.savefig(outfile)
    pprint(f"Wrote {outfile}")
    plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(target_preds, alpha=0.5, bins=100, label="Targets")
    ax[0].hist(decoy_preds, alpha=0.5, bins=100, label="Decoys")
    # ax[0].axvline(x=score_at_onepct, color="k", linestyle="--", alpha=0.5)
    ax[0].legend()

    ax[1].hist(target_preds, alpha=0.5, bins=100, label="Targets")
    ax[1].hist(decoy_preds, alpha=0.5, bins=100, label="Decoys")
    # ax[1].axvline(x=score_at_onepct, color="k", linestyle="--", alpha=0.5)
    ax[1].set_yscale("log")
    ax[1].legend()
    plt.title("Histogram of 'rescoring' score.")
    outfile = output_dir / "plot_hist.png"
    plt.savefig(outfile)
    pprint(f"Wrote {outfile}")
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


def mokapot_stuff(data, outdir):
    ds = to_mokapot(data)
    pprint("Brewing Mokapot")
    models, scores = mokapot.brew([ds])
    qvals = mokapot.qvalues.qvalues_from_scores(scores[0], ds.targets)
    for ct in [0.01, 0.05, 0.1, 0.5, 1.0]:
        num_at_thresh = np.sum(ds.targets[qvals < ct])
        ssc = scores[0][qvals < ct]
        if len(ssc) == 0:
            pprint(f"No scores at {ct}")
            continue
        score_at_thresh = np.min(ssc)
        pprint(
            f"Mokapot Number of targets at {ct}: {num_at_thresh}; Score: {score_at_thresh}"
        )


@dataclass
class Report:
    targets_at_1: int
    targets_at_5: int
    targets_at_10: int

    def save_to_toml(self, path: Path):
        with open(path, "w") as f:
            f.write("[report]\n")
            f.write(f"targets_at_1 = {self.targets_at_1}\n")
            f.write(f"targets_at_5 = {self.targets_at_5}\n")
            f.write(f"targets_at_10 = {self.targets_at_10}\n")


@dataclass
class KFoldModel:
    """
    The idea here is to have N folds that will be used in this way:
    1. Train on 1 fold.
    2. Use 1 fold for early stopping.
    3. For prediction, use all other folds (models where the fold was not used for training or early stopping).
    """

    folds: list[xgb.DMatrix]
    models: list[xgb.Booster | None]
    scores: list[np.ndarray | None]

    @staticmethod
    def from_folds(folds: list[xgb.DMatrix]):
        return KFoldModel(folds, [None] * len(folds), [None] * len(folds))

    def train(self):
        for i in range(len(self.folds)):
            pprint(f"Training model {i}/{len(self.folds)}")
            train = self.folds[i]
            early_stop = self.folds[(i + 1) % len(self.folds)]
            model = xgb.train(
                {"objective": "binary:logistic", "scale_pos_weight": 0.2},
                train,
                num_boost_round=200,
                evals=[
                    (early_stop, "validation"),
                ],
                early_stopping_rounds=5,
                verbose_eval=1,
            )
            self.models[i] = model

    def score(self):
        for i in range(len(self.folds)):
            pprint(f"Scoring fold {i}/{len(self.folds)}")
            scores = []
            for j, model in enumerate(self.models):
                if j == i:
                    continue
                if j == (i + 1) % len(self.folds):
                    continue
                scores.append(model.predict(self.folds[i]))
            self.scores[i] = np.mean(scores, axis=0)

    def get_importances(self):
        imps = [model.get_score(importance_type="gain") for model in self.models]
        imps_order = sorted(imps[0].items(), key=lambda x: x[1], reverse=True)
        out = {k[0]: [w.get(k[0], 0) for w in imps] for k in imps_order}
        return out

    def concat_scores(self):
        if self.scores[0] is None:
            raise ValueError("Scores not computed")
        return np.concatenate(self.scores)

    def concat_targets(self):
        targets = [x.get_label() for x in self.folds]
        return np.concatenate(targets)


######## Neural network based rescoring ########


class AsymmetricMarginBCELoss(nn.Module):
    def __init__(self, *, margin_0=0.1, margin_1=0.4):
        """
        margin_0: margin for negative class (0s)
        margin_1: margin for positive class (1s)

        This loss pushes negative predictions further from decision boundary

        Larger margin_0 makes the model more conservative about predicting 1s
            (reduces false positives)
        Smaller margin_1 means we're more lenient about false negatives
        Both margins are clamped to keep predictions in [0,1] range
        """
        super().__init__()
        self.margin_0 = margin_0
        self.margin_1 = margin_1

    def forward(self, predictions, targets):
        # Add margins to predictions based on true class
        adjusted_preds = torch.where(
            targets == 1, predictions + self.margin_1, predictions - self.margin_0
        )

        # Clamp to valid probability range
        adjusted_preds = torch.clamp(adjusted_preds, 0, 1)

        # Compute BCE loss
        loss = F.binary_cross_entropy(adjusted_preds, targets, reduction="mean")

        return loss


class WeightedBCELoss(nn.Module):
    def __init__(self, pos_weight=0.2, neg_weight=1.0):
        """
        pos_weight: weight for positive class (1s)
        neg_weight: weight for negative class (0s)
        """
        super().__init__()
        self.pos_weight = pos_weight
        self.neg_weight = neg_weight

    def forward(self, predictions, targets):
        # Create weight tensor based on targets
        weights = torch.where(
            targets == 1,
            torch.tensor(self.pos_weight, device=targets.device),
            torch.tensor(self.neg_weight, device=targets.device),
        )

        # Standard BCE loss
        bce_loss = F.binary_cross_entropy(predictions, targets, reduction="none")

        # Apply weights
        weighted_loss = weights * bce_loss

        return weighted_loss.mean()


class FocalLoss3(nn.Module):
    def __init__(self, alpha=0.25, gamma=2.0, reduce=True):
        """
        Focal Loss: (1 - p)^gamma * log(p) for positive class
                   p^gamma * log(1-p) for negative class

        alpha: weighing factor for positive class
        gamma: focusing parameter that reduces the loss contribution from easy examples
        """
        super().__init__()
        if alpha < 0 or alpha > 1:
            raise ValueError("Alpha must be in [0, 1]")
        self.alpha = alpha
        self.gamma = gamma
        self.reduce = reduce

    def forward(self, predictions, targets):
        # BCE loss
        bce_loss = F.binary_cross_entropy(predictions, targets, reduction="none")

        # Focal term
        pt = torch.where(targets == 1, predictions, 1 - predictions)
        focal_term = (1 - pt) ** self.gamma

        # Alpha weighing
        alpha_weight = torch.where(
            targets == 1,
            torch.tensor(self.alpha, device=targets.device),
            torch.tensor(1 - self.alpha, device=targets.device),
        )

        loss = alpha_weight * focal_term * bce_loss
        if self.reduce:
            return loss.mean()
        else:
            return loss


class BinaryClassifier(nn.Module):
    def __init__(
        self,
        input_dim: int,
        nhidden_layers: int = 4,
        hidden_dims: int = 64,
        dropout: float = 0.1,
    ):
        super().__init__()

        # ACTIVATION = nn.ReLU
        # Selu seems to be critical to train deeper networks.
        ACTIVATION = nn.SELU

        layers = []
        layers.append(nn.BatchNorm1d(input_dim))
        layers.append(nn.Linear(input_dim, hidden_dims))
        layers.append(ACTIVATION())
        self.input_layer = nn.Sequential(*layers)

        layers = []
        for _ in range(nhidden_layers):
            layers.append(nn.Linear(hidden_dims, hidden_dims))
            layers.append(ACTIVATION())
            if dropout > 0.0:
                layers.append(nn.Dropout(dropout))

        self.hidden_layers = nn.ModuleList(layers)

        layers = []
        layers.append(nn.Linear(hidden_dims, 1))
        layers.append(nn.Sigmoid())
        self.output_layer = nn.Sequential(*layers)

        pprint(self)

    def forward(self, x):
        x = self.input_layer(x)
        for layer in self.hidden_layers:
            # x = x + layer(layer_norm(x))
            x = layer(x)

        return self.output_layer(x)


@dataclass
class MLPKFoldModel:
    """
    PyTorch implementation of K-Fold cross validation.
    Each fold contains:
    1. One fold for training
    2. One fold for validation (early stopping)
    3. Remaining folds for inference
    """

    folds: List[tuple[torch.Tensor, torch.Tensor]]  # List of (features, targets) tuples
    models: List[Optional[BinaryClassifier]]
    scores: List[Optional[torch.Tensor]]
    device: torch.device

    @staticmethod
    def from_folds(
        folds: List[tuple[torch.Tensor, torch.Tensor]], device: torch.device
    ):
        return MLPKFoldModel(
            folds=folds,
            models=[None] * len(folds),
            scores=[None] * len(folds),
            device=device,
        )

    def train(
        self,
        batch_size: int = 124,
        epochs: int = 20,
        learning_rate: float = 1e-4,
        pos_weight: float = 0.2,
        **kwargs,
    ):
        for i in range(len(self.folds)):
            print(f"Training model {i}/{len(self.folds)}")

            # Prepare data
            train_data = self.folds[i]
            val_data = self.folds[(i + 1) % len(self.folds)]

            train_dataset = TensorDataset(train_data[0], train_data[1])
            val_dataset = TensorDataset(val_data[0], val_data[1])

            train_loader = DataLoader(
                train_dataset, batch_size=batch_size, shuffle=True
            )
            val_loader = DataLoader(val_dataset, batch_size=batch_size)

            # Initialize model
            model = BinaryClassifier(input_dim=train_data[0].shape[1], **kwargs).to(
                self.device
            )
            optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate)
            # Choose one of the following loss functions based on your needs:
            # criterion = WeightedBCELoss(pos_weight=pos_weight, neg_weight=1.0)
            # criterion = AsymmetricMarginBCELoss(margin_0=0.5, margin_1=0.2)
            # criterion = WeightedBCELoss2(fneg_weight=pos_weight)

            # Usually gamma is positive bc the desire is to emphasize well classified
            # Examples, whilst we actually want to focus in misclassified examples
            # Where they are kind of "hard" to distinguish
            criterion = FocalLoss3(alpha=pos_weight, gamma=0.5)

            # Training loop
            best_val_loss = float("inf")
            patience = 5
            patience_counter = 0
            best_model = None

            for epoch in range(epochs):
                # Training
                model.train()
                train_losses = []
                for x_batch, y_batch in train_loader:
                    x_batch, y_batch = x_batch.to(self.device), y_batch.to(self.device)

                    optimizer.zero_grad()
                    y_pred = model(x_batch)
                    loss = criterion(y_pred, y_batch.view(-1, 1))
                    loss.backward()
                    optimizer.step()

                    train_losses.append(loss.item())

                # Validation
                model.eval()
                val_losses = []
                with torch.no_grad():
                    for x_batch, y_batch in val_loader:
                        x_batch, y_batch = (
                            x_batch.to(self.device),
                            y_batch.to(self.device),
                        )
                        y_pred = model(x_batch)
                        val_loss = criterion(y_pred, y_batch.view(-1, 1))
                        val_losses.append(val_loss.item())

                avg_val_loss = np.mean(val_losses)

                # Early stopping
                print(
                    f"Epoch {epoch}: train_loss = {np.mean(train_losses):.4f}, val_loss = {avg_val_loss:.4f}"
                )
                if avg_val_loss < best_val_loss:
                    best_val_loss = avg_val_loss
                    best_model = model.state_dict()
                    patience_counter = 0
                else:
                    patience_counter += 1
                    if patience_counter >= patience:
                        print(
                            f"Early stopping at epoch {epoch}. Best validation loss: {best_val_loss:.4f}"
                        )
                        break

            # Save best model
            model.load_state_dict(best_model)
            self.models[i] = model

    def score(self, batch_size: int = 32):
        for i in range(len(self.folds)):
            print(f"Scoring fold {i}/{len(self.folds)}")
            fold_data = self.folds[i]
            dataset = TensorDataset(fold_data[0], fold_data[1])
            loader = DataLoader(dataset, batch_size=batch_size)

            scores_list = []

            for j, model in enumerate(self.models):
                if j == i or j == (i + 1) % len(self.folds):
                    continue

                model.eval()
                fold_scores = []

                with torch.no_grad():
                    for x_batch, _ in loader:
                        x_batch = x_batch.to(self.device)
                        predictions = model(x_batch)
                        fold_scores.append(predictions.cpu())

                scores_list.append(torch.cat(fold_scores))

            self.scores[i] = torch.stack(scores_list).mean(dim=0)

    def get_importances(self, feat_names: list[str] | None = None):
        # Note: Feature importance is not as straightforward in neural networks
        # This is a simple implementation using gradient-based importance
        importances = []

        for model in self.models:
            importance_dict = {}

            # Compute average gradient magnitude for each feature
            for i, (features, targets) in enumerate(self.folds):
                features.requires_grad_(True)
                output = model(features.to(self.device))
                output.sum().backward()

                grad_magnitude = features.grad.abs().mean(dim=0)
                importance_dict = {}
                for j in range(features.shape[1]):
                    feat_name = (
                        feat_names[j] if feat_names is not None else f"feature_{j}"
                    )
                    importance_dict[feat_name] = grad_magnitude[j].item()

                features.requires_grad_(False)
                model.zero_grad()

            importances.append(importance_dict)

        # Format similar to original
        imps_order = sorted(importances[0].items(), key=lambda x: x[1], reverse=True)
        out = {k[0]: [w.get(k[0], 0) for w in importances] for k in imps_order}
        return out

    def get_importances2(self, feat_names: list[str] | None = None):
        """
        Calculate feature importance using gradients after BatchNorm layer.
        Uses hooks to capture intermediate gradients.
        """
        importances = []

        for model in self.models:
            importance_dict = {}
            post_bn_gradients = []

            # Register hook to capture gradients after BatchNorm
            def hook_fn(module, grad_input, grad_output):
                # post_bn_gradients.append(grad_input[0].detach().cpu())
                post_bn_gradients.append(grad_output[0].detach().cpu())

            # Get the BatchNorm layer
            bn_layer = model.input_layer[0]  # First layer is BatchNorm
            hook = bn_layer.register_backward_hook(hook_fn)

            # Compute gradients for each fold
            for i, (features, targets) in enumerate(self.folds):
                features = features.to(self.device)
                features.requires_grad_(True)
                post_bn_gradients = []  # Reset for each batch

                # Forward and backward pass
                output = model(features)
                output.sum().backward()

                # Average gradients across samples in the fold
                fold_grads = post_bn_gradients[0].abs().mean(dim=0)

                # Update importance dictionary
                if not importance_dict:
                    if feat_names is None:
                        feat_names = [f"feature_{j}" for j in range(features.shape[1])]
                    importance_dict = {
                        k: fold_grads[j].item() for j, k in enumerate(feat_names)
                    }
                else:
                    for j, k in enumerate(feat_names):
                        importance_dict[k] += fold_grads[j].item()

                features.requires_grad_(False)
                model.zero_grad()

            # Average importance across folds
            for key in importance_dict:
                importance_dict[key] /= len(self.folds)

            importances.append(importance_dict)

            # Remove the hook
            hook.remove()

        # Format similar to original
        imps_order = sorted(importances[0].items(), key=lambda x: x[1], reverse=True)
        out = {k[0]: [w.get(k[0], 0) for w in importances] for k in imps_order}
        return out

    def concat_scores(self):
        if self.scores[0] is None:
            raise ValueError("Scores not computed")
        return torch.cat(self.scores)

    def concat_targets(self):
        return torch.cat([fold[1] for fold in self.folds])


def to_torch_folds(shuffled_df: pl.LazyFrame, num_folds: int, cols):
    tmp = to_folds(shuffled_df=shuffled_df, cols=cols, num_folds=num_folds)
    out = [
        (torch.from_numpy(x[0]).float(), torch.from_numpy(x[1]).float()) for x in tmp
    ]
    return out


def mlp_stuff(
    df_use: pl.LazyFrame, cols, output_dir: Path, pos_weight: float = 0.05, **kwargs
):
    pprint("Shuffling")
    df_use = df_use.sample(frac=1).reset_index(drop=True, inplace=False)
    folds = to_torch_folds(shuffled_df=df_use, num_folds=5, cols=cols)
    fold_model = MLPKFoldModel.from_folds(folds, device="cpu")
    fold_model.train(pos_weight=pos_weight, **kwargs)
    fold_model.score()
    # importances = fold_model.get_importances(feat_names=cols.feature_columns)
    importances = fold_model.get_importances2(feat_names=cols.feature_columns)
    pprint(importances)
    fig = plot_importances(importances)
    outfile = output_dir / "importances_nn.png"
    fig.savefig(outfile)
    pprint(f"Wrote {outfile}")
    plt.close()

    ctargs = fold_model.concat_targets().numpy().flatten()
    cscores = fold_model.concat_scores().numpy().flatten()
    df_use["rescore_score"] = cscores
    df_use["qvalue"] = mokapot.qvalues.qvalues_from_scores(cscores, ctargs == 1)
    df_use = df_use.sort_values("rescore_score", ascending=False)
    outfile = output_dir / "rescored_values_nn.parquet"
    df_use.to_parquet(outfile, index=False)
    pprint(f"Wrote {outfile}")
    order = np.argsort(-cscores)
    ctargs = ctargs[order]
    cscores = cscores[order]
    qvals = mokapot.qvalues.qvalues_from_scores(cscores, ctargs == 1)
    for ct in [0.01, 0.05, 0.1, 0.5, 1.0]:
        num_at_thresh = int(np.sum(ctargs[qvals < ct]))
        ssc = cscores[qvals < ct]
        if len(ssc) == 0:
            pprint(f"No scores at {ct}")
            continue
        score_at_thresh = np.min(ssc)
        pprint(f"Score at {ct}: {score_at_thresh}")
        pprint(f"Number of targets at {ct}: {num_at_thresh}")

    target_preds = cscores[ctargs == 1]
    decoy_preds = cscores[ctargs == 0]
    histogram(target_preds, title="Target scores")
    histogram(decoy_preds, title="Decoy scores")

    report = Report(
        targets_at_1=np.sum(ctargs[qvals < 0.01]).item(),
        targets_at_5=np.sum(ctargs[qvals < 0.05]).item(),
        targets_at_10=np.sum(ctargs[qvals < 0.1]).item(),
    )
    pprint(report)
    outfile = output_dir / "report_nn.toml"
    report.save_to_toml(outfile)
    pprint(f"Wrote {outfile}")
    return np.sum(ctargs[qvals < 0.01]).item()


####### End neural network based rescoring #######


def main(args):
    paths = [Path(p) for p in args.results_dir]
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(f"Path {p} does not exist")
    data = read_files(paths)
    data = data.filter(pl.col("obs_mobility").is_not_nan())

    outdir = Path(args.output_dir)
    if not outdir.exists():
        outdir.mkdir(parents=True)

    main_score_hist(data, outdir)
    xgboost_stuff(data, outdir)
    mokapot_stuff(data, outdir)
    data, cols = to_mokapot_df(data)
    score = mlp_stuff(data, cols=cols, output_dir=outdir)

    exit()
    # Leaving the code I used for hparam tuning of the MLP
    # (0.2, 4, 512, 0)
    # (0.05, 3, 64, 0.1)
    scores = {}
    best_score = score
    best_params = None
    combs = []
    for pw in [0.05, 0.2, 0.5, 0.8]:
        # for do in [0, 0.05, 0.1, 0.2, 0.5]:
        for do in [0, 0.05, 0.1]:
            for hl in [3, 4, 5]:
                for hd in [64, 128]:
                    combs.append((hl, hd, pw, do))
            for hl in [1, 2, 3, 4, 5]:
                for hd in [256, 512]:
                    combs.append((hl, hd, pw, do))

    # shuffle combs
    import random

    random.shuffle(combs)
    for hl, hd, pw, do in combs:
        tmp_outdir = outdir / f"hl-{hl}_hd-{hd}-pw-{pw}-do-{do}"
        tmp_outdir.mkdir(parents=True, exist_ok=True)
        score = mlp_stuff(
            data,
            cols=cols,
            output_dir=tmp_outdir,
            pos_weight=pw,
            nhidden_layers=hl,
            hidden_dims=hd,
            dropout=do,
        )
        scores[(hl, hd, do)] = score
        if score > best_score:
            pprint(f"New best score: {score} at {pw} {hl} {hd} {do}")
            best_score = score
            best_params = (pw, hl, hd, do)
        else:
            pprint("Keeping old score")
            pprint(f"Current score {score}")
            pprint(
                f"Old best score: {best_score} at {best_params if best_params is not None else 'BASELINE'}"
            )

    sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    pprint(sorted_scores)
    pprint(best_params)


if __name__ == "__main__":
    parser = build_parser()
    args, unkargs = parser.parse_known_args()
    if unkargs:
        raise ValueError(f"Unknown arguments: {unkargs}")

    main(args)
