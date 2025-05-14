from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mokapot
import numpy as np
import polars as pl
import xgboost as xgb
from rich.pretty import pprint
from uniplot import histogram

from ..datamodels import Report
from ..feateng import to_mokapot_df
from ..folding import to_folds_xgb
from ..plotting import plot_importances


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
