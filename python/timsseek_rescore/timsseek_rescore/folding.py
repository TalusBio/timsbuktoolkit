import numpy as np
import pandas as pd
import xgboost as xgb
from mokapot.column_defs import ColumnGroups


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
