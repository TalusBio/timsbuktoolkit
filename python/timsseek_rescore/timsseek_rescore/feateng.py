import math
from pathlib import Path

import mokapot
import pandas as pd
import polars as pl
from mokapot.column_defs import ColumnGroups, OptionalColumns
from rich.pretty import pprint
from tqdm.auto import tqdm


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
        exprs.append(pl.col(col).fill_nan(0).log1p().fill_nan(0))
    out = df.with_columns(exprs)
    return out


def add_id(df: pl.LazyFrame) -> pl.LazyFrame:
    # Add an id column ... pretty simple incremental integer
    df = df.with_row_index(name="id")
    return df


def check_noninf(df: pl.DataFrame, columns: list[str]):
    pprint("Checking for infinite values")
    any_inf = False
    df_inf = df.select(columns).filter(pl.any_horizontal(pl.all().is_infinite()))
    if len(df_inf) > 0:
        pprint(df_inf)
        pprint(df_inf[0].to_dict(as_series=False))
    for col in columns:
        if df_inf[col].is_infinite().any():
            ninf = df_inf[col].is_infinite().sum()
            pprint(f"Column {col} has infinite ({ninf}/{len(df)}) values")
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

    return df_use


def to_mokapot_df(
    df: pl.LazyFrame,
    make_nonmissing: bool = True,
    make_monotonic: bool = True,
    scale_cols: bool = True,
) -> tuple[pd.DataFrame, ColumnGroups]:
    pprint("Starting to_mokapot")
    loggable_cols = (
        # Log
        "npeaks",
        "main_score",
        "delta_next",
        "delta_second_next",
        "apex_lazyerscore",
        "apex_lazyerscore_vs_baseline",
        "ms2_isotope_lazyerscore",
        "ms2_lazyerscore",
        "ms2_isotope_lazyerscore_ratio",
        "ms1_summed_precursor_intensity",
        "ms2_summed_transition_intensity",
        # TODO: consider clamping instead of logging here.
        "sq_delta_theo_rt",
        "calibrated_sq_delta_theo_rt",
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
    generated_cols = []

    if scale_cols:
        df_use = log_cols(df, loggable_cols)
    if make_monotonic:
        df_use = lazy_abs_and_maxfill(df_use, imputable_cols)
    if make_nonmissing:
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
            "ms1_corr_v_gauss",
            "ms2_corr_v_gauss",
            "nqueries",
            # Intensity ratios
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
            # Cycle counts
            "raising_cycles",
            "falling_cycles",
            "apex_norm_lazyerscore_vs_baseline",
            # ...
            "delta_group_ratio",
            "delta_group",
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
    if make_nonmissing:
        check_nonnan(df_use, feat_cols)
    if scale_cols:
        check_nonexp(df_use, feat_cols)

    pprint("Converting to pandas")
    df_use = df_use.to_pandas()
    cols = ColumnGroups(
        columns=df_use.columns,
        target_column="is_target",
        peptide_column="sequence",
        # spectrum_columns=("id", "td_id"),
        spectrum_columns=("id",),
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


def read_files(results_dirs: list[Path]) -> pl.LazyFrame:
    files = set()
    for results_dir in results_dirs:
        files.update(results_dir.glob("results.parquet"))

    files = list(files)
    pprint(f"Scanning {len(files)} files -> {files}")
    data = pl.scan_parquet(files)
    pprint("Done scanning")
    return data


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
