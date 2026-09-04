#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["polars", "numpy", "scikit-learn"]
# ///
"""Refit `supersimpleprediction` 1/k0 model from a timsseek results.parquet.

Mirrors the feature set in
rust/timsseek/src/fragment_mass/elution_group_converter.rs:supersimpleprediction.
Prints the new intercept + coefs ready to paste back into the Rust source.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.linear_model import HuberRegressor, LinearRegression


def build_features(mz: np.ndarray, z: np.ndarray) -> np.ndarray:
    sq_mz_over_z = mz**2 / z
    return np.column_stack([
        np.log1p(mz),
        mz,
        np.log1p(sq_mz_over_z),
        sq_mz_over_z,
        z,
    ])


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("parquet", type=Path, help="path to results.parquet")
    p.add_argument(
        "--min-main-score",
        type=float,
        default=0.0,
        help="filter rows with main_score < threshold (default: 0.0)",
    )
    p.add_argument(
        "--robust",
        action="store_true",
        help="use HuberRegressor instead of OLS",
    )
    p.add_argument(
        "--holdout",
        type=float,
        default=0.2,
        help="fraction of rows reserved for holdout MAPE (default: 0.2)",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
    )
    args = p.parse_args()

    if not args.parquet.exists():
        print(f"parquet not found: {args.parquet}", file=sys.stderr)
        return 1

    df = pl.read_parquet(args.parquet).filter(
        pl.col("obs_mobility").is_finite()
        & pl.col("precursor_mz").is_finite()
        & (pl.col("main_score") > args.min_main_score)
    )
    print(f"rows after filter: {df.height}")
    if df.height < 100:
        print("not enough rows to fit", file=sys.stderr)
        return 1

    mz = df["precursor_mz"].to_numpy().astype(np.float64)
    z = df["precursor_charge"].cast(pl.Float64).to_numpy()
    y = df["obs_mobility"].to_numpy().astype(np.float64)

    X = build_features(mz, z)

    rng = np.random.default_rng(args.seed)
    idx = rng.permutation(len(y))
    cut = int(len(y) * (1 - args.holdout))
    tr, ho = idx[:cut], idx[cut:]

    model = HuberRegressor(max_iter=500) if args.robust else LinearRegression()
    model.fit(X[tr], y[tr])

    def mape(a: np.ndarray, b: np.ndarray) -> float:
        return float(np.mean(np.abs((a - b) / b)) * 100)

    tr_mape = mape(model.predict(X[tr]), y[tr])
    ho_mape = mape(model.predict(X[ho]), y[ho])

    feats = [
        "log1p_mz",
        "mz",
        "log1p_sq_mz_over_charge",
        "sq_mz_over_charge",
        "charge",
    ]
    print()
    print(f"intercept: {model.intercept_:.6e}")
    for name, coef in zip(feats, model.coef_):
        print(f"  {name:>24s}: {coef:+.6e}")
    print()
    print(f"train MAPE: {tr_mape:.4f}%")
    print(f"holdout MAPE: {ho_mape:.4f}%")

    print()
    print("--- Rust paste block ---")
    print(f"    let intercept_ = {model.intercept_:.3e};")
    print("    let log1p_mz = (mz + 1.).ln();")
    print("    let sq_mz_over_charge = mz.powi(2) / charge as f64;")
    print("    let log1p_sq_mz_over_charge = (sq_mz_over_charge + 1.).ln();")
    print()
    print("    intercept_")
    c = model.coef_
    print(f"        + ({c[0]:+.3e} * log1p_mz)")
    print(f"        + ({c[1]:+.3e} * mz)")
    print(f"        + ({c[2]:+.3e} * log1p_sq_mz_over_charge)")
    print(f"        + ({c[3]:+.3e} * sq_mz_over_charge)")
    print(f"        + ({c[4]:+.3e} * charge as f64)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
