import numpy as np


def supersimpleprediction(mz, charge):
    intercept_ = -1.660e00
    log1p_mz = np.log1p(mz)
    sq_mz_over_charge = (mz**2) / charge
    log1p_sq_mz_over_charge = np.log1p(sq_mz_over_charge)

    out = (
        intercept_
        + (-3.798e-01 * log1p_mz)
        + (-2.389e-04 * mz)
        + (3.957e-01 * log1p_sq_mz_over_charge)
        + (4.157e-07 * sq_mz_over_charge)
        + (1.417e-01 * charge)
    )
    return out
