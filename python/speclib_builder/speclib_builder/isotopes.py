# Credit where credit is due ...
# This is essentially a port of the original code from Sage
# and its immplementation of the isotope distribution calculation

import numpy as np
from rustyms import MolecularFormula


def convolve(a: list[float], b: list[float]) -> tuple[float, float, float, float]:
    """
    Performs a custom convolution operation on two arrays of length 4.

    Args:
        a: First array of 4 floating point numbers
        b: Second array of 4 floating point numbers

    Returns:
        List of 4 floating point numbers representing the convolution result
    """
    return (
        a[0] * b[0],
        a[0] * b[1] + a[1] * b[0],
        a[0] * b[2] + a[1] * b[1] + a[2] * b[0],
        a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0],
    )


def carbon_isotopes(count: int) -> list[float]:
    """
    Calculates carbon isotope distributions.

    Args:
        count: Number of carbon atoms

    Returns:
        List of 4 floating point numbers representing isotope distributions
    """
    lambda_val = float(count) * 0.011
    c13 = [0.0] * 4
    fact = [1, 1, 2, 6]

    for k in range(4):
        c13[k] = pow(lambda_val, k) * np.exp(-lambda_val) / float(fact[k])

    return c13


def sulfur_isotopes(count: int) -> tuple[float, float, float, float]:
    """
    Calculates sulfur isotope distributions.

    Args:
        count: Number of sulfur atoms

    Returns:
        List of 4 floating point numbers representing convolved isotope distributions
    """
    lambda33 = float(count) * 0.0076
    lambda35 = float(count) * 0.044
    s33 = [0.0] * 4
    s35 = [
        pow(lambda35, 0) * np.exp(-lambda35),
        0.0,
        pow(lambda35, 1) * np.exp(-lambda35),
        0.0,
    ]

    fact = [1, 1, 2, 6]
    for k in range(4):
        s33[k] = pow(lambda33, k) * np.exp(-lambda33) / float(fact[k])

    return convolve(s33, s35)


def peptide_isotopes(carbons: int, sulfurs: int) -> tuple[float, float, float]:
    """
    Calculates peptide isotope distributions based on number of carbon and sulfur atoms.

    Args:
        carbons: Number of carbon atoms
        sulfurs: Number of sulfur atoms

    Returns:
        List of 3 floating point numbers representing normalized isotope distributions
    """
    c = carbon_isotopes(carbons)
    s = sulfur_isotopes(sulfurs)
    result = convolve(c, s)
    max_val = max(result[:3]).item()  # Only consider first 3 values for normalization

    # Normalize first 3 values
    return [val.item() / max_val for val in result[:3]]


def peptide_formula_dist(formula: MolecularFormula) -> tuple[float, float, float]:
    c_count = 0
    s_count = 0

    for elem, _, count in formula.elements():
        str_elem = str(elem)
        if str_elem == "C":
            c_count = count
        elif str_elem == "S":
            s_count = count
        else:
            continue

    if c_count == 0:
        raise ValueError("No carbons found in formula")

    return peptide_isotopes(c_count, s_count)
