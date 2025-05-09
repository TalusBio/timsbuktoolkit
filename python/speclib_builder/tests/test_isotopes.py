from speclib_builder.isotopes import peptide_isotopes


def test_peptide_isotopes():
    """
    Test function to verify peptide isotope calculations.
    """
    iso = peptide_isotopes(60, 5)
    expected = [0.3972, 0.2824, 0.1869, 0.0846]
    expected = [val / 0.3972 for val in expected[:3]]  # Normalize first 3 values

    # Check if all differences are within tolerance
    tolerance = 0.02
    matched = all(abs(a - b) <= tolerance for a, b in zip(iso, expected, strict=True))

    assert matched, f"Test failed: {iso} != {expected}"
    print("Test passed successfully!")
