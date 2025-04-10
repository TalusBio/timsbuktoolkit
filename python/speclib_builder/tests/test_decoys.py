import pytest
from speclib_builder.decoys import (
    DecoyStrategy,
    as_decoy,
    as_decoy_proforma,
)

SAMPLES = [
    (
        "MYPEPTIDEK",
        {
            DecoyStrategy.REVERSE: "MEDITPEPYK",
            DecoyStrategy.MUTATE: "LSLDLSVEDL",
            DecoyStrategy.EDGE_MUTATE: "MSPEPTIDDK",
        },
    ),
    (
        "MYPEPTIDEK/2",
        {
            DecoyStrategy.REVERSE: "MEDITPEPYK/2",
            DecoyStrategy.MUTATE: "LSLDLSVEDL/2",
            DecoyStrategy.EDGE_MUTATE: "MSPEPTIDDK/2",
        },
    ),
    (
        "MYPEPT[U:21]IDEK",
        {
            DecoyStrategy.REVERSE: "MEDIT[U:21]PEPYK",
            DecoyStrategy.MUTATE: "LSLDLSVEDL",
            DecoyStrategy.EDGE_MUTATE: "MSPEPT[U:21]IDDK",
        },
    ),
]


@pytest.mark.parametrize("peptide, expected", SAMPLES)
def test_as_decoy(peptide, expected):
    for strat, expect_seq in expected.items():
        assert as_decoy(peptide, strat) == expect_seq, f"{peptide}, {strat}"
        assert as_decoy_proforma(peptide, strat) == expect_seq, f"{peptide}, {strat}"
