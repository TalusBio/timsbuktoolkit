import copy

from speclib_builder.base import (
    ElutionGroup,
    PrecursorEntry,
    SpeclibElement,
)
from speclib_builder.decoys import DecoyStrategy


def _sample_element():
    entry = {
        "precursor": {"sequence": "SIIQSAQQDSIKK/2", "charge": 2, "decoy": 0},
        "elution_group": {
            "id": 4,
            "mobility": 1.0175,
            "rt_seconds": 0,
            "precursors": [
                (-1, 722.8980675419999),
                (1, 723.4024),
                (2, 723.906732458),
                (3, 724.411064916),
            ],
            "fragments": [
                ("y9^1", 1004.5382),
                ("y11^2", 623.34436),
                ("y10^1", 1132.5968),
                ("y8^1", 917.50616),
                ("y11^1", 1245.6808),
                ("y10^2", 566.80231),
                ("b3^1", 314.20798),
                ("y4^1", 475.32495),
                ("y7^1", 846.46906),
                ("y6^1", 718.41046),
                ("y5^1", 590.35193),
                ("b4^1", 442.26657),
            ],
            "precursor_intensities": [
                0.001,
                1.0,
                0.6709999999999999,
                0.22512049999999995,
            ],
            "fragment_intensities": {
                "y9^1": 1.0,
                "y11^2": 0.72318065,
                "y10^1": 0.50104403,
                "y8^1": 0.32450068,
                "y11^1": 0.29743919,
                "y10^2": 0.23268574,
                "b3^1": 0.20915699,
                "y4^1": 0.18289551,
                "y7^1": 0.17597583,
                "y6^1": 0.14646909,
                "y5^1": 0.11264818,
                "b4^1": 0.060880262,
            },
        },
    }

    elem = SpeclibElement(**entry)

    return elem


def test_peptide_swap_same():
    elem = _sample_element()
    elem_bkp = copy.deepcopy(elem)
    new_elem = elem.swap_peptide(
        proforma=elem.precursor.sequence,
        decoy=False,
        id=0,
    )
    assert elem.precursor == new_elem.precursor
    assert new_elem.elution_group.id == 0

    # Make sure no modification in-place has happened
    assert elem_bkp.elution_group.id == 4
    assert elem.elution_group.id == 4

    assert elem.elution_group.precursors == new_elem.elution_group.precursors
    assert set(x[0] for x in elem.elution_group.fragments) == set(
        x[0] for x in new_elem.elution_group.fragments
    )
    orig_dict = {x[0]: x[1] for x in elem.elution_group.fragments}
    new_dict = {x[0]: x[1] for x in new_elem.elution_group.fragments}
    for k, _ in elem.elution_group.fragments:
        orig = orig_dict[k]
        new = new_dict[k]
        diff = orig - new

        assert abs(diff) < 0.002, f"{k}: {diff}; original: {orig}, new: {new}"


def test_peptide_swap_decoy():
    elem = _sample_element()
    new_pep = elem.precursor.sequence.replace("I", "A")
    new_elem = elem.swap_peptide(
        proforma=new_pep,
        decoy=True,
        id=0,
    )
    assert new_elem.precursor.decoy
    assert new_elem.precursor.sequence == new_pep
    assert new_elem.precursor.sequence != elem.precursor.sequence

    assert set(x[0] for x in elem.elution_group.fragments) == set(
        x[0] for x in new_elem.elution_group.fragments
    )
    sames = {}
    diffs = {}

    orig_dict = {x[0]: x[1] for x in elem.elution_group.fragments}
    new_dict = {x[0]: x[1] for x in new_elem.elution_group.fragments}
    for i, x in enumerate(elem.elution_group.fragments):
        orig = orig_dict[x[0]]
        new = new_dict[x[0]]
        diff = orig - new

        if abs(diff) < 0.002:
            sames[x[0]] = diff
        else:
            diffs[x[0]] = diff

    assert len(sames) == 0
    V_TO_A_MASS = 42.04807294523107
    expect_one_diff = {
        "y9^1",
        "y11^2",
        "y10^1",
        "y8^1",
        "y4^1",
        "y7^1",
        "y6^1",
        "y5^1",
    }
    expect_two_diff = {
        "y11^1",
        "b3^1",
        "b4^1",
    }
    expect_half_diff = {
        "y10^2",
    }

    for k in expect_one_diff:
        assert k in diffs
        assert (abs(diffs[k]) - V_TO_A_MASS) < 0.002

    for k in expect_two_diff:
        assert k in diffs
        assert (abs(diffs[k]) - (2 * V_TO_A_MASS)) < 0.002

    for k in expect_half_diff:
        assert k in diffs
        assert (abs(diffs[k]) - (V_TO_A_MASS / 2)) < 0.002


def test_peptide_decoy_pseudorev():
    elem = _sample_element()
    new_elem = elem.generate_massshift_decoy(id=0, decoy_strategy=DecoyStrategy.REVERSE)
    assert new_elem.precursor.decoy
    assert new_elem.precursor.sequence != elem.precursor.sequence
    assert new_elem.precursor.sequence == "SKISDQQASQIIK/2"

    # http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet
    expect = [
        (
            1445.79587,
            "y13^1",
        ),
        (
            1358.76385,
            "y12^1",
        ),
        (
            1230.66888,
            "y11^1",
        ),
        (
            1117.58482,
            "y10^1",
        ),
        (
            1030.55279,
            "y9^1",
        ),
        (
            915.52585,
            "y8^1",
        ),
        (
            787.46727,
            "y7^1",
        ),
        (
            659.40869,
            "y6^1",
        ),
        (
            588.37158,
            "y5^1",
        ),
        (
            501.33955,
            "y4^1",
        ),
        (
            373.28098,
            "y3^1",
        ),
        (
            260.19691,
            "y2^1",
        ),
        (
            147.11285,
            "y1^1",
        ),
    ]

    for mz, k in expect:
        # Note that here I am using the original element to check
        # for the presence of the key. BUT the new one for the mass.
        if k in elem.elution_group.fragments:
            assert abs(new_elem.elution_group.fragments[k] - mz) < 0.01
