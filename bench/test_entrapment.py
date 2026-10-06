import csv

import pytest

from bench.entrapment import (
    Hit,
    estimate,
    generate,
    load_pairs,
    run_search,
    stripped_sequence,
)


def test_generate_writes_distinct_matched_pairs(tmp_path):
    parent = tmp_path / "parent.fasta"
    parent.write_text(">one\nACDEFGHKLMNPQSTVWYK\n")
    result = generate(parent, tmp_path / "experiment", seed=7)

    assert result["pairs"] == 2
    labels, pairing = load_pairs(result["pair_file"])
    assert len(labels) == 4
    assert len(pairing) == 2
    for entrapment, target in pairing.items():
        assert entrapment != target
        assert entrapment[0] == target[0]
        assert entrapment[-1] == target[-1]
        assert sorted(entrapment) == sorted(target)
    with result["pair_file"].open(newline="") as stream:
        assert len(list(csv.DictReader(stream, delimiter="\t"))) == 4


def test_paired_fdp_tracks_an_entrapment_before_its_target():
    labels = {
        "AAAAAAK": "target",
        "AAAADAK": "p_target",
        "CCCCCCK": "target",
        "CCCCDCK": "p_target",
    }
    pairing = {"AAAADAK": "AAAAAAK", "CCCCDCK": "CCCCCCK"}
    hits = [
        Hit("AAAAAAK", "AAAAAAK", 2, 0.01, 10),
        Hit("AAAADAK", "AAAADAK", 2, 0.02, 9),
        Hit("CCCCDCK", "CCCCDCK", 2, 0.03, 8),
        Hit("CCCCCCK", "CCCCCCK", 2, 0.04, 7),
    ]
    report = estimate(hits, labels, pairing)

    assert [row["paired_fdp"] for row in report] == [0, 0.5, 1, 1]
    assert report[2]["lower_bound_fdp"] == pytest.approx(2 / 3)
    assert report[2]["combined_fdp"] == pytest.approx(4 / 3)


def test_equal_q_values_share_the_end_of_group_estimate():
    labels = {"AAAAAAK": "target", "AAAADAK": "p_target"}
    pairing = {"AAAADAK": "AAAAAAK"}
    hits = [
        Hit("AAAAAAK", "AAAAAAK", 2, 0.01, 9),
        Hit("AAAADAK", "AAAADAK", 2, 0.01, 10),
    ]

    assert [row["paired_fdp"] for row in estimate(hits, labels, pairing)] == [1.5, 1.5]


def test_stripped_sequence_accepts_timsseek_proforma():
    assert stripped_sequence("_PEPC[UniMod:4]TIDEK_") == "PEPCTLDEK"
    with pytest.raises(ValueError):
        stripped_sequence("not a peptide")


def test_run_preserves_generated_length_range_and_full_q_values(tmp_path, monkeypatch):
    fasta = tmp_path / "entrapment.fasta"
    fasta.write_text(f">short\n{'A' * 7}\n>long\n{'A' * 34}\n")
    calls = []
    monkeypatch.setattr(
        "bench.entrapment.subprocess.run", lambda cmd, check: calls.append(cmd)
    )

    run_search(
        tmp_path / "timsseek",
        fasta,
        [tmp_path / "one.d", tmp_path / "two.d"],
        tmp_path / "search",
        tmp_path / "entrapment.mzspeclib.txt.gz",
        None,
    )

    build, search = calls
    assert build[build.index("--min-length") + 1] == "7"
    assert build[build.index("--max-length") + 1] == "34"
    assert search[search.index("--max-qvalue") + 1] == "1"
    assert search.count("--raw-inputs") == 2
