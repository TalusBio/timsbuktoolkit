import polars as pl

from bench.entrapment import (
    PeptideClass,
    classify_peptides,
    parse_fasta,
    strip_mods,
)


def test_strip_mods():
    assert strip_mods("PEPTIDEK") == "PEPTIDEK"
    assert strip_mods("PEPC[U:4]TIDEK") == "PEPCTIDEK"
    assert strip_mods("PEP(Phospho)TIDEK") == "PEPTIDEK"
    assert strip_mods("123.45PEPTIDEK") == "PEPTIDEK"
    assert strip_mods("n[42]PEPTIDEK") == "nPEPTIDEK"  # keep alpha n-term marker


def test_parse_fasta(tmp_path):
    p = tmp_path / "p.fasta"
    p.write_text(">sp|P1|A\nMKLAA\nDDDD\n>sp|P2|B\nLLLL\n")
    out = parse_fasta(p)
    assert out == {"sp|P1|A": "MKLAADDDD", "sp|P2|B": "LLLL"}


def test_classify_peptides(tmp_path):
    target = tmp_path / "t.fasta"
    target.write_text(">T1\nAAAAPEPTIDEKBBBB\n>T2\nMMMMSHAREDXXXX\n")
    entrap = tmp_path / "e.fasta"
    entrap.write_text(">E1\nQQQQENTRAPEPTKZZZZ\n>E2\nMMMMSHAREDYYYY\n")

    df = pl.DataFrame(
        {"sequence": ["PEPTIDEK", "ENTRAPEPTK", "SHARED", "GHOSTAA"]}
    )
    classified = classify_peptides(df, target, entrap)

    classes = dict(zip(classified["sequence"], classified["class"]))
    assert classes["PEPTIDEK"] == PeptideClass.TARGET.value
    assert classes["ENTRAPEPTK"] == PeptideClass.ENTRAPMENT.value
    assert classes["SHARED"] == PeptideClass.SHARED_DROPPED.value
    assert classes["GHOSTAA"] == PeptideClass.UNKNOWN.value
    # is_entrapment column: True only for ENTRAPMENT
    is_e = dict(zip(classified["sequence"], classified["is_entrapment"]))
    assert is_e["ENTRAPEPTK"] is True
    assert is_e["PEPTIDEK"] is False
    assert is_e["SHARED"] is False
    assert is_e["GHOSTAA"] is False


def test_classify_peptides_strips_mods_before_match(tmp_path):
    target = tmp_path / "t.fasta"
    target.write_text(">T1\nAAAAPEPTIDEKBBBB\n")
    entrap = tmp_path / "e.fasta"
    entrap.write_text(">E1\nQQQQ\n")

    df = pl.DataFrame({"sequence": ["PEPC[U:4]PTIDEK"]})
    # Stripped form is PEPCPTIDEK which is NOT in target. Confirm we end up unknown.
    classified = classify_peptides(df, target, entrap)
    assert classified["class"][0] == PeptideClass.UNKNOWN.value

    # But a real target match works after stripping
    df2 = pl.DataFrame({"sequence": ["PEPT[U:4]IDEK"]})
    classified2 = classify_peptides(df2, target, entrap)
    assert classified2["class"][0] == PeptideClass.TARGET.value
