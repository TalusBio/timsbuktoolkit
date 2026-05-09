import polars as pl

from bench.entrapment import (
    PeptideClass,
    analyse,
    classify_peptides,
    compute_fdr_curve,
    parse_fasta,
    plot_fdr_curve,
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


def test_compute_fdr_curve_basic():
    classified = pl.DataFrame(
        {
            "qvalue": [0.001, 0.005, 0.01, 0.02, 0.05],
            "class": ["target", "target", "entrapment", "target", "entrapment"],
        }
    )
    curve = compute_fdr_curve(classified)
    # Sorted ascending by qvalue
    assert curve["qvalue"].to_list() == [0.001, 0.005, 0.01, 0.02, 0.05]
    # n_target cumulative
    assert curve["n_target"].to_list() == [1, 2, 2, 3, 3]
    # n_entrap cumulative
    assert curve["n_entrap"].to_list() == [0, 0, 1, 1, 2]
    # empirical_fdr = n_e / (n_t + n_e)
    last = curve.row(-1, named=True)
    assert last["empirical_fdr"] == 2 / 5


def test_compute_fdr_curve_excludes_shared_and_unknown():
    classified = pl.DataFrame(
        {
            "qvalue": [0.01, 0.01, 0.01, 0.01],
            "class": ["target", "shared_dropped", "unknown", "entrapment"],
        }
    )
    curve = compute_fdr_curve(classified)
    # Only one target + one entrapment row contribute
    assert curve.height == 2
    assert sorted(curve["class"].to_list()) == ["entrapment", "target"]


def test_plot_fdr_curve_writes_png(tmp_path):
    curve = pl.DataFrame(
        {
            "qvalue": [0.001, 0.01, 0.05],
            "n_target": [10, 50, 100],
            "n_entrap": [0, 1, 5],
            "empirical_fdr": [0.0, 1 / 51, 5 / 105],
        }
    )
    out = tmp_path / "fdr.png"
    plot_fdr_curve(curve, out, title="test")
    assert out.exists()
    assert out.stat().st_size > 1000  # not an empty or stub file


def test_analyse_end_to_end(tmp_path):
    target = tmp_path / "t.fasta"
    target.write_text(">T1\nAAAAPEPTIDEKBBBB\n")
    entrap = tmp_path / "e.fasta"
    entrap.write_text(">E1\nQQQQENTRAPEPTKZZZZ\n")

    results = pl.DataFrame(
        {
            "sequence": ["PEPTIDEK", "ENTRAPEPTK", "PEPTIDEK", "ENTRAPEPTK"],
            "qvalue": [0.001, 0.02, 0.005, 0.04],
        }
    )
    results_path = tmp_path / "results.parquet"
    results.write_parquet(results_path)

    out = analyse(
        results_parquet=results_path,
        target_fasta=target,
        entrapment_fasta=entrap,
        out_parquet=tmp_path / "classified.parquet",
        out_plot=tmp_path / "fdr.png",
    )

    # Returned scalars
    assert out["entrap/n_target_at_q01"] == 2  # both PEPTIDEK rows have q <= 0.01
    assert out["entrap/n_entrap_at_q01"] == 0
    assert out["entrap/empirical_fdr_at_q01"] == 0.0
    assert out["entrap/n_target_at_q05"] == 2
    assert out["entrap/n_entrap_at_q05"] == 2

    # Outputs
    assert (tmp_path / "classified.parquet").exists()
    assert (tmp_path / "fdr.png").exists()

    classified = pl.read_parquet(tmp_path / "classified.parquet")
    assert "class" in classified.columns and "is_entrapment" in classified.columns
