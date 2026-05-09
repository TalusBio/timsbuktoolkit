import polars as pl

from bench.entrapment import (
    PeptideClass,
    analyse,
    classify_peptides,
    compute_fdr_curve,
    count_kmers,
    kmer_normalization_factor,
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

    df = pl.DataFrame({"sequence": ["PEPTIDEK", "ENTRAPEPTK", "SHARED", "GHOSTAA"]})
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
    classified = pl.DataFrame({
        "qvalue": [0.001, 0.005, 0.01, 0.02, 0.05],
        "class": ["target", "target", "entrapment", "target", "entrapment"],
    })
    curve = compute_fdr_curve(classified)
    # Sorted ascending by qvalue
    assert curve["qvalue"].to_list() == [0.001, 0.005, 0.01, 0.02, 0.05]
    # n_target cumulative
    assert curve["n_target"].to_list() == [1, 2, 2, 3, 3]
    # n_entrap cumulative
    assert curve["n_entrap"].to_list() == [0, 0, 1, 1, 2]
    # empirical_fdr = n_e / (n_t + n_e)
    last = curve.row(-1, named=True)
    assert last["empirical_fdr_raw"] == 2 / 5
    assert last["empirical_fdr_norm"] == 2 / 5  # factor defaults to 1.0


def test_compute_fdr_curve_excludes_shared_and_unknown():
    classified = pl.DataFrame({
        "qvalue": [0.01, 0.01, 0.01, 0.01],
        "class": ["target", "shared_dropped", "unknown", "entrapment"],
    })
    curve = compute_fdr_curve(classified)
    # Only one target + one entrapment row contribute
    assert curve.height == 2
    assert sorted(curve["class"].to_list()) == ["entrapment", "target"]


def test_plot_fdr_curve_writes_png(tmp_path):
    curve = pl.DataFrame({
        "qvalue": [0.001, 0.01, 0.05],
        "n_target": [10, 50, 100],
        "n_entrap": [0, 1, 5],
        "empirical_fdr_raw": [0.0, 1 / 51, 5 / 105],
        "empirical_fdr_norm": [0.0, 1 / 51, 5 / 105],
    })
    out = tmp_path / "fdr.png"
    plot_fdr_curve(curve, out, title="test")
    assert out.exists()
    assert out.stat().st_size > 1000  # not an empty or stub file


def test_analyse_end_to_end(tmp_path):
    target = tmp_path / "t.fasta"
    target.write_text(">T1\nAAAAPEPTIDEKBBBB\n")
    entrap = tmp_path / "e.fasta"
    entrap.write_text(">E1\nQQQQENTRAPEPTKZZZZ\n")

    results = pl.DataFrame({
        "sequence": ["PEPTIDEK", "ENTRAPEPTK", "PEPTIDEK", "ENTRAPEPTK"],
        "qvalue": [0.001, 0.02, 0.005, 0.04],
    })
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
    assert out["entrap/empirical_fdr_raw_at_q01"] == 0.0
    assert out["entrap/n_target_at_q05"] == 2
    assert out["entrap/n_entrap_at_q05"] == 2

    # Outputs
    assert (tmp_path / "classified.parquet").exists()
    assert (tmp_path / "fdr.png").exists()

    classified = pl.read_parquet(tmp_path / "classified.parquet")
    assert "class" in classified.columns and "is_entrapment" in classified.columns


def test_count_kmers_basic(tmp_path):
    p = tmp_path / "p.fasta"
    p.write_text(">A\nABCDEFG\n>B\nABCDEFGH\n")
    kmers = count_kmers(p, k=7)
    # First protein contributes "ABCDEFG" exactly (length 7 → one kmer).
    # Second protein contributes "ABCDEFG" and "BCDEFGH".
    # Set dedupes the shared kmer.
    assert kmers == {"ABCDEFG", "BCDEFGH"}


def test_count_kmers_skips_too_short(tmp_path):
    p = tmp_path / "p.fasta"
    p.write_text(">A\nABCDEF\n")  # length 6, smaller than k=7
    assert count_kmers(p, k=7) == set()


def test_kmer_normalization_factor(tmp_path):
    target = tmp_path / "t.fasta"
    # 4 kmers of length 7: ABCDEFG, BCDEFGH, CDEFGHI, DEFGHIJ
    target.write_text(">T\nABCDEFGHIJ\n")
    entrap = tmp_path / "e.fasta"
    entrap.write_text(">E\nABCDEFG\n")  # 1 kmer: ABCDEFG (shared with target)

    # After dropping shared kmers: target has {BCDEFGH, CDEFGHI, DEFGHIJ} (3);
    # entrap has {} (0). Factor = target / max(1, entrap) → 3/1 = 3.0.
    f = kmer_normalization_factor(target, entrap, k=7)
    assert f == 3.0


def test_kmer_normalization_factor_balanced(tmp_path):
    target = tmp_path / "t.fasta"
    target.write_text(">T\nAAAAAAAA\n")  # 2 kmers: AAAAAAA, AAAAAAA → set: {AAAAAAA}
    entrap = tmp_path / "e.fasta"
    entrap.write_text(">E\nBBBBBBBB\n")  # 2 kmers: BBBBBBB, BBBBBBB → set: {BBBBBBB}
    # No shared kmers; both have 1 → factor = 1.0
    f = kmer_normalization_factor(target, entrap, k=7)
    assert f == 1.0


def test_compute_fdr_curve_with_normalization():
    """Plain raw FDR uses n_e / (n_e + n_t); normalized scales n_e by factor."""
    classified = pl.DataFrame(
        {
            "qvalue": [0.001, 0.005, 0.01, 0.02, 0.05],
            "class": ["target", "target", "entrapment", "target", "entrapment"],
        }
    )
    curve = compute_fdr_curve(classified, normalization_factor=3.0)

    # Raw counts unchanged
    assert curve["n_target"].to_list() == [1, 2, 2, 3, 3]
    assert curve["n_entrap"].to_list() == [0, 0, 1, 1, 2]
    # Raw fdr unchanged
    assert curve["empirical_fdr_raw"].to_list() == [0.0, 0.0, 1 / 3, 1 / 4, 2 / 5]
    # Normalized fdr at last row: n_e * 3 / (n_t + n_e * 3) = 6 / (3 + 6) = 2/3
    last = curve.row(-1, named=True)
    assert last["empirical_fdr_norm"] == 2 / 3
    assert last["n_entrap_norm"] == 6.0  # 2 * 3.0


def test_compute_fdr_curve_default_factor_is_one():
    """Default factor=1.0 keeps backward-compatible raw == norm behavior."""
    classified = pl.DataFrame(
        {
            "qvalue": [0.001, 0.01],
            "class": ["target", "entrapment"],
        }
    )
    curve = compute_fdr_curve(classified)
    assert curve["empirical_fdr_norm"].to_list() == curve["empirical_fdr_raw"].to_list()


def test_analyse_includes_normalization_factor(tmp_path):
    """analyse() computes and applies kmer normalization automatically."""
    target = tmp_path / "t.fasta"
    # homopolymer of length 50 → exactly 1 unique kmer "AAAAAAA"
    target.write_text(">T1\n" + "A" * 50 + "\n")
    entrap = tmp_path / "e.fasta"
    entrap.write_text(">E1\n" + "B" * 50 + "\n")

    results = pl.DataFrame(
        {"sequence": ["AAAAAAA", "BBBBBBB"], "qvalue": [0.001, 0.001]}
    )
    results_path = tmp_path / "r.parquet"
    results.write_parquet(results_path)
    out = analyse(
        results_parquet=results_path,
        target_fasta=target,
        entrapment_fasta=entrap,
        out_parquet=tmp_path / "c.parquet",
        out_plot=tmp_path / "f.png",
    )
    # Both proteins are pure homopolymer of length 50 → 1 unique kmer each → factor=1.0
    # Returned scalars include the factor
    assert "entrap/normalization_factor" in out
    assert out["entrap/normalization_factor"] == 1.0
