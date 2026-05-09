import polars as pl

from bench.entrapment import (
    PeptideClass,
    analyse,
    classify_peptides,
    compute_fdr_curve,
    load_pairing,
    load_peptide_set,
    plot_fdr_curve,
    plot_score_histogram,
    strip_mods,
)


def test_strip_mods():
    assert strip_mods("PEPTIDEK") == "PEPTIDEK"
    assert strip_mods("PEPC[U:4]TIDEK") == "PEPCTIDEK"
    assert strip_mods("PEP(Phospho)TIDEK") == "PEPTIDEK"
    assert strip_mods("123.45PEPTIDEK") == "PEPTIDEK"
    assert strip_mods("n[42]PEPTIDEK") == "nPEPTIDEK"


def test_load_peptide_set(tmp_path):
    p = tmp_path / "peps.txt"
    p.write_text("PEPTIDEK\nLAGEPRVK\n\nMRSEQGLAR\n")
    out = load_peptide_set(p)
    assert out == {"PEPTIDEK", "LAGEPRVK", "MRSEQGLAR"}


def test_load_pairing(tmp_path):
    p = tmp_path / "pairs.tsv"
    p.write_text("target_peptide\tentrap_peptide\nPEPTIDEK\tEDPEKTIK\nLAGEPRVK\tGAEPLRVK\n")
    out = load_pairing(p)
    assert out == {"PEPTIDEK": "EDPEKTIK", "LAGEPRVK": "GAEPLRVK"}


def test_classify_peptides_set_membership(tmp_path):
    target = {"PEPTIDEK", "SHARED", "ANOTHER"}
    entrap = {"ENTRAPEPTK", "SHARED"}

    df = pl.DataFrame({"sequence": ["PEPTIDEK", "ENTRAPEPTK", "SHARED", "GHOSTAA"]})
    classified = classify_peptides(df, target, entrap)

    classes = dict(zip(classified["sequence"], classified["class"]))
    assert classes["PEPTIDEK"] == PeptideClass.TARGET.value
    assert classes["ENTRAPEPTK"] == PeptideClass.ENTRAPMENT.value
    assert classes["SHARED"] == PeptideClass.SHARED_DROPPED.value
    assert classes["GHOSTAA"] == PeptideClass.UNKNOWN.value

    is_e = dict(zip(classified["sequence"], classified["is_entrapment"]))
    assert is_e["ENTRAPEPTK"] is True
    assert is_e["PEPTIDEK"] is False
    assert is_e["SHARED"] is False
    assert is_e["GHOSTAA"] is False


def test_classify_strips_mods_before_match():
    target = {"PEPTIDEK"}
    entrap = {"AAAAAAAA"}
    df = pl.DataFrame({"sequence": ["PEPT[U:4]IDEK", "PEPC[U:4]PTIDEK"]})
    classified = classify_peptides(df, target, entrap)
    classes = classified["class"].to_list()
    # stripped → "PEPTIDEK" ∈ target
    assert classes[0] == PeptideClass.TARGET.value
    # stripped → "PEPCPTIDEK" not in either set
    assert classes[1] == PeptideClass.UNKNOWN.value


def test_classify_filters_to_targets():
    """is_target=False rows are excluded from FDR walk regardless of class column.

    With the new design, classify_peptides classifies ALL rows by sequence
    membership; the FDR walk in compute_fdr_curve filters is_target=True.
    """
    target = {"PEP1", "PEP2", "PEP3"}
    entrap = {"ENT1", "ENT2"}
    df = pl.DataFrame({
        "sequence": ["PEP1", "ENT1", "PEP2", "ENT2"],
        "qvalue":   [0.001, 0.005, 0.01, 0.02],
        "is_target":[True,  True,  False, True],
    })
    classified = classify_peptides(df, target, entrap)
    curve = compute_fdr_curve(classified, ratio=1.0)
    # Only the 3 is_target=True rows enter the curve
    assert curve.height == 3
    # is_target=True: PEP1 (target, q=.001), ENT1 (entrap, q=.005), ENT2 (entrap, q=.02)
    last = curve.row(-1, named=True)
    assert last["n_target"] == 1
    assert last["n_entrap"] == 2


def test_compute_fdr_curve_lower_combined(tmp_path):
    classified = pl.DataFrame({
        "qvalue":    [0.001, 0.005, 0.01, 0.02, 0.05],
        "class":     ["target", "target", "entrapment", "target", "entrapment"],
        "is_target": [True, True, True, True, True],
        "main_score":[100.0, 90.0, 80.0, 70.0, 60.0],
    })
    curve = compute_fdr_curve(classified, ratio=1.0)
    # n_target cumulative
    assert curve["n_target"].to_list() == [1, 2, 2, 3, 3]
    # n_entrap cumulative
    assert curve["n_entrap"].to_list() == [0, 0, 1, 1, 2]
    # Lower at last: 2/(3+2) = 0.4
    last = curve.row(-1, named=True)
    assert abs(last["empirical_fdr_lower"] - 2/5) < 1e-9
    # Combined at last: 2 * (1+1/1) / (3+2) = 0.8 (factor 2 for r=1)
    assert abs(last["empirical_fdr_combined"] - 4/5) < 1e-9


def test_compute_fdr_curve_combined_with_ratio_above_one():
    classified = pl.DataFrame({
        "qvalue":    [0.01, 0.02],
        "class":     ["target", "entrapment"],
        "is_target": [True, True],
        "main_score":[100.0, 90.0],
    })
    curve = compute_fdr_curve(classified, ratio=2.0)
    # Combined: n_e * (1 + 1/2) / (n_e + n_t) = 1 * 1.5 / 2 = 0.75
    last = curve.row(-1, named=True)
    assert abs(last["empirical_fdr_combined"] - 0.75) < 1e-9


def test_compute_fdr_curve_matched_with_pairing():
    """Matched estimator requires pairing dict and main_score column.

    Construct: 2 targets, both at q≤s; their paired entrap peptides also score
    well (one beats its target, one loses).
    """
    classified = pl.DataFrame({
        "sequence": ["PEP1", "ENT1", "PEP2", "ENT2"],
        "qvalue":   [0.001, 0.002, 0.003, 0.004],  # all under threshold
        "class":    ["target", "entrapment", "target", "entrapment"],
        "is_target":[True, True, True, True],
        "main_score":[100.0, 200.0, 150.0, 50.0],  # ENT1 beats PEP1; ENT2 loses to PEP2
    })
    pairing = {"PEP1": "ENT1", "PEP2": "ENT2"}
    curve = compute_fdr_curve(classified, ratio=1.0, pairing=pairing)
    # n_e = 2, n_t = 2 at the end
    # n_p_t_s (entrap > paired_target ≥ s): ENT1 wins over PEP1, both above s → 1
    # n_p_s_t (entrap ≥ s > paired_target): both targets are above s, so this is 0
    last = curve.row(-1, named=True)
    expected = (2 + 0 + 2*1) / (2 + 2)  # = 1.0
    assert abs(last["empirical_fdr_matched"] - expected) < 1e-9


def test_compute_fdr_curve_matched_target_below_threshold():
    """Entrap discovered, paired target NOT under the chosen threshold.

    Walk only includes rows that ARE discovered (q ≤ threshold). To exercise
    the n_p_s_t branch, place the entrap row in the curve and the paired target
    OUT (e.g., qvalue > threshold). compute_fdr_curve walks ALL rows; we read
    off the row at the entrap's qvalue.
    """
    classified = pl.DataFrame({
        "sequence":  ["PEP1", "ENT1"],
        "qvalue":    [0.5, 0.001],   # PEP1 has high qvalue, ENT1 low
        "class":     ["target", "entrapment"],
        "is_target": [True, True],
        "main_score":[10.0, 100.0],
    })
    pairing = {"PEP1": "ENT1"}
    curve = compute_fdr_curve(classified, ratio=1.0, pairing=pairing)
    # Sorted by qvalue: ENT1 (q=0.001), PEP1 (q=0.5)
    # At ENT1's row (first row): n_t=0, n_e=1
    # ENT1's paired target PEP1 has q=0.5 > 0.001 (not yet discovered)
    # → n_p_s_t = 1 (entrap discovered, paired target NOT yet)
    # → n_p_t_s = 0
    first = curve.row(0, named=True)
    assert first["n_target"] == 0
    assert first["n_entrap"] == 1
    expected = (1 + 1 + 0) / (0 + 1)
    assert abs(first["empirical_fdr_matched"] - expected) < 1e-9


def test_compute_fdr_curve_no_pairing_no_matched_column():
    classified = pl.DataFrame({
        "qvalue":    [0.01, 0.02],
        "class":     ["target", "entrapment"],
        "is_target": [True, True],
        "main_score":[100.0, 90.0],
    })
    curve = compute_fdr_curve(classified, ratio=1.0)
    assert "empirical_fdr_matched" not in curve.columns
    assert "empirical_fdr_lower" in curve.columns
    assert "empirical_fdr_combined" in curve.columns


def test_compute_fdr_curve_excludes_shared_and_unknown():
    classified = pl.DataFrame({
        "qvalue":    [0.01, 0.01, 0.01, 0.01],
        "class":     ["target", "shared_dropped", "unknown", "entrapment"],
        "is_target": [True,    True,             True,      True],
        "main_score":[100.0,   90.0,             80.0,      70.0],
    })
    curve = compute_fdr_curve(classified, ratio=1.0)
    assert curve.height == 2  # only target + entrap survive
    assert sorted(curve["class"].to_list()) == ["entrapment", "target"]


def test_plot_fdr_curve_writes_png(tmp_path):
    curve = pl.DataFrame({
        "qvalue":                  [0.001, 0.01, 0.05],
        "n_target":                [10, 50, 100],
        "n_entrap":                [0, 1, 5],
        "empirical_fdr_lower":     [0.0, 1/51, 5/105],
        "empirical_fdr_combined":  [0.0, 2/51, 10/105],
    })
    out = tmp_path / "fdr.png"
    plot_fdr_curve(curve, out, title="test")
    assert out.exists() and out.stat().st_size > 1000


def test_plot_score_histogram_writes_png(tmp_path):
    classified = pl.DataFrame({
        "main_score": [1e3, 1e4, 1e5, 1e6, 1e7] * 5,
        "class":      ["target","entrapment","target","entrapment","target"] * 5,
        "is_target":  [True, True, False, False, True] * 5,
    })
    out = tmp_path / "hist.png"
    plot_score_histogram(classified, out, title="test")
    assert out.exists() and out.stat().st_size > 1000


def test_analyse_end_to_end(tmp_path):
    target_peps = tmp_path / "t.txt"
    target_peps.write_text("PEP1\nPEP2\nPEP3\n")
    entrap_peps = tmp_path / "e.txt"
    entrap_peps.write_text("ENT1\nENT2\n")

    results = pl.DataFrame({
        "sequence":  ["PEP1", "ENT1", "PEP2", "ENT2"],
        "qvalue":    [0.001, 0.02, 0.005, 0.04],
        "is_target": [True, True, True, True],
        "main_score":[100.0, 80.0, 95.0, 70.0],
    })
    results_path = tmp_path / "r.parquet"
    results.write_parquet(results_path)
    out = analyse(
        results_parquet=results_path,
        target_peptides=target_peps,
        entrapment_peptides=entrap_peps,
        ratio=1.0,
        pairing_path=None,
        out_parquet=tmp_path / "c.parquet",
        out_fdr_plot=tmp_path / "fdr.png",
        out_hist_plot=tmp_path / "hist.png",
    )
    # Required scalars
    assert "entrap/n_target_at_q01" in out
    assert "entrap/n_entrap_at_q01" in out
    assert "entrap/empirical_fdr_lower_at_q01" in out
    assert "entrap/empirical_fdr_combined_at_q01" in out
    assert "entrap/ratio" in out
    assert out["entrap/ratio"] == 1.0
    # Files exist
    assert (tmp_path / "c.parquet").exists()
    assert (tmp_path / "fdr.png").exists()
    assert (tmp_path / "hist.png").exists()


def test_analyse_with_pairing_emits_matched(tmp_path):
    target_peps = tmp_path / "t.txt"
    target_peps.write_text("PEP1\nPEP2\n")
    entrap_peps = tmp_path / "e.txt"
    entrap_peps.write_text("ENT1\nENT2\n")
    pairs_path = tmp_path / "pairs.tsv"
    pairs_path.write_text("target_peptide\tentrap_peptide\nPEP1\tENT1\nPEP2\tENT2\n")

    results = pl.DataFrame({
        "sequence":  ["PEP1", "ENT1", "PEP2", "ENT2"],
        "qvalue":    [0.001, 0.002, 0.003, 0.004],
        "is_target": [True, True, True, True],
        "main_score":[100.0, 200.0, 150.0, 50.0],
    })
    results_path = tmp_path / "r.parquet"
    results.write_parquet(results_path)
    out = analyse(
        results_parquet=results_path,
        target_peptides=target_peps,
        entrapment_peptides=entrap_peps,
        ratio=1.0,
        pairing_path=pairs_path,
        out_parquet=tmp_path / "c.parquet",
        out_fdr_plot=tmp_path / "fdr.png",
        out_hist_plot=tmp_path / "hist.png",
    )
    assert "entrap/empirical_fdr_matched_at_q01" in out
