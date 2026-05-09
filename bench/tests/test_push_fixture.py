"""Tests for bench/push_fixture.py.

S3 + subprocess + resolve_dbs are always mocked via the fake_runtime fixture.
Digest + shuffle helpers run for real so algorithm correctness is verified.
"""

from __future__ import annotations

import textwrap
import tomllib
import warnings
from pathlib import Path
from unittest.mock import patch

import pytest

from bench._fixture_schema import load_fixture
from bench.push_fixture import (
    build_fixture_toml,
    parse_args,
    run_pipeline,
)

# ---------------------------------------------------------------------------
# parse_args
# ---------------------------------------------------------------------------


def test_parse_args_minimal():
    args = parse_args([
        "--name", "hela",
        "--bucket", "bk",
        "--prefix", "fx",
        "--db", "UP000005640",
        "--raw", "/tmp/sample.d",
        "--config", "/tmp/cfg.toml",
    ])
    assert args.name == "hela"
    assert args.bucket == "bk"
    assert args.prefix == "fx"
    assert args.db == ["UP000005640"]
    assert args.entrap_db == []
    assert args.calib_db == []
    assert args.dry_run is False
    assert args.overwrite is False
    assert args.entrap_ratio == 1.0
    assert args.peptide_min_len == 7
    assert args.peptide_max_len == 30
    assert args.missed_cleavages == 1
    assert args.seed == 42


def test_parse_args_all_flags():
    args = parse_args([
        "--name", "hy",
        "--bucket", "bk",
        "--prefix", "fx",
        "--db", "UP000005640",
        "--db", "s3://bkt/extra.fasta",
        "--entrap-db", "UP000002311",
        "--calib-db", "P12345",
        "--raw", "/tmp/sample.d",
        "--config", "/tmp/cfg.toml",
        "--entrap-ratio", "2.0",
        "--peptide-min-len", "6",
        "--peptide-max-len", "25",
        "--missed-cleavages", "2",
        "--seed", "7",
        "--request-delay-ms", "200",
        "--dry-run",
        "--overwrite",
        "--force",
    ])
    assert args.db == ["UP000005640", "s3://bkt/extra.fasta"]
    assert args.entrap_db == ["UP000002311"]
    assert args.calib_db == ["P12345"]
    assert args.entrap_ratio == 2.0
    assert args.peptide_min_len == 6
    assert args.peptide_max_len == 25
    assert args.missed_cleavages == 2
    assert args.seed == 7
    assert args.request_delay_ms == 200
    assert args.dry_run is True
    assert args.overwrite is True
    assert args.force is True


def test_parse_args_shuffled_entrap_db():
    args = parse_args([
        "--name", "s",
        "--bucket", "bk",
        "--prefix", "fx",
        "--db", "UP000005640",
        "--raw", "/tmp/sample.d",
        "--config", "/tmp/cfg.toml",
        "--entrap-db", "SHUFFLED",
    ])
    assert args.entrap_db == ["SHUFFLED"]


def test_parse_args_ratio_lt_1_rejected():
    with pytest.raises(SystemExit):
        parse_args([
            "--name", "x",
            "--bucket", "bk",
            "--prefix", "fx",
            "--db", "UP000005640",
            "--raw", "/tmp/s.d",
            "--config", "/tmp/cfg.toml",
            "--entrap-ratio", "0.5",
        ])


def test_parse_args_shuffled_mixed_with_foreign_rejected():
    with pytest.raises(SystemExit):
        parse_args([
            "--name", "x",
            "--bucket", "bk",
            "--prefix", "fx",
            "--db", "UP000005640",
            "--raw", "/tmp/s.d",
            "--config", "/tmp/cfg.toml",
            "--entrap-db", "SHUFFLED",
            "--entrap-db", "UP000002311",
        ])


# ---------------------------------------------------------------------------
# build_fixture_toml
# ---------------------------------------------------------------------------


def test_build_fixture_toml_minimal(tmp_path: Path):
    cfg = tmp_path / "cfg.toml"
    cfg.write_text("[analysis]\nchunk_size = 20000\n")
    out = build_fixture_toml(
        name="hela",
        description="200ng HeLa",
        config_path=cfg,
        target_peptides_uri="s3://bk/fx/hela/target.peptides.txt",
        speclib_uri="s3://bk/fx/hela/lib.msgpack.zst",
        raw_uri="s3://bk/fx/hela/sample.d",
    )
    target_path = tmp_path / "fx.toml"
    target_path.write_text(out)
    fx = load_fixture(target_path)
    assert fx.name == "hela"
    assert fx.inputs.entrapment_peptides is None
    assert fx.inputs.entrapment_ratio is None
    assert fx.inputs.entrapment_mode is None
    assert fx.inputs.pairing is None
    assert fx.config["analysis"]["chunk_size"] == 20000


def test_build_fixture_toml_shuffled_with_pairing(tmp_path: Path):
    cfg = tmp_path / "cfg.toml"
    cfg.write_text("[analysis]\nchunk_size = 1\n")
    out = build_fixture_toml(
        name="s",
        description="",
        config_path=cfg,
        target_peptides_uri="s3://b/s/target.peptides.txt",
        speclib_uri="s3://b/s/lib.msgpack.zst",
        raw_uri="s3://b/s/sample.d",
        entrapment_peptides_uri="s3://b/s/entrap.peptides.txt",
        entrapment_ratio=1.0,
        entrapment_mode="shuffled",
        pairing_uri="s3://b/s/pairing.tsv",
    )
    p = tmp_path / "fx.toml"
    p.write_text(out)
    fx = load_fixture(p)
    assert fx.has_entrapment()
    assert fx.has_pairing()
    assert fx.inputs.entrapment_mode == "shuffled"
    assert fx.inputs.entrapment_ratio == 1.0


def test_build_fixture_toml_foreign_with_calib(tmp_path: Path):
    cfg = tmp_path / "cfg.toml"
    cfg.write_text("[analysis]\nchunk_size = 1\n")
    out = build_fixture_toml(
        name="x",
        description="x",
        config_path=cfg,
        target_peptides_uri="s3://b/x/target.peptides.txt",
        speclib_uri="s3://b/x/lib.msgpack.zst",
        raw_uri="s3://b/x/sample.d",
        entrapment_peptides_uri="s3://b/x/entrap.peptides.txt",
        entrapment_ratio=2.0,
        entrapment_mode="foreign",
        calibration_speclib_uri="s3://b/x/calib_lib.msgpack.zst",
    )
    p = tmp_path / "fx.toml"
    p.write_text(out)
    fx = load_fixture(p)
    assert fx.has_entrapment()
    assert not fx.has_pairing()
    assert fx.inputs.entrapment_mode == "foreign"
    assert fx.has_calibration_speclib()


# ---------------------------------------------------------------------------
# fake_runtime fixture
# ---------------------------------------------------------------------------


def _stub_fasta_content() -> str:
    """Three proteins that digest into distinct tryptic peptides."""
    return textwrap.dedent("""\
        >sp|P1|PROT1
        MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEK
        >sp|P2|PROT2
        MAKQADSVSVKAEQYLSAELREQNLAKMSAAEERNRIAESQRQLAEQQKQLEQLKQKLEQLKQKLEQLK
        >sp|P3|PROT3
        MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY
    """)


@pytest.fixture
def fake_runtime(tmp_path):
    """Patch S3 + speclib_build_cli. resolve_dbs writes a real stub fasta so
    digest helpers run for real. Also patches generate_shuffled_entrapment
    is NOT mocked — the real implementation runs with deterministic seed."""
    fx_dir = tmp_path / "bench_fixtures"
    fx_dir.mkdir()

    def _resolve(specs, out_path):
        out_path.write_text(_stub_fasta_content())

    with (
        patch("bench.push_fixture.s3_upload_file") as up_file,
        patch("bench.push_fixture.s3_upload_dir") as up_dir,
        patch("bench.push_fixture.run_speclib_build") as build,
        patch("bench.push_fixture.resolve_dbs", side_effect=_resolve) as res,
    ):
        yield {
            "up_file": up_file,
            "up_dir": up_dir,
            "build": build,
            "res": res,
            "fx_dir": fx_dir,
        }


def _common_args(tmp_path):
    cfg = tmp_path / "cfg.toml"
    cfg.write_text("[analysis]\nchunk_size = 20000\n")
    raw = tmp_path / "sample.d"
    raw.mkdir()
    (raw / "metadata").write_bytes(b"x")
    return cfg, raw


# ---------------------------------------------------------------------------
# run_pipeline — no entrapment
# ---------------------------------------------------------------------------


def test_run_pipeline_minimal(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "hela.toml"
    run_pipeline(
        name="hela",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=[],
        calib_db=[],
        speclib_uri=None,
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=output_toml,
        overwrite=False,
        dry_run=False,
    )
    # resolve_dbs called once for target
    assert fake_runtime["res"].call_count == 1
    # target.peptides.txt + database.peptides.txt uploaded
    uploaded_uris = [c.args[1] for c in fake_runtime["up_file"].call_args_list]
    assert any("target.peptides.txt" in u for u in uploaded_uris)
    assert any("database.peptides.txt" in u for u in uploaded_uris)
    # entrap.peptides.txt NOT uploaded
    assert not any("entrap.peptides.txt" in u for u in uploaded_uris)
    # Raw dir uploaded
    fake_runtime["up_dir"].assert_called_once()
    assert fake_runtime["up_dir"].call_args.args[1] == "s3://bk/fx/hela/sample.d"
    # speclib build called once
    fake_runtime["build"].assert_called_once()
    # speclib built from database.peptides.txt
    build_call_args = fake_runtime["build"].call_args.args
    assert "database.peptides.txt" in build_call_args[0]
    # fixture written
    assert output_toml.exists()
    fx = load_fixture(output_toml)
    assert fx.name == "hela"
    assert fx.inputs.entrapment_peptides is None


def test_run_pipeline_skip_build_when_speclib_provided(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "hela.toml"
    run_pipeline(
        name="hela",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=[],
        calib_db=[],
        speclib_uri="s3://other/lib.msgpack.zst",
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=output_toml,
        overwrite=False,
        dry_run=False,
    )
    fake_runtime["build"].assert_not_called()
    assert "s3://other/lib.msgpack.zst" in output_toml.read_text()


# ---------------------------------------------------------------------------
# run_pipeline — foreign entrapment (Algorithm 2)
# ---------------------------------------------------------------------------


def test_run_pipeline_foreign_entrap(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "x.toml"

    # Use distinct proteins so entrap peptides survive the target-subtraction step
    # and actual_r stays >= 1.0 (passes the schema validator).
    target_fasta = textwrap.dedent("""\
        >sp|T1|TARGET
        MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFK
    """)
    entrap_fasta = textwrap.dedent("""\
        >sp|E1|ENTRAP
        MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY
        >sp|E2|ENTRAP2
        MAKQADSVSVKAEQYLSAELREQNLAKMSAAEERNRIAESQR
        >sp|E3|ENTRAP3
        MALPVTALLLPLALLLHAARPSFSLVKRGELKPAPKALPEPKPAPKALPEPKPAPKALPEPKPVSKMAPP
    """)

    call_count = [0]

    def _resolve_distinct(specs, out_path):
        call_count[0] += 1
        if call_count[0] == 1:
            out_path.write_text(target_fasta)
        else:
            out_path.write_text(entrap_fasta)

    fake_runtime["res"].side_effect = _resolve_distinct

    run_pipeline(
        name="x",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=["UP000002311"],
        calib_db=[],
        speclib_uri=None,
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=output_toml,
        overwrite=False,
        dry_run=False,
    )
    # resolve_dbs called for target + entrap = 2
    assert fake_runtime["res"].call_count == 2
    uploaded_uris = [c.args[1] for c in fake_runtime["up_file"].call_args_list]
    assert any("entrap.peptides.txt" in u for u in uploaded_uris)
    # pairing.tsv NOT uploaded (foreign mode, not shuffled)
    assert not any("pairing.tsv" in u for u in uploaded_uris)
    fx = load_fixture(output_toml)
    assert fx.inputs.entrapment_mode == "foreign"
    assert fx.inputs.pairing is None


def test_run_pipeline_with_calib_db(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "x.toml"
    run_pipeline(
        name="x",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=[],
        calib_db=["P12345"],
        speclib_uri=None,
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=output_toml,
        overwrite=False,
        dry_run=False,
    )
    # resolve_dbs: target + calib = 2
    assert fake_runtime["res"].call_count == 2
    # build called twice: main + calib
    assert fake_runtime["build"].call_count == 2
    uploaded_uris = [c.args[1] for c in fake_runtime["up_file"].call_args_list]
    assert any("calib.peptides.txt" in u for u in uploaded_uris)
    fx = load_fixture(output_toml)
    assert fx.has_calibration_speclib()


def test_run_pipeline_foreign_and_calib(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "x.toml"
    run_pipeline(
        name="x",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=["UP000002311"],
        calib_db=["P12345"],
        speclib_uri=None,
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=output_toml,
        overwrite=False,
        dry_run=False,
    )
    # resolve_dbs: target + entrap + calib = 3
    assert fake_runtime["res"].call_count == 3
    # build called twice
    assert fake_runtime["build"].call_count == 2
    body = output_toml.read_text()
    assert "entrap" in body and "calibration_speclib" in body


# ---------------------------------------------------------------------------
# run_pipeline — shuffled entrapment (Algorithm 1)
# ---------------------------------------------------------------------------


def test_run_pipeline_shuffled(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "s.toml"
    run_pipeline(
        name="s",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=["SHUFFLED"],
        calib_db=[],
        speclib_uri=None,
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=output_toml,
        overwrite=False,
        dry_run=False,
        entrap_ratio=1.0,
        seed=42,
    )
    # resolve_dbs called only for target (SHUFFLED is not resolved via resolve_dbs)
    assert fake_runtime["res"].call_count == 1
    uploaded_uris = [c.args[1] for c in fake_runtime["up_file"].call_args_list]
    # pairing.tsv uploaded (r=1, shuffled)
    assert any("pairing.tsv" in u for u in uploaded_uris)
    assert any("entrap.peptides.txt" in u for u in uploaded_uris)
    fx = load_fixture(output_toml)
    assert fx.inputs.entrapment_mode == "shuffled"
    assert fx.inputs.entrapment_ratio == 1.0
    assert fx.has_pairing()


def test_run_pipeline_shuffled_r2_no_pairing(tmp_path, fake_runtime):
    """r=2 shuffled: two shuffles per target; pairing.tsv NOT emitted."""
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "s2.toml"
    run_pipeline(
        name="s2",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=["SHUFFLED"],
        calib_db=[],
        speclib_uri=None,
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=output_toml,
        overwrite=False,
        dry_run=False,
        entrap_ratio=2.0,
        seed=42,
    )
    uploaded_uris = [c.args[1] for c in fake_runtime["up_file"].call_args_list]
    assert not any("pairing.tsv" in u for u in uploaded_uris)
    fx = load_fixture(output_toml)
    assert fx.inputs.entrapment_mode == "shuffled"
    assert not fx.has_pairing()


def test_run_pipeline_shuffled_mixing_foreign_rejected(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "bad.toml"
    with pytest.raises(ValueError, match="SHUFFLED"):
        run_pipeline(
            name="bad",
            bucket="bk",
            prefix="fx",
            db=["UP000005640"],
            raw=str(raw),
            config=str(cfg),
            entrap_db=["SHUFFLED", "UP000002311"],
            calib_db=[],
            speclib_uri=None,
            calibration_speclib_uri=None,
            koina_url=None,
            fixture_target=output_toml,
            overwrite=False,
            dry_run=False,
        )


# ---------------------------------------------------------------------------
# run_pipeline — overwrite / dry_run / force guards
# ---------------------------------------------------------------------------


def test_run_pipeline_refuses_overwrite(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    target = fake_runtime["fx_dir"] / "hela.toml"
    target.write_text("# pre-existing")

    with pytest.raises(FileExistsError):
        run_pipeline(
            name="hela",
            bucket="bk",
            prefix="fx",
            db=["UP000005640"],
            raw=str(raw),
            config=str(cfg),
            entrap_db=[],
            calib_db=[],
            speclib_uri=None,
            calibration_speclib_uri=None,
            koina_url=None,
            fixture_target=target,
            overwrite=False,
            dry_run=False,
        )


def test_run_pipeline_dry_run(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    target = fake_runtime["fx_dir"] / "hela.toml"
    run_pipeline(
        name="hela",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=[],
        calib_db=[],
        speclib_uri=None,
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=target,
        overwrite=False,
        dry_run=True,
    )
    fake_runtime["res"].assert_not_called()
    fake_runtime["up_file"].assert_not_called()
    fake_runtime["up_dir"].assert_not_called()
    fake_runtime["build"].assert_not_called()
    assert not target.exists()


def test_run_pipeline_default_skips_existing_uploads(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "hela.toml"
    run_pipeline(
        name="hela",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=[],
        calib_db=[],
        speclib_uri=None,
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=output_toml,
        overwrite=False,
        dry_run=False,
    )
    for c in fake_runtime["up_file"].call_args_list:
        assert c.kwargs.get("skip_if_exists") is True
    assert fake_runtime["up_dir"].call_args.kwargs.get("idempotent") is True


def test_run_pipeline_force_overrides_skip(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "hela.toml"
    run_pipeline(
        name="hela",
        bucket="bk",
        prefix="fx",
        db=["UP000005640"],
        raw=str(raw),
        config=str(cfg),
        entrap_db=[],
        calib_db=[],
        speclib_uri=None,
        calibration_speclib_uri=None,
        koina_url=None,
        fixture_target=output_toml,
        overwrite=False,
        dry_run=False,
        force=True,
    )
    for c in fake_runtime["up_file"].call_args_list:
        assert c.kwargs.get("skip_if_exists") is False
    assert fake_runtime["up_dir"].call_args.kwargs.get("idempotent") is False


# ---------------------------------------------------------------------------
# Algorithm verification — actual r recorded correctly
# ---------------------------------------------------------------------------


def test_run_pipeline_foreign_insufficient_warns(tmp_path, fake_runtime):
    """When foreign peptides are fewer than needed, warn and record actual r.

    Note: if actual_r < 1.0 the schema validator rejects it, so we check the
    raw TOML text rather than loading through the schema.
    """
    cfg, raw = _common_args(tmp_path)
    output_toml = fake_runtime["fx_dir"] / "short.toml"

    # target gets a long protein with many tryptic peptides
    target_fasta = textwrap.dedent("""\
        >sp|T1|PROT
        MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEK
    """)
    # Tiny entrap — too short for 7-mer tryptic peptide → 0 foreign after filter
    entrap_fasta = textwrap.dedent("""\
        >sp|E1|ENTRAP
        MKWVTFISLLL
    """)

    call_count = [0]

    def _resolve_custom(specs, out_path):
        call_count[0] += 1
        if call_count[0] == 1:
            out_path.write_text(target_fasta)
        else:
            out_path.write_text(entrap_fasta)

    fake_runtime["res"].side_effect = _resolve_custom

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        run_pipeline(
            name="short",
            bucket="bk",
            prefix="fx",
            db=["UP000005640"],
            raw=str(raw),
            config=str(cfg),
            entrap_db=["UP000002311"],
            calib_db=[],
            speclib_uri=None,
            calibration_speclib_uri=None,
            koina_url=None,
            fixture_target=output_toml,
            overwrite=False,
            dry_run=False,
            entrap_ratio=5.0,
        )
    # Should have warned about insufficient peptides
    assert any("Not enough foreign" in str(warning.message) for warning in w)
    # actual_r in the raw TOML should be < 5.0
    toml_text = output_toml.read_text()
    assert "entrapment_ratio" in toml_text
    data = tomllib.loads(toml_text)
    actual_r = data["inputs"]["entrapment_ratio"]
    assert actual_r < 5.0
