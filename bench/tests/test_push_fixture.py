import textwrap
from pathlib import Path
from unittest.mock import ANY, patch

import pytest

from bench.push_fixture import build_fixture_toml, parse_args


def test_parse_args_minimal():
    args = parse_args(
        [
            "--name", "hela",
            "--bucket", "bk",
            "--prefix", "fx",
            "--db", "UP000005640",
            "--raw", "/tmp/sample.d",
            "--config", "/tmp/cfg.toml",
        ]
    )
    assert args.name == "hela"
    assert args.bucket == "bk"
    assert args.prefix == "fx"
    assert args.db == ["UP000005640"]
    assert args.entrap_db == []
    assert args.calib_db == []
    assert args.dry_run is False
    assert args.overwrite is False


def test_parse_args_multiple_db_and_entrap():
    args = parse_args(
        [
            "--name", "hy",
            "--bucket", "bk",
            "--prefix", "fx",
            "--db", "UP000005640",
            "--db", "/tmp/extra.fasta",
            "--entrap-db", "UP000002311",
            "--raw", "/tmp/sample.d",
            "--config", "/tmp/cfg.toml",
            "--dry-run",
        ]
    )
    assert args.db == ["UP000005640", "/tmp/extra.fasta"]
    assert args.entrap_db == ["UP000002311"]
    assert args.dry_run is True


def test_build_fixture_toml(tmp_path: Path):
    cfg = tmp_path / "cfg.toml"
    cfg.write_text(
        textwrap.dedent(
            """
            [analysis]
            chunk_size = 20000
            """
        ).strip()
    )
    out = build_fixture_toml(
        name="hela",
        description="200ng HeLa",
        config_path=cfg,
        fasta_uri="s3://bk/fx/hela/proteome.fasta",
        speclib_uri="s3://bk/fx/hela/lib.msgpack.zst",
        raw_uri="s3://bk/fx/hela/sample.d",
        entrapment_fasta_uri=None,
        calibration_speclib_uri=None,
    )
    # Round-trip via the schema loader to verify validity
    target_path = tmp_path / "fx.toml"
    target_path.write_text(out)
    from bench._fixture_schema import load_fixture
    fx = load_fixture(target_path)
    assert fx.name == "hela"
    assert fx.inputs.entrapment_fasta is None
    assert fx.config["analysis"]["chunk_size"] == 20000


def test_build_fixture_toml_with_entrap_and_calib(tmp_path: Path):
    cfg = tmp_path / "cfg.toml"
    cfg.write_text("[analysis]\nchunk_size = 1\n")
    out = build_fixture_toml(
        name="x",
        description="x",
        config_path=cfg,
        fasta_uri="s3://b/x/proteome.fasta",
        speclib_uri="s3://b/x/lib.msgpack.zst",
        raw_uri="s3://b/x/sample.d",
        entrapment_fasta_uri="s3://b/x/entrap.fasta",
        calibration_speclib_uri="s3://b/x/calib.msgpack.zst",
    )
    p = tmp_path / "fx.toml"
    p.write_text(out)
    from bench._fixture_schema import load_fixture
    fx = load_fixture(p)
    assert fx.has_entrapment()
    assert fx.has_calibration_speclib()


def _common_args(tmp_path):
    cfg = tmp_path / "cfg.toml"
    cfg.write_text("[analysis]\nchunk_size = 20000\n")
    raw = tmp_path / "sample.d"
    raw.mkdir()
    (raw / "metadata").write_bytes(b"x")
    return cfg, raw


@pytest.fixture
def fake_runtime(tmp_path):
    """Patch S3 + speclib_build_cli + resolve_dbs for the run_pipeline tests."""
    fx_dir = tmp_path / "bench_fixtures"
    fx_dir.mkdir()

    with (
        patch("bench.push_fixture.s3_upload_file") as up_file,
        patch("bench.push_fixture.s3_upload_dir") as up_dir,
        patch("bench.push_fixture.run_speclib_build") as build,
        patch("bench.push_fixture.resolve_dbs") as res,
    ):
        # resolve_dbs writes a stub fasta to the requested output path
        def _resolve(specs, out):
            out.write_text(">x\nMK\n")

        res.side_effect = _resolve

        yield {
            "up_file": up_file,
            "up_dir": up_dir,
            "build": build,
            "res": res,
            "fx_dir": fx_dir,
        }


def test_run_pipeline_minimal(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    from bench.push_fixture import run_pipeline

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

    # Resolved + uploaded the target fasta
    assert fake_runtime["res"].call_count == 1
    fake_runtime["up_file"].assert_any_call(ANY, "s3://bk/fx/hela/proteome.fasta")
    # Uploaded the raw directory
    fake_runtime["up_dir"].assert_called_once()
    args = fake_runtime["up_dir"].call_args.args
    assert args[0] == str(raw)
    assert args[1] == "s3://bk/fx/hela/sample.d"
    # Built the speclib
    fake_runtime["build"].assert_called_once()
    # Wrote fixture TOML
    assert output_toml.exists()
    body = output_toml.read_text()
    assert "s3://bk/fx/hela/lib.msgpack.zst" in body


def test_run_pipeline_skip_build_when_speclib_provided(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    from bench.push_fixture import run_pipeline

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
    # No speclib build
    fake_runtime["build"].assert_not_called()
    # Fixture TOML references the user-provided URI
    assert "s3://other/lib.msgpack.zst" in output_toml.read_text()


def test_run_pipeline_with_entrap_and_calib_db(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    from bench.push_fixture import run_pipeline

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
    # resolve_dbs called for db, entrap_db, calib_db = 3 times
    assert fake_runtime["res"].call_count == 3
    # speclib_build called twice: main + calibration
    assert fake_runtime["build"].call_count == 2
    body = output_toml.read_text()
    assert "entrapment_fasta" in body and "calibration_speclib" in body


def test_run_pipeline_refuses_overwrite(tmp_path, fake_runtime):
    cfg, raw = _common_args(tmp_path)
    from bench.push_fixture import run_pipeline

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
    from bench.push_fixture import run_pipeline

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
    # No side effects
    fake_runtime["res"].assert_not_called()
    fake_runtime["up_file"].assert_not_called()
    fake_runtime["up_dir"].assert_not_called()
    fake_runtime["build"].assert_not_called()
    assert not target.exists()
