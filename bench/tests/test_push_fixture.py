import textwrap
from pathlib import Path

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
