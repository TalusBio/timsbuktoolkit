import textwrap
from pathlib import Path
from unittest.mock import patch

import pytest


def _write_fixture(
    dir: Path,
    name: str,
    *,
    with_entrap: bool = False,
    with_pairing: bool = False,
    with_calib: bool = False,
) -> Path:
    extras = ""
    if with_entrap:
        extras += '\nentrapment_peptides = "s3://b/entrap.peptides.txt"'
        extras += '\nentrapment_ratio = 1.0'
        mode = "shuffled" if with_pairing else "foreign"
        extras += f'\nentrapment_mode = "{mode}"'
    if with_pairing:
        extras += '\npairing = "s3://b/pairs.tsv"'
    if with_calib:
        extras += '\ncalibration_speclib = "s3://b/calib.msgpack.zst"'
    p = dir / f"{name}.toml"
    p.write_text(
        textwrap.dedent(
            f"""
            name = "{name}"
            description = "x"

            [inputs]
            target_peptides = "s3://b/target.peptides.txt"
            speclib = "s3://b/lib.msgpack.zst"
            raw = "s3://b/sample.d"{extras}

            [config.analysis]
            chunk_size = 20000
            """
        ).strip()
    )
    return p


@pytest.fixture
def fake_s3(tmp_path):
    """Patch S3 ops; downloads create stub files at the requested local paths."""
    def fake_download_file(uri, dst):
        Path(dst).parent.mkdir(parents=True, exist_ok=True)
        Path(dst).write_text(f"# stub of {uri}\n")

    def fake_sync(uri, dst):
        Path(dst).mkdir(parents=True, exist_ok=True)
        (Path(dst) / "metadata").write_text(f"# stub of {uri}\n")

    with (
        patch(
            "bench.stage_fixture.s3_download_file", side_effect=fake_download_file
        ) as dl,
        patch("bench.stage_fixture.s3_sync_dir", side_effect=fake_sync) as sync,
    ):
        yield {"download": dl, "sync": sync}


def test_stage_minimal(tmp_path, fake_s3):
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    _write_fixture(fx_dir, "hela")
    cache = tmp_path / "cache"
    out = tmp_path / "staged" / "hela.toml"

    from bench.stage_fixture import stage

    stage(
        name="hela",
        fixtures_dir=fx_dir,
        cache_dir=cache,
        out=out,
        overwrite=False,
        force=False,
    )

    # Files were "downloaded"
    assert (cache / "hela" / "target.peptides.txt").exists()
    assert (cache / "hela" / "lib.msgpack.zst").exists()
    assert (cache / "hela" / "sample.d" / "metadata").exists()

    # download_file called for target_peptides + speclib (2x); sync called for raw (1x)
    assert fake_s3["download"].call_count == 2
    assert fake_s3["sync"].call_count == 1

    # Output TOML exists and is valid against schema
    assert out.exists()
    from bench._fixture_schema import load_fixture
    fx = load_fixture(out)
    assert fx.name == "hela"
    assert fx.inputs.target_peptides == str(cache / "hela" / "target.peptides.txt")
    assert fx.inputs.speclib == str(cache / "hela" / "lib.msgpack.zst")
    assert fx.inputs.raw == str(cache / "hela" / "sample.d")


def test_stage_with_entrap_and_calib(tmp_path, fake_s3):
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    _write_fixture(fx_dir, "hy", with_entrap=True, with_calib=True)
    cache = tmp_path / "cache"
    out = tmp_path / "staged" / "hy.toml"

    from bench.stage_fixture import stage

    stage(
        name="hy",
        fixtures_dir=fx_dir,
        cache_dir=cache,
        out=out,
        overwrite=False,
        force=False,
    )
    # 4 downloads (target_peptides, speclib, entrap_peptides, calib); 1 sync (raw)
    assert fake_s3["download"].call_count == 4
    assert fake_s3["sync"].call_count == 1

    from bench._fixture_schema import load_fixture
    fx = load_fixture(out)
    assert fx.has_entrapment()
    assert fx.has_calibration_speclib()
    assert fx.inputs.entrapment_peptides == str(cache / "hy" / "entrap.peptides.txt")
    assert fx.inputs.entrapment_ratio == 1.0
    assert fx.inputs.entrapment_mode == "foreign"
    assert fx.inputs.calibration_speclib == str(cache / "hy" / "calib.msgpack.zst")


def test_stage_with_pairing(tmp_path, fake_s3):
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    _write_fixture(fx_dir, "sh", with_entrap=True, with_pairing=True)
    cache = tmp_path / "cache"
    out = tmp_path / "staged" / "sh.toml"

    from bench.stage_fixture import stage

    stage(
        name="sh",
        fixtures_dir=fx_dir,
        cache_dir=cache,
        out=out,
        overwrite=False,
        force=False,
    )
    # 4 file downloads: target_peptides, entrap_peptides, pairing, speclib; 1 sync (raw)
    assert fake_s3["download"].call_count == 4
    assert fake_s3["sync"].call_count == 1

    from bench._fixture_schema import load_fixture
    fx = load_fixture(out)
    assert fx.has_entrapment()
    assert fx.has_pairing()
    assert fx.inputs.entrapment_mode == "shuffled"
    assert fx.inputs.pairing == str(cache / "sh" / "pairing.tsv")


def test_stage_skips_existing_local_files(tmp_path, fake_s3):
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    _write_fixture(fx_dir, "hela")
    cache = tmp_path / "cache"
    # Pre-create target_peptides + speclib (raw is always synced — sync is idempotent)
    (cache / "hela").mkdir(parents=True)
    (cache / "hela" / "target.peptides.txt").write_text("preexisting")
    (cache / "hela" / "lib.msgpack.zst").write_text("preexisting")
    out = tmp_path / "staged" / "hela.toml"

    from bench.stage_fixture import stage

    stage(
        name="hela",
        fixtures_dir=fx_dir,
        cache_dir=cache,
        out=out,
        overwrite=False,
        force=False,
    )
    # No download calls because local files already exist
    assert fake_s3["download"].call_count == 0
    # Sync still runs (raw .d) — sync itself is idempotent
    assert fake_s3["sync"].call_count == 1
    # Local target_peptides content was NOT overwritten
    assert (cache / "hela" / "target.peptides.txt").read_text() == "preexisting"


def test_stage_force_redownloads(tmp_path, fake_s3):
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    _write_fixture(fx_dir, "hela")
    cache = tmp_path / "cache"
    (cache / "hela").mkdir(parents=True)
    (cache / "hela" / "target.peptides.txt").write_text("preexisting")
    out = tmp_path / "staged" / "hela.toml"

    from bench.stage_fixture import stage

    stage(
        name="hela",
        fixtures_dir=fx_dir,
        cache_dir=cache,
        out=out,
        overwrite=False,
        force=True,
    )
    # Force re-downloads even if local exists
    assert fake_s3["download"].call_count == 2
    # Stub overwrites the preexisting file
    assert (cache / "hela" / "target.peptides.txt").read_text().startswith("# stub")


def test_stage_preserves_local_paths(tmp_path, fake_s3):
    """Inputs that are already absolute local paths are referenced as-is, no copy."""
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    p = fx_dir / "lo.toml"
    local_peptides = tmp_path / "abs_target.peptides.txt"
    local_peptides.write_text("MK\nLL\n")
    p.write_text(
        textwrap.dedent(
            f"""
            name = "lo"
            description = "x"

            [inputs]
            target_peptides = "{local_peptides}"
            speclib = "s3://b/lib.msgpack.zst"
            raw = "s3://b/sample.d"

            [config.analysis]
            chunk_size = 1
            """
        ).strip()
    )
    cache = tmp_path / "cache"
    out = tmp_path / "staged" / "lo.toml"

    from bench.stage_fixture import stage

    stage(
        name="lo",
        fixtures_dir=fx_dir,
        cache_dir=cache,
        out=out,
        overwrite=False,
        force=False,
    )
    # Only the speclib gets downloaded; target_peptides is referenced as-is
    assert fake_s3["download"].call_count == 1
    from bench._fixture_schema import load_fixture
    fx = load_fixture(out)
    assert fx.inputs.target_peptides == str(local_peptides)


def test_stage_refuses_overwrite(tmp_path, fake_s3):
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    _write_fixture(fx_dir, "hela")
    cache = tmp_path / "cache"
    out = tmp_path / "staged" / "hela.toml"
    out.parent.mkdir(parents=True)
    out.write_text("# preexisting")

    from bench.stage_fixture import stage

    with pytest.raises(FileExistsError):
        stage(
            name="hela",
            fixtures_dir=fx_dir,
            cache_dir=cache,
            out=out,
            overwrite=False,
            force=False,
        )


def test_parse_args_defaults(monkeypatch, tmp_path):
    monkeypatch.delenv("BENCH_CACHE_DIR", raising=False)
    from bench.stage_fixture import parse_args
    args = parse_args(["hela"])
    assert args.name == "hela"
    assert args.fixtures_dir == Path("bench/fixtures")
    assert args.cache_dir == Path("bench_out/cache")
    assert args.out == Path("bench_out/staged/hela.toml")
    assert args.overwrite is False
    assert args.force is False


def test_parse_args_env_var_for_cache(monkeypatch):
    monkeypatch.setenv("BENCH_CACHE_DIR", "/tmp/my_cache")
    from bench.stage_fixture import parse_args
    args = parse_args(["hela"])
    assert args.cache_dir == Path("/tmp/my_cache")
