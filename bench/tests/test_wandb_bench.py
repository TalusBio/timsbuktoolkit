import textwrap
from pathlib import Path

import pytest


def _write_fx(dir: Path, name: str) -> Path:
    p = dir / f"{name}.toml"
    p.write_text(
        textwrap.dedent(
            f"""
            name = "{name}"
            description = "x"

            [inputs]
            fasta = "s3://b/p.fasta"
            speclib = "s3://b/lib.msgpack.zst"
            raw = "s3://b/sample.d"

            [config.analysis]
            chunk_size = 20000
            """
        ).strip()
    )
    return p


def test_select_positional(tmp_path):
    _write_fx(tmp_path, "hela")
    _write_fx(tmp_path, "yeast")
    from bench.wandb_bench import select_fixtures
    out = select_fixtures(["hela"], all_=False, match=None, fixtures_dir=tmp_path)
    assert [f.name for f in out] == ["hela"]


def test_select_all(tmp_path):
    _write_fx(tmp_path, "a")
    _write_fx(tmp_path, "b")
    from bench.wandb_bench import select_fixtures
    out = select_fixtures([], all_=True, match=None, fixtures_dir=tmp_path)
    assert sorted(f.name for f in out) == ["a", "b"]


def test_select_match_glob(tmp_path):
    _write_fx(tmp_path, "hela_a")
    _write_fx(tmp_path, "hela_b")
    _write_fx(tmp_path, "yeast_c")
    from bench.wandb_bench import select_fixtures
    out = select_fixtures([], all_=False, match="hela*", fixtures_dir=tmp_path)
    assert sorted(f.name for f in out) == ["hela_a", "hela_b"]


def test_select_unknown_name_errors(tmp_path):
    _write_fx(tmp_path, "hela")
    from bench.wandb_bench import select_fixtures
    with pytest.raises(SystemExit) as exc:
        select_fixtures(["nope"], all_=False, match=None, fixtures_dir=tmp_path)
    assert "nope" in str(exc.value)


def test_select_combinations_error(tmp_path):
    _write_fx(tmp_path, "hela")
    from bench.wandb_bench import select_fixtures
    with pytest.raises(SystemExit):
        select_fixtures(["hela"], all_=True, match=None, fixtures_dir=tmp_path)
    with pytest.raises(SystemExit):
        select_fixtures(["hela"], all_=False, match="*", fixtures_dir=tmp_path)
    with pytest.raises(SystemExit):
        select_fixtures([], all_=True, match="*", fixtures_dir=tmp_path)


def test_select_no_args_lists_available(tmp_path, capsys):
    _write_fx(tmp_path, "hela")
    _write_fx(tmp_path, "yeast")
    from bench.wandb_bench import select_fixtures
    with pytest.raises(SystemExit):
        select_fixtures([], all_=False, match=None, fixtures_dir=tmp_path)
    err = capsys.readouterr().err
    assert "hela" in err and "yeast" in err
