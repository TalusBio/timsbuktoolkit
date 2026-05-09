import json
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from bench._fixture_schema import load_fixture


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


def _write_perf_report(out_dir: Path, raw_stem: str, payload: dict) -> None:
    """`out_dir` is timsseek's `--output-dir` (i.e., the runner's `res_dir`).
    Timsseek writes its outputs under `<output-dir>/<raw_stem>/`."""
    sub = out_dir / raw_stem
    sub.mkdir(parents=True, exist_ok=True)
    (sub / "performance_report.json").write_text(json.dumps(payload))


@pytest.fixture
def fake_wandb():
    with patch("bench.wandb_bench.wandb") as w:
        run = MagicMock()
        w.init.return_value = run
        yield {"wandb": w, "run": run}


def test_run_one_fixture_logs_perf_report(tmp_path, fake_wandb):
    """timsseek subprocess is mocked; we drop a fake performance_report.json
    where the runner expects it, then assert wandb.run.log got the payload."""
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    fixture = _write_fx(fx_dir, "hela")
    out_root = tmp_path / "out"

    def fake_subprocess(cmd, *a, **kw):
        # The runner passes --output-dir <bench_out/logs/hela-...>/res; locate it
        idx = cmd.index("--output-dir")
        out = Path(cmd[idx + 1])
        # Raw file basename is also passed via --dotd-files
        raw_idx = cmd.index("--dotd-files")
        raw_stem = Path(cmd[raw_idx + 1]).stem
        _write_perf_report(out, raw_stem, {"runtime_s": 12.3, "n_targets_q01": 1000})
        return MagicMock(returncode=0)

    from bench.wandb_bench import run_one
    with patch("bench.wandb_bench.subprocess.run", side_effect=fake_subprocess):
        run_one(
            load_fixture(fixture),
            out_root=out_root,
            notes="hi",
            dry_run=False,
        )

    fake_wandb["wandb"].init.assert_called_once()
    # log should have been invoked with the payload that came from the json file
    fake_wandb["run"].log.assert_called()
    logged = fake_wandb["run"].log.call_args.args[0]
    assert logged["runtime_s"] == 12.3
    assert logged["n_targets_q01"] == 1000


def test_run_one_fixture_runs_entrapment_when_field_present(tmp_path, fake_wandb):
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    p = fx_dir / "hy.toml"
    p.write_text(
        textwrap.dedent(
            """
            name = "hy"
            description = "x"

            [inputs]
            fasta = "s3://b/p.fasta"
            speclib = "s3://b/lib.msgpack.zst"
            raw = "s3://b/sample.d"
            entrapment_fasta = "s3://b/entrap.fasta"

            [config.analysis]
            chunk_size = 20000
            """
        ).strip()
    )
    out_root = tmp_path / "out"

    def fake_subprocess(cmd, *a, **kw):
        idx = cmd.index("--output-dir")
        out = Path(cmd[idx + 1])
        raw_idx = cmd.index("--dotd-files")
        raw_stem = Path(cmd[raw_idx + 1]).stem
        _write_perf_report(out, raw_stem, {"runtime_s": 1.0})
        # Also drop a results.parquet so analyse() can read it
        import polars as pl
        pl.DataFrame({"sequence": ["MK"], "qvalue": [0.001]}).write_parquet(
            out / raw_stem / "results.parquet"
        )
        return MagicMock(returncode=0)

    from bench.wandb_bench import run_one
    with (
        patch("bench.wandb_bench.subprocess.run", side_effect=fake_subprocess),
        patch("bench.wandb_bench.analyse") as analyse_mock,
        patch("bench.wandb_bench.s3_download_file") as s3_dl,
    ):
        # entrapment.analyse returns scalars; s3_download fetches the two fastas locally
        analyse_mock.return_value = {"entrap/empirical_fdr_at_q01": 0.012}

        def _dl(uri, dst):
            Path(dst).write_text(">x\nMK\n")
        s3_dl.side_effect = _dl

        run_one(load_fixture(p), out_root=out_root, notes=None, dry_run=False)

    analyse_mock.assert_called_once()
    # wandb.run.log should have been called with the entrapment scalars at some point
    log_payloads = [c.args[0] for c in fake_wandb["run"].log.call_args_list]
    assert any("entrap/empirical_fdr_at_q01" in p for p in log_payloads)


def test_run_one_dry_run_no_subprocess(tmp_path, fake_wandb):
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    fixture = _write_fx(fx_dir, "hela")

    from bench.wandb_bench import run_one
    with patch("bench.wandb_bench.subprocess.run") as sp:
        run_one(
            load_fixture(fixture),
            out_root=tmp_path / "out",
            notes=None,
            dry_run=True,
        )
    sp.assert_not_called()
    fake_wandb["wandb"].init.assert_not_called()
