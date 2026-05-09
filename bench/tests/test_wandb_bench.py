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
            target_peptides = "s3://b/target.peptides.txt"
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
            target_peptides = "s3://b/t.peptides.txt"
            entrapment_peptides = "s3://b/e.peptides.txt"
            entrapment_ratio = 1.0
            entrapment_mode = "foreign"
            speclib = "s3://b/lib.msgpack.zst"
            raw = "s3://b/sample.d"

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
        sub = out / raw_stem
        sub.mkdir(parents=True, exist_ok=True)
        (sub / "performance_report.json").write_text(json.dumps({"runtime_s": 1.0}))
        import polars as pl

        pl.DataFrame({
            "sequence":  ["MK", "LL"],
            "qvalue":    [0.001, 0.002],
            "is_target": [True, True],
            "main_score": [100.0, 90.0],
        }).write_parquet(sub / "results.parquet")
        return MagicMock(returncode=0)

    from bench.wandb_bench import run_one

    with (
        patch("bench.wandb_bench.subprocess.run", side_effect=fake_subprocess),
        patch("bench.wandb_bench.analyse") as analyse_mock,
        patch("bench.wandb_bench.s3_download_file") as s3_dl,
    ):
        analyse_mock.return_value = {"entrap/empirical_fdr_combined_at_q01": 0.012}

        def _dl(uri, dst):
            Path(dst).write_text("MK\nLL\n")

        s3_dl.side_effect = _dl

        run_one(load_fixture(p), out_root=out_root, notes=None, dry_run=False)

    analyse_mock.assert_called_once()
    # Verify analyse was called with peptide paths (not fasta) AND ratio
    kwargs = analyse_mock.call_args.kwargs
    assert "target_peptides" in kwargs
    assert "entrapment_peptides" in kwargs
    assert kwargs["ratio"] == 1.0
    assert kwargs["pairing_path"] is None
    assert "out_hist_plot" in kwargs
    # wandb.run.log was called with the entrap scalars + with the histogram image
    log_payloads = [c.args[0] for c in fake_wandb["run"].log.call_args_list]
    assert any("entrap/empirical_fdr_combined_at_q01" in pp for pp in log_payloads)
    assert any("entrap/mainscore_hist" in pp for pp in log_payloads)


def test_run_one_fixture_with_pairing(tmp_path, fake_wandb):
    """When pairing is set, the runner downloads it and forwards path to analyse."""
    fx_dir = tmp_path / "fx"
    fx_dir.mkdir()
    p = fx_dir / "shuffle.toml"
    p.write_text(
        textwrap.dedent(
            """
            name = "shuffle"
            description = "x"

            [inputs]
            target_peptides = "s3://b/t.peptides.txt"
            entrapment_peptides = "s3://b/e.peptides.txt"
            entrapment_ratio = 1.0
            entrapment_mode = "shuffled"
            pairing = "s3://b/pairs.tsv"
            speclib = "s3://b/lib.msgpack.zst"
            raw = "s3://b/sample.d"

            [config.analysis]
            chunk_size = 20000
            """
        ).strip()
    )

    def fake_subprocess(cmd, *a, **kw):
        idx = cmd.index("--output-dir")
        out = Path(cmd[idx + 1])
        raw_idx = cmd.index("--dotd-files")
        raw_stem = Path(cmd[raw_idx + 1]).stem
        sub = out / raw_stem
        sub.mkdir(parents=True, exist_ok=True)
        (sub / "performance_report.json").write_text(json.dumps({"runtime_s": 1.0}))
        import polars as pl

        pl.DataFrame({
            "sequence":  ["MK"],
            "qvalue":    [0.001],
            "is_target": [True],
            "main_score": [100.0],
        }).write_parquet(sub / "results.parquet")
        return MagicMock(returncode=0)

    from bench.wandb_bench import run_one

    with (
        patch("bench.wandb_bench.subprocess.run", side_effect=fake_subprocess),
        patch("bench.wandb_bench.analyse") as analyse_mock,
        patch("bench.wandb_bench.s3_download_file") as s3_dl,
    ):
        analyse_mock.return_value = {"entrap/empirical_fdr_matched_at_q01": 0.005}
        s3_dl.side_effect = lambda uri, dst: Path(dst).write_text("MK\n")
        run_one(load_fixture(p), out_root=tmp_path / "out", notes=None, dry_run=False)

    # 3 downloads: target, entrap, pairing
    assert s3_dl.call_count == 3
    # analyse received a non-None pairing_path
    assert analyse_mock.call_args.kwargs["pairing_path"] is not None


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


def test_parse_args_fixtures_dir_flag():
    from bench.wandb_bench import parse_args
    args = parse_args([
        "--fixtures-dir", "bench_out/staged",
        "myfix",
    ])
    assert args.fixtures_dir == Path("bench_out/staged")
    assert args.fixtures == ["myfix"]


def test_parse_args_fixtures_dir_default():
    from bench.wandb_bench import DEFAULT_FIXTURES_DIR, parse_args
    args = parse_args(["myfix"])
    assert args.fixtures_dir == DEFAULT_FIXTURES_DIR
