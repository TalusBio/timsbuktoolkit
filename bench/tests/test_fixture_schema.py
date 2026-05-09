import textwrap

import pytest

from bench._fixture_schema import load_fixture


def _write(tmp_path, content: str):
    p = tmp_path / "fix.toml"
    p.write_text(textwrap.dedent(content).strip())
    return p


def test_minimal_fixture_loads(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "hela"
        description = "test"

        [inputs]
        fasta = "s3://b/p.fasta"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.name == "hela"
    assert f.inputs.fasta == "s3://b/p.fasta"
    assert f.inputs.entrapment_fasta is None
    assert f.inputs.calibration_speclib is None
    assert not f.has_entrapment()
    assert f.config["analysis"]["chunk_size"] == 20000


def test_entrapment_and_calib_optional_present(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "hela_entrap"
        description = "test"

        [inputs]
        fasta = "s3://b/p.fasta"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"
        entrapment_fasta = "s3://b/entrap.fasta"
        calibration_speclib = "s3://b/calib.msgpack.zst"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.has_entrapment()
    assert f.inputs.calibration_speclib == "s3://b/calib.msgpack.zst"
    assert f.has_calibration_speclib()


def test_relative_path_in_inputs_rejected(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "bad"
        description = "test"

        [inputs]
        fasta = "relative/path.fasta"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    with pytest.raises(ValueError, match="absolute local path"):
        load_fixture(p)


def test_absolute_local_path_in_inputs_accepted(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "ok"
        description = "test"

        [inputs]
        fasta = "/abs/path/proteome.fasta"
        speclib = "/abs/path/lib.msgpack.zst"
        raw = "/abs/path/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.inputs.fasta == "/abs/path/proteome.fasta"


def test_tilde_path_expanded(tmp_path):
    import os
    home = os.path.expanduser("~")
    p = _write(
        tmp_path,
        f"""
        name = "ok"
        description = "test"

        [inputs]
        fasta = "~/proteome.fasta"
        speclib = "{home}/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.inputs.fasta == f"{home}/proteome.fasta"
    assert f.inputs.speclib == f"{home}/lib.msgpack.zst"


def test_missing_required_input_rejected(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "missing"
        description = "test"

        [inputs]
        fasta = "s3://b/p.fasta"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    with pytest.raises(ValueError, match="speclib"):
        load_fixture(p)
