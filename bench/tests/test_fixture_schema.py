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
        target_peptides = "s3://b/peps.txt"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.name == "hela"
    assert f.inputs.target_peptides == "s3://b/peps.txt"
    assert f.inputs.entrapment_peptides is None
    assert f.inputs.entrapment_ratio is None
    assert f.inputs.entrapment_mode is None
    assert f.inputs.pairing is None
    assert f.inputs.calibration_speclib is None
    assert not f.has_entrapment()
    assert not f.has_pairing()


def test_full_entrapment_fixture(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "h_y"
        description = "test"

        [inputs]
        target_peptides = "s3://b/t.txt"
        entrapment_peptides = "s3://b/e.txt"
        entrapment_ratio = 1.0
        entrapment_mode = "shuffled"
        pairing = "s3://b/pairs.tsv"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.has_entrapment()
    assert f.has_pairing()
    assert f.inputs.entrapment_ratio == 1.0
    assert f.inputs.entrapment_mode == "shuffled"


def test_foreign_mode_no_pairing(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "h_y"
        description = "test"

        [inputs]
        target_peptides = "s3://b/t.txt"
        entrapment_peptides = "s3://b/e.txt"
        entrapment_ratio = 1.0
        entrapment_mode = "foreign"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.has_entrapment()
    assert not f.has_pairing()
    assert f.inputs.entrapment_mode == "foreign"


def test_entrapment_peptides_without_ratio_rejected(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "bad"
        description = "test"

        [inputs]
        target_peptides = "s3://b/t.txt"
        entrapment_peptides = "s3://b/e.txt"
        entrapment_mode = "foreign"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    with pytest.raises(ValueError, match="entrapment_ratio"):
        load_fixture(p)


def test_entrapment_peptides_without_mode_rejected(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "bad"
        description = "test"

        [inputs]
        target_peptides = "s3://b/t.txt"
        entrapment_peptides = "s3://b/e.txt"
        entrapment_ratio = 1.0
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    with pytest.raises(ValueError, match="entrapment_mode"):
        load_fixture(p)


def test_orphan_ratio_rejected(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "bad"
        description = "test"

        [inputs]
        target_peptides = "s3://b/t.txt"
        entrapment_ratio = 1.0
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    with pytest.raises(ValueError, match="without entrapment_peptides"):
        load_fixture(p)


def test_pairing_only_for_shuffled_mode(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "bad"
        description = "test"

        [inputs]
        target_peptides = "s3://b/t.txt"
        entrapment_peptides = "s3://b/e.txt"
        entrapment_ratio = 1.0
        entrapment_mode = "foreign"
        pairing = "s3://b/pairs.tsv"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    with pytest.raises(ValueError, match="pairing"):
        load_fixture(p)


def test_ratio_below_one_rejected(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "bad"
        description = "test"

        [inputs]
        target_peptides = "s3://b/t.txt"
        entrapment_peptides = "s3://b/e.txt"
        entrapment_ratio = 0.5
        entrapment_mode = "foreign"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    with pytest.raises(ValueError, match="entrapment_ratio must be >= 1.0"):
        load_fixture(p)


def test_relative_path_in_inputs_rejected(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "bad"
        description = "test"

        [inputs]
        target_peptides = "relative/path.txt"
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    with pytest.raises(ValueError, match="absolute local path"):
        load_fixture(p)


def test_absolute_local_path_accepted(tmp_path):
    p = _write(
        tmp_path,
        """
        name = "ok"
        description = "test"

        [inputs]
        target_peptides = "/abs/path/peps.txt"
        speclib = "/abs/path/lib.msgpack.zst"
        raw = "/abs/path/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.inputs.target_peptides == "/abs/path/peps.txt"


def test_tilde_path_expanded(tmp_path):
    import os
    home = os.path.expanduser("~")
    p = _write(
        tmp_path,
        f"""
        name = "ok"
        description = "test"

        [inputs]
        target_peptides = "~/peps.txt"
        speclib = "{home}/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.inputs.target_peptides == f"{home}/peps.txt"


def test_minimal_fixture_without_target_peptides(tmp_path):
    """Non-entrap fixtures don't need target_peptides."""
    p = _write(
        tmp_path,
        """
        name = "hela"
        description = "test"

        [inputs]
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"

        [config.analysis]
        chunk_size = 20000
        """,
    )
    f = load_fixture(p)
    assert f.inputs.target_peptides is None
    assert not f.has_entrapment()


def test_entrapment_without_target_peptides_rejected(tmp_path):
    """If entrapment_peptides is set, target_peptides must be too."""
    p = _write(
        tmp_path,
        """
        name = "bad"
        description = "test"

        [inputs]
        speclib = "s3://b/lib.msgpack.zst"
        raw = "s3://b/sample.d"
        entrapment_peptides = "s3://b/e.txt"
        entrapment_ratio = 1.0
        entrapment_mode = "foreign"

        [config.analysis]
        chunk_size = 1
        """,
    )
    with pytest.raises(ValueError, match="target_peptides"):
        load_fixture(p)
