from pathlib import Path
from unittest.mock import patch

import pytest

from bench._db_resolver import DbSpec, DbSpecKind, classify_db_spec, resolve_dbs


def test_classify_local_fasta(tmp_path):
    p = tmp_path / "x.fasta"
    p.write_text(">a\nMK\n")
    spec = classify_db_spec(str(p))
    assert spec == DbSpec(kind=DbSpecKind.LOCAL_FASTA, value=str(p))


def test_classify_local_fasta_gz(tmp_path):
    p = tmp_path / "x.fasta.gz"
    p.write_bytes(b"\x1f\x8b")
    spec = classify_db_spec(str(p))
    assert spec.kind == DbSpecKind.LOCAL_FASTA


def test_classify_local_accession_list(tmp_path):
    p = tmp_path / "ids.txt"
    p.write_text("P12345\nQ67890\n")
    spec = classify_db_spec(str(p))
    assert spec.kind == DbSpecKind.ACCESSION_LIST_FILE


def test_classify_s3():
    spec = classify_db_spec("s3://bkt/proteome.fasta")
    assert spec.kind == DbSpecKind.S3_FASTA


def test_classify_uniprot_proteome():
    spec = classify_db_spec("UP000005640")
    assert spec.kind == DbSpecKind.UNIPROT_PROTEOME


def test_classify_uniprot_accession():
    spec = classify_db_spec("P12345")
    assert spec.kind == DbSpecKind.UNIPROT_ACCESSION


def test_classify_unknown_raises():
    with pytest.raises(ValueError, match="unrecognised"):
        classify_db_spec("not-a-real-thing")


def test_classify_local_missing_file_raises(tmp_path):
    with pytest.raises(ValueError, match="unrecognised"):
        classify_db_spec(str(tmp_path / "nope.fasta"))


def test_classify_local_bad_extension_raises(tmp_path):
    p = tmp_path / "sequences.parquet"
    p.write_bytes(b"PAR1")
    with pytest.raises(ValueError, match="unrecognised"):
        classify_db_spec(str(p))


def test_resolve_local_fasta_only(tmp_path):
    a = tmp_path / "a.fasta"
    a.write_text(">a\nMK\n")
    b = tmp_path / "b.fasta"
    b.write_text(">b\nLL\n")
    out = tmp_path / "merged.fasta"
    resolve_dbs([str(a), str(b)], out)
    text = out.read_text()
    assert ">a" in text and ">b" in text
    assert text.count(">") == 2


def test_resolve_uniprot_proteome_concat(tmp_path):
    out = tmp_path / "merged.fasta"
    with patch("bench._db_resolver.fetch_proteome", return_value=">P\nMK\n") as m:
        resolve_dbs(["UP000005640"], out)
    m.assert_called_once_with("UP000005640")
    assert out.read_text() == ">P\nMK\n"


def test_resolve_uniprot_accession_batched(tmp_path):
    out = tmp_path / "merged.fasta"
    with patch(
        "bench._db_resolver.fetch_accession_batch", return_value=">A\nM\n>B\nL\n"
    ) as m:
        resolve_dbs(["P12345", "Q67890"], out)
    # Single batch call carries both accessions
    m.assert_called_once_with(["P12345", "Q67890"])
    assert ">A" in out.read_text() and ">B" in out.read_text()


def test_resolve_accession_list_file(tmp_path):
    ids = tmp_path / "ids.txt"
    ids.write_text("P12345\nQ67890\n\n")
    out = tmp_path / "merged.fasta"
    with patch("bench._db_resolver.fetch_accession_batch", return_value=">x\nM\n") as m:
        resolve_dbs([str(ids)], out)
    m.assert_called_once_with(["P12345", "Q67890"])


def test_resolve_s3_uses_aws_cp(tmp_path):
    out = tmp_path / "merged.fasta"

    def fake_s3_download(uri, dst):
        Path(dst).write_text(">s3\nMK\n")

    with patch(
        "bench._db_resolver.s3_download_file", side_effect=fake_s3_download
    ) as m:
        resolve_dbs(["s3://bkt/p.fasta"], out)
    m.assert_called_once()
    assert ">s3" in out.read_text()
