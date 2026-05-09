import pytest

from bench._db_resolver import DbSpec, DbSpecKind, classify_db_spec


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


def test_classify_s3(tmp_path):
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
