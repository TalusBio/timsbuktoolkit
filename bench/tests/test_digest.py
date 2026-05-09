from bench._digest import digest_protein, digest_proteins, length_filter, parse_fasta


def test_parse_fasta_basic(tmp_path):
    p = tmp_path / "p.fasta"
    p.write_text(">A\nMKL\nAAR\n>B\nPPP\n")
    out = parse_fasta(p)
    assert out == {"A": "MKLAAR", "B": "PPP"}


def test_digest_protein_no_missed_cleavage():
    """Trypsin: cleave after K/R, never before P. missed_cleavages=0."""
    # MKR|TPK|GCD - cleave after K, R, K (KP rule: K-P would not cleave)
    pep = digest_protein("MKRTPKGCD", missed_cleavages=0)
    assert pep == ["MK", "R", "TPK", "GCD"]


def test_digest_protein_kp_rp_rule():
    """K-P and R-P bonds are NOT cleaved."""
    # AAKPBBR - the K-P bond stays, R at end cleaves
    pep = digest_protein("AAKPBBR", missed_cleavages=0)
    assert pep == ["AAKPBBR"]


def test_digest_protein_one_missed_cleavage():
    pep = digest_protein("MKRTPKGCD", missed_cleavages=1)
    # Base segments: ["MK", "R", "TPK", "GCD"]
    # Plus contiguous merges of length 2: ["MKR", "RTPK", "TPKGCD"]
    assert "MK" in pep and "R" in pep and "TPK" in pep and "GCD" in pep
    assert "MKR" in pep and "RTPK" in pep and "TPKGCD" in pep


def test_digest_proteins_dedupes():
    # Two proteins sharing a peptide
    pep = digest_proteins({"A": "MKR", "B": "MKR"}, missed_cleavages=0)
    # set ordering doesn't matter; check membership
    assert "MK" in pep
    assert "R" in pep


def test_length_filter():
    peps = {"A", "AA", "AAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "AAAAAA"}
    out = length_filter(peps, min_len=7, max_len=30)
    # 7-mer kept; 6-mer dropped; 31-mer dropped; 1-mer/2-mer dropped
    assert out == {"AAAAAAA"}
