import responses

from bench._uniprot import fetch_accession_batch, fetch_proteome


@responses.activate
def test_fetch_proteome():
    responses.add(
        responses.GET,
        "https://rest.uniprot.org/uniprotkb/stream",
        body=">sp|P12345|TEST_HUMAN test\nMKLAA\n",
        status=200,
        content_type="text/plain",
    )
    fasta = fetch_proteome("UP000005640")
    assert fasta.startswith(">sp|P12345")
    assert "MKLAA" in fasta


@responses.activate
def test_fetch_accession_batch():
    responses.add(
        responses.GET,
        "https://rest.uniprot.org/uniprotkb/stream",
        body=">sp|P12345|A\nMK\n>sp|Q67890|B\nLL\n",
        status=200,
        content_type="text/plain",
    )
    fasta = fetch_accession_batch(["P12345", "Q67890"])
    assert "P12345" in fasta and "Q67890" in fasta
    # Verify the request used the OR'd accession query
    assert len(responses.calls) == 1
    qs = responses.calls[0].request.url or ""
    assert "accession%3AP12345" in qs and "accession%3AQ67890" in qs
    assert "format=fasta" in qs


@responses.activate
def test_fetch_accession_batch_chunks_long_lists():
    """Uniprot stream URLs have practical length limits; verify we batch."""
    # Stub one response per chunk regardless of how many calls happen
    responses.add(
        responses.GET,
        "https://rest.uniprot.org/uniprotkb/stream",
        body=">sp|X|x\nM\n",
        status=200,
    )
    accs = [f"P{n:05d}" for n in range(250)]  # 250 → at least 2 chunks at 100/chunk
    fetch_accession_batch(accs)
    assert len(responses.calls) == 3  # ceil(250/100)
