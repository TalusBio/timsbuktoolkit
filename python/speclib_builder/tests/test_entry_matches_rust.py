import json
import subprocess
from pathlib import Path

import pytest
from speclib_builder.base import ElutionGroup, SpeclibElement
from speclib_builder.decoys import DecoyStrategy
from speclib_builder.fasta_cli import (
    DummyAnnotator,
    OnnxPeptideTransformerAnnotator,
    _main,
)


@pytest.fixture
def rust_json():
    # Calls the rust side to generate a json serialized version
    # of the library file, we will use that to compare
    # against the python generated version (so if we can serialize-deserialize)
    # and get the same result, we are good.
    args = [
        "cargo",
        "run",
        "--bin",
        "timsseek_sample_speclib",
        "sample",
    ]

    outs = subprocess.run(args, check=True, capture_output=True)
    _json_sample = json.loads(outs.stdout)
    # files = {k.stem: k for k in list(Path(tmpdir).rglob("*.json"))}
    yield outs.stdout


def test_speclib_ser(rust_json):
    speclib_list = json.loads(rust_json)
    # Will raise a validation error if the
    # python model does not match the example from the rust model.
    _data = ElutionGroup(**speclib_list[0]["elution_group"])
    _data = SpeclibElement(**speclib_list[0])
    pass


@pytest.mark.parametrize(
    "decoy_strategy",
    [
        DecoyStrategy.REVERSE,
        DecoyStrategy.MUTATE,
        DecoyStrategy.EDGE_MUTATE,
    ],
)
@pytest.mark.parametrize(
    "annotator",
    [
        DummyAnnotator,
        OnnxPeptideTransformerAnnotator,
    ],
)
def test_speclib_deser(tmpdir, annotator, decoy_strategy):
    fasta_path = Path(tmpdir) / "test.fasta"
    outfile = Path(tmpdir) / "test.ndjson"
    fasta_contents = [
        ">sp|P62805|H4_HUMAN Histone H4 OS=Homo sapiens OX=9606 GN=H4C1 PE=1 SV=2",
        "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLK",
        "VFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG",
    ]
    with open(fasta_path, "w") as f:
        f.writelines([f"{x}\n" for x in fasta_contents])

    if hasattr(annotator, "get_default"):
        ann = annotator.get_default()
    else:
        ann = annotator()

    _main(
        fasta_file=str(fasta_path),
        outfile=str(outfile),
        max_keep=20,
        decoy_strategy=decoy_strategy,
        annotator=ann,
    )

    with open(outfile) as f:
        lines = f.readlines()
        assert lines, "No speclib was written."
        for line in lines:
            line_con = json.loads(line)
            # Make sure it passes validation
            _ = SpeclibElement(**line_con)
