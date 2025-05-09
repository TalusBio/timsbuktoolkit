import enum
import json

from rich.pretty import pprint
from tqdm.auto import tqdm

from .builder import DummyAnnotator, EntryBuilder
from .decoys import DecoyStrategy
from .fasta import PeptideBuilder
from .onxx_predictor import OnnxPeptideTransformerAnnotator


class ModelType(enum.Enum):
    ONNX = "onnx"
    DUMMY = "dummy"


def build_parser():
    import argparse

    # parser = argparse.ArgumentParser()
    # Same as above but shows default values and such ... also prettyer using rich
    parser = argparse.ArgumentParser(
        description="Build a speclib from a fasta file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--fasta_file",
        type=str,
        default="/Users/sebastianpaez/fasta/20231030_UP000005640_9606.fasta",
        help="Fasta file to use",
    )
    parser.add_argument(
        "--decoy_strategy",
        type=str,
        default="REVERSE",
        choices=["REVERSE", "MUTATE", "EDGE_MUTATE"],
        help="Decoy strategy to use",
    )
    parser.add_argument(
        "--max_ions",
        type=int,
        default=10,
        help="Maximum number of ions to keep per precursor",
    )
    parser.add_argument(
        "--outfile", type=str, default="FUUUUU.ndjson", help="Output file"
    )
    parser.add_argument(
        "--model",
        type=str,
        default="onnx",
        help="Model to use, the dummy model just sets all intensities to 1",
        choices=["onnx", "dummy"],
    )
    return parser


def main():
    args = build_parser().parse_args()

    if args.decoy_strategy == "REVERSE":
        decoy_strategy = DecoyStrategy.REVERSE
    elif args.decoy_strategy == "MUTATE":
        decoy_strategy = DecoyStrategy.MUTATE
    elif args.decoy_strategy == "EDGE_MUTATE":
        decoy_strategy = DecoyStrategy.EDGE_MUTATE
    else:
        raise NotImplementedError(f"Unknown decoy strategy {args.decoy_strategy}")

    if args.model == "onnx":
        annotator = OnnxPeptideTransformerAnnotator.get_default()
    elif args.model == "dummy":
        annotator = DummyAnnotator()
    else:
        raise NotImplementedError(f"Unknown model {args.model}")

    fasta_file = args.fasta_file
    outfile = args.outfile
    max_keep = args.max_ions

    _main(
        fasta_file=fasta_file,
        outfile=outfile,
        max_keep=max_keep,
        decoy_strategy=decoy_strategy,
        annotator=annotator,
    )


def _main(
    *,
    fasta_file: str,
    outfile: str,
    max_keep: int,
    decoy_strategy: DecoyStrategy,
    annotator: OnnxPeptideTransformerAnnotator | DummyAnnotator,
):
    pretty_outfile = f"{outfile}.pretty.json"

    peptide_builder = PeptideBuilder(
        fasta_file=fasta_file,
        min_charge=2,
        max_charge=4,
        decoy_strategy=decoy_strategy,
    )

    # # outfile = "FUUUUU_small.ndjson"
    # outfile = "FUUUUU.ndjson"
    # fasta_file = "/Users/sebastianpaez/fasta/20231030_UP000005640_9606.fasta"
    # # fasta_file = "/Users/sebastianpaez/git/timsseek/data/HeLa_cannonical_proteins.fasta"
    entry_builder = EntryBuilder(
        min_ions=3,
        max_ions_keep=max_keep,
        min_mz=400,
        max_mz=2000,
        min_ion_mz=250,
        max_ion_mz=2000,
    )

    pretty_outs = []
    is_first_n = 10
    id = 0

    pprint(f"Writing output to file: {outfile}")
    with open(outfile, "w") as file:
        targ_use = peptide_builder.get_modified_targets()
        for x in tqdm(
            annotator.batched_model(targ_use),
            desc="Targets",
            total=len(targ_use),
        ):
            id += 1
            elem = entry_builder.as_entry(
                peptide=x.peptide,
                decoy=False,
                id=id,
                ion_dict=x.ion_dict,
                rt_seconds=x.rt_seconds,
            )
            if elem is None:
                continue
            if is_first_n > 0:
                pprint(elem)
                pretty_outs.append(elem)
                is_first_n -= 1

            file.write(elem.model_dump_json() + "\n")
            file.flush()

        is_first_n = 10

        decoys = peptide_builder.get_modified_decoys()
        for x in tqdm(
            annotator.batched_model(decoys),
            desc="Decoys",
            total=len(decoys),
        ):
            id += 1
            elem = entry_builder.as_entry(
                peptide=x.peptide,
                ion_dict=x.ion_dict,
                decoy=True,
                id=id,
                rt_seconds=x.rt_seconds,
            )
            if elem is None:
                continue
            if is_first_n > 0:
                pprint(elem)
                pretty_outs.append(elem)
                is_first_n -= 1

            file.write(elem.model_dump_json() + "\n")
            file.flush()

    pprint(f"Writing pretty output to file: {pretty_outfile}")
    with open(pretty_outfile, "w") as file:
        file.write(json.dumps([x.model_dump() for x in pretty_outs], indent=4))
        file.flush()


if __name__ == "__main__":
    # Run the test
    main()
