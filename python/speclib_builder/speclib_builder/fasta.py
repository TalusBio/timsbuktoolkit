from dataclasses import dataclass
from typing import Generator, Iterable

from pyteomics import fasta, parser

from .base import PeptideElement
from .decoys import DecoyStrategy, yield_with_decoys


def get_peptides(fasta_file: str) -> list[str]:
    print("Cleaving the proteins with trypsin...")
    unique_peptides = set()
    nseqs = 0
    with open(fasta_file) as file:
        for description, sequence in fasta.FASTA(file):
            nseqs += 1
            new_peptides = parser.cleave(
                sequence,
                "trypsin",
                min_length=6,
                max_length=20,
                missed_cleavages=1,
            )
            unique_peptides.update(new_peptides)

    unique_peptides = list(unique_peptides)
    print(unique_peptides[:5])
    print(f"Done, {len(unique_peptides)} sequences obtained from {nseqs} proteins!")
    return unique_peptides


def yield_with_mods(
    peptides: Iterable[tuple[str, bool, int]],
) -> Generator[tuple[str, bool, int], None, None]:
    # TODO: implement more mods ...
    for peptide in peptides:
        yield (peptide[0].replace("C", "C[UNIMOD:4]"), peptide[1], peptide[2])


def yield_with_charges(
    peptides: Iterable[tuple[str, bool, int]],
    min_charge,
    max_charge,
) -> Generator[tuple[str, bool, int, int, float], None, None]:
    for peptide in peptides:
        for charge in range(min_charge, max_charge + 1):
            yield (peptide[0], peptide[1], peptide[2], charge, 30.0)


# This is more a factory than a builder, but anyway ...
@dataclass
class PeptideBuilder:
    fasta_file: str
    min_charge: int
    max_charge: int
    decoy_strategy: DecoyStrategy

    def get_targets(self) -> list[str]:
        peps = get_peptides(fasta_file=self.fasta_file)
        return peps

    def get_modified_target_decoys(self) -> list[PeptideElement]:
        targ_use = list(
            yield_with_charges(
                yield_with_mods(
                    yield_with_decoys(self.get_targets(), self.decoy_strategy)
                ),
                self.min_charge,
                self.max_charge,
            )
        )
        return [
            PeptideElement(
                peptide=pep,
                charge=charge,
                nce=nce,
                decoy=decoy,
                decoy_group=decoy_group,
            )
            for pep, decoy, decoy_group, charge, nce in targ_use
        ]
