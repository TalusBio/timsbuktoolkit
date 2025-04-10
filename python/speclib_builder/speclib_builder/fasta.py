from dataclasses import dataclass
from typing import Generator

from pyteomics import fasta, parser

from .base import PeptideElement
from .decoys import DecoyStrategy, yield_as_decoys


def get_peptides(fasta_file: str) -> list[str]:
    print("Cleaving the proteins with trypsin...")
    unique_peptides = set()
    with open(fasta_file) as file:
        for description, sequence in fasta.FASTA(file):
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
    print("Done, {0} sequences obtained!".format(len(unique_peptides)))
    return unique_peptides


def yield_with_mods(peptides: list[str]) -> Generator[str, None, None]:
    for peptide in peptides:
        yield peptide.replace("C", "C[UNIMOD:4]")


def yield_with_charges(
    peptides,
    min_charge,
    max_charge,
) -> Generator[tuple[str, int, float], None, None]:
    for peptide in peptides:
        for charge in range(min_charge, max_charge + 1):
            yield (peptide, charge, 30.0)


@dataclass
class PeptideBuilder:
    fasta_file: str
    min_charge: int
    max_charge: int
    decoy_strategy: DecoyStrategy

    def get_targets(self) -> list[str]:
        peps = get_peptides(fasta_file=self.fasta_file)
        return peps

    def get_modified_targets(self) -> list[PeptideElement]:
        targ_use = list(
            yield_with_charges(
                yield_with_mods(self.get_targets()), self.min_charge, self.max_charge
            )
        )
        return [
            PeptideElement(pep, charge, nce, False) for pep, charge, nce in targ_use
        ]

    def get_decoys(self) -> list[str]:
        peps = get_peptides(fasta_file=self.fasta_file)
        decoys = list(yield_as_decoys(peps, self.decoy_strategy))
        return decoys

    def get_modified_decoys(self) -> list[PeptideElement]:
        decoys = self.get_decoys()
        decoys_use = list(
            yield_with_charges(
                yield_with_mods(decoys), self.min_charge, self.max_charge
            )
        )
        return [
            PeptideElement(pep, charge, nce, True) for pep, charge, nce in decoys_use
        ]
