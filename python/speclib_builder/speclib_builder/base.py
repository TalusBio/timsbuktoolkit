from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass

from rustyms import FragmentationModel, LinearPeptide

from .decoys import DecoyStrategy, build_massshift_dict

PROTON_MASS = 1.007276466
NEUTRON_MASS = 1.008664916


@dataclass(slots=True, frozen=True)
class PeptideElement:
    peptide: str
    charge: int
    nce: float
    decoy: bool


@dataclass(slots=True, frozen=True)
class MzIntPair:
    mz: float
    intensity: float


@dataclass(slots=True, frozen=True)
class PrecursorEntry:
    sequence: str
    charge: int
    decoy: bool


@dataclass
class ElutionGroup:
    id: int
    mobility: float
    rt_seconds: float
    precursor_mzs: list[float]
    fragment_mzs: dict[str, float]
    precursor_intensities: list[float]
    fragment_intensities: dict[str, float]


@dataclass
class SpeclibElement:
    precursor: PrecursorEntry
    elution_group: ElutionGroup

    def as_rt_entry(self) -> dict:
        return {
            "sequence": self.precursor.sequence,
            "charge": self.precursor.charge,
            "elution_group": {
                "id": self.elution_group.id,
                "mobility": self.elution_group.mobility,
                "rt_seconds": self.elution_group.rt_seconds,
                "precursor_mzs": self.elution_group.precursor_mzs,
                "fragment_mzs": self.elution_group.fragment_mzs,
            },
            "expected_intensities": {
                "fragment_intensities": self.elution_group.fragment_intensities,
                "precursor_intensities": self.elution_group.precursor_intensities,
            },
        }

    def as_dict(self) -> dict:
        return asdict(self)

    def swap_peptide(
        self,
        proforma: str,
        decoy: bool,
        id: int,
        fragmentation_model: FragmentationModel = FragmentationModel.CidHcd,
    ) -> SpeclibElement:
        """Swap the peptide for an elution group.

        This function also swaps the masses of the ions that are contained in the
        current elution group.

        Notably, Right now it does not change the precursor mass.... or anything
        other than the fragment m/z values.
        """
        peptide = LinearPeptide(proforma)
        curr_fragment_mzs = self.elution_group.fragment_mzs
        frags = peptide.generate_theoretical_fragments(
            peptide.charge - 1, fragmentation_model
        )
        frags = [x for x in frags if x.neutral_loss is None]
        ion_mass_dict = {
            f"{f.ion}^{f.charge}": f.formula.mass() / f.charge for f in frags
        }

        keep = {}
        for k, v in curr_fragment_mzs.items():
            if k in ion_mass_dict:
                keep[k] = ion_mass_dict[k]
            else:
                warnings.warn(f"Fragment {k} not found in new peptide {proforma}.")

        out = SpeclibElement(
            precursor=PrecursorEntry(
                sequence=proforma,
                charge=self.precursor.charge,
                decoy=decoy,
            ),
            elution_group=ElutionGroup(
                id=id,
                mobility=self.elution_group.mobility,
                rt_seconds=self.elution_group.rt_seconds,
                precursor_mzs=self.elution_group.precursor_mzs,
                fragment_mzs=keep,
                precursor_intensities=self.elution_group.precursor_intensities,
                fragment_intensities=self.elution_group.fragment_intensities,
            ),
        )

        return out

    def generate_massshift_decoy(
        self,
        id: int,
        decoy_strategy: DecoyStrategy,
        fragmentation_model: FragmentationModel = FragmentationModel.CidHcd,
    ) -> SpeclibElement:
        massshift_dict, new_seq = build_massshift_dict(
            self.precursor.sequence,
            decoy_strategy,
            max_charge=self.precursor.charge,
            fragmentation_model=fragmentation_model,
        )

        new_fragment_mzs = {}
        for k, v in self.elution_group.fragment_mzs.items():
            new_fragment_mzs[k] = v + massshift_dict[k]

        return SpeclibElement(
            precursor=PrecursorEntry(
                sequence=new_seq,
                charge=self.precursor.charge,
                decoy=True,
            ),
            elution_group=ElutionGroup(
                id=id,
                mobility=self.elution_group.mobility,
                rt_seconds=self.elution_group.rt_seconds,
                precursor_mzs=self.elution_group.precursor_mzs,
                fragment_mzs=new_fragment_mzs,
                precursor_intensities=self.elution_group.precursor_intensities,
                fragment_intensities=self.elution_group.fragment_intensities,
            ),
        )


@dataclass
class EntryElements:
    peptide: LinearPeptide
    ion_dict: dict[str, MzIntPair]
    rt_seconds: float
    decoy: bool
    id: int
