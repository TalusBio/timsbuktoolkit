from __future__ import annotations

import warnings

from pydantic import BaseModel, ConfigDict
from rustyms import FragmentationModel, LinearPeptide

from .decoys import DecoyStrategy, build_massshift_dict

PROTON_MASS = 1.007276466
NEUTRON_MASS = 1.008664916


class PeptideElement(BaseModel):
    peptide: str
    charge: int
    nce: float
    decoy: bool
    decoy_group: int


class MzIntPair(BaseModel):
    mz: float
    intensity: float


class PrecursorEntry(BaseModel):
    sequence: str
    charge: int
    decoy: bool
    decoy_group: int


class ElutionGroup(BaseModel):
    id: int
    mobility: float
    rt_seconds: float
    precursor_mz: float
    precursor_charge: int
    precursor_labels: list[int]
    fragment_mzs: list[float]
    fragment_labels: list[str]

    # Note that for the speclib entry, the intensities are required
    # but not for timsquery/as an output
    precursor_intensities: list[float] | None = None
    fragment_intensities: list[float] | None = None


class SpeclibElement(BaseModel):
    precursor: PrecursorEntry
    elution_group: ElutionGroup

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
        curr_fragment_labels = self.elution_group.fragment_labels
        frags = peptide.generate_theoretical_fragments(
            peptide.charge - 1, fragmentation_model
        )
        frags = [x for x in frags if x.neutral_loss is None]
        ion_mass_dict = {
            f"{f.ion}^{f.charge}": f.formula.mass() / f.charge for f in frags
        }

        keep = {}
        for k, v in zip(curr_fragment_labels, curr_fragment_mzs):
            if k in ion_mass_dict:
                keep[k] = ion_mass_dict[k]
            else:
                warnings.warn(f"Fragment {k} not found in new peptide {proforma}.")

        out = SpeclibElement(
            precursor=PrecursorEntry(
                sequence=proforma,
                charge=self.precursor.charge,
                decoy=decoy,
                decoy_group=self.precursor.decoy_group,
            ),
            elution_group=ElutionGroup(
                id=id,
                mobility=self.elution_group.mobility,
                rt_seconds=self.elution_group.rt_seconds,
                fragment_mzs=[v for v in keep.values()],
                fragment_labels=[k for k in keep.keys()],
                precursor_labels=self.elution_group.precursor_labels,
                precursor_mz=self.elution_group.precursor_mz,
                precursor_charge=self.elution_group.precursor_charge,
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

        new_fragment_mzs = []
        for k, v in zip(
            self.elution_group.fragment_labels,
            self.elution_group.fragment_mzs,
        ):
            new_fragment_mzs.append(v + massshift_dict[k])

        return SpeclibElement(
            precursor=PrecursorEntry(
                sequence=new_seq,
                charge=self.precursor.charge,
                decoy=True,
                decoy_group=self.precursor.decoy_group,
            ),
            elution_group=ElutionGroup(
                id=id,
                mobility=self.elution_group.mobility,
                rt_seconds=self.elution_group.rt_seconds,
                precursor_mz=self.elution_group.precursor_mz,
                precursor_charge=self.elution_group.precursor_charge,
                fragment_mzs=new_fragment_mzs,
                precursor_labels=self.elution_group.precursor_labels,
                fragment_labels=self.elution_group.fragment_labels,
                precursor_intensities=self.elution_group.precursor_intensities,
                fragment_intensities=self.elution_group.fragment_intensities,
            ),
        )


class EntryElements(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    peptide: LinearPeptide
    ion_dict: dict[str, MzIntPair]
    rt_seconds: float
    decoy: bool
    id: int
    decoy_group: int
