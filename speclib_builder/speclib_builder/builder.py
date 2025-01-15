from dataclasses import dataclass
from typing import Generator

import rustyms
from rustyms import LinearPeptide

from .base import (
    NEUTRON_MASS,
    PROTON_MASS,
    ElutionGroup,
    EntryElements,
    MzIntPair,
    PeptideElement,
    PrecursorEntry,
    SpeclibElement,
)
from .ims import supersimpleprediction
from .isotopes import peptide_formula_dist


@dataclass
class PeptideAnnotator:
    def model(self, PeptideElement) -> EntryElements:
        raise NotImplementedError

    def batched_model(
        self, elements: list[PeptideElement]
    ) -> Generator[EntryElements, None, None]:
        for elem in elements:
            yield self.model(elem)


class DummyAnnotator(PeptideAnnotator):
    fragmentation_model: rustyms.FragmentationModel = rustyms.FragmentationModel.CidHcd

    def model(self, elem: PeptideElement) -> EntryElements:
        seq = elem.peptide
        charge = elem.charge
        decoy = elem.decoy

        pep = rustyms.LinearPeptide(f"{seq}/{charge}")
        frags = pep.generate_theoretical_fragments(
            pep.charge - 1, self.fragmentation_model
        )
        frags = [x for x in frags if x.neutral_loss is None]
        frags = [x for x in frags if x.charge <= (pep.charge - 1)]
        ion_dict = {
            f"{f.ion}^{f.charge}": MzIntPair(f.formula.mass() / f.charge, 1.0)
            for f in frags
        }
        return EntryElements(pep, ion_dict, decoy, 0)


@dataclass
class EntryBuilder:
    max_ions_keep: int = 20
    min_mz: float = 400
    max_mz: float = 2000
    min_ion_mz: float = 250
    max_ion_mz: float = 2000
    min_ions: int = 3

    def build_entry(self, elem: EntryElements) -> dict | None:
        return self.as_entry(elem.peptide, elem.ion_dict, elem.decoy, elem.id)

    def as_entry(
        self,
        peptide: LinearPeptide,
        ion_dict: dict[str, MzIntPair],
        decoy: bool,
        id: int,
    ) -> SpeclibElement | None:
        pep_formula = peptide.formula()[0]
        isotope_dist = peptide_formula_dist(pep_formula)
        precursor_mz = (
            pep_formula.mass() + (PROTON_MASS * peptide.charge)
        ) / peptide.charge
        if precursor_mz < self.min_mz or precursor_mz > self.max_mz:
            return None

        ims = float(supersimpleprediction(precursor_mz, peptide.charge))

        neutron_fraction = NEUTRON_MASS / peptide.charge
        precursor_mzs = [
            float(precursor_mz + (neutron_fraction * isotope))
            for isotope in [-1, 0, 1, 2]
        ]

        intensities = [v.intensity for v in ion_dict.values()]
        intensities.sort(reverse=True)
        max_intensity = intensities[0]
        min_inten_keep = max_intensity * 0.02
        intensities = [x for x in intensities if x > min_inten_keep]
        if len(intensities) > self.max_ions_keep:
            intensities = intensities[: self.max_ions_keep]

        min_inten_keep = intensities[-1]

        ion_dict = {
            k: v
            for k, v in ion_dict.items()
            if (v.mz > self.min_ion_mz)
            and (v.intensity >= min_inten_keep)
            and (v.mz < self.max_ion_mz)
        }

        ion_mzs = {k: v.mz for k, v in ion_dict.items()}
        ion_intensities = {k: v.intensity for k, v in ion_dict.items()}
        if len(ion_mzs) < self.min_ions:
            return None

        precursor_entry = {
            "sequence": str(peptide),
            "charge": peptide.charge,
            "decoy": decoy,
        }
        elution_group = {
            "id": id,
            "mobility": ims,
            "rt_seconds": 0.0,
            "precursor_mzs": precursor_mzs,
            "fragment_mzs": ion_mzs,
        }
        expected_intensities = {
            "precursor_intensities": [0.001] + list(isotope_dist),
            "fragment_intensities": ion_intensities,
        }

        entry = SpeclibElement(
            PrecursorEntry(**precursor_entry),
            ElutionGroup(**elution_group, **expected_intensities),
        )
        return entry
