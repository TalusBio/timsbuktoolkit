import argparse
from dataclasses import dataclass
from typing import Generator

import polars as pl
import rustyms
from tqdm.auto import tqdm

from speclib_builder.base import (
    NEUTRON_MASS,
    ElutionGroup,
    MzIntPair,
    PrecursorEntry,
    SpeclibElement,
)

from .decoys import DecoyStrategy
from .isotopes import peptide_formula_dist
from .writer import SpeclibWriter


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, default="speclib.txt")
    # parser.add_argument("--input_format", type=str, default="diann_speclib")
    parser.add_argument("--add_decoys", action="store_true")
    parser.add_argument("--output", type=str, default="output.ndjson")
    return parser


@dataclass
class DiannTransitionGroup:
    transition_group_id: str
    file_name: str
    precursor_mz: float
    ion_mobility: float
    decoy: bool
    stripped_sequence: str
    modified_sequence: str
    charge: int
    ion_dict: dict[str, MzIntPair]

    def as_speclib_element(self, int_id: int):
        neutron_mass_frac = NEUTRON_MASS / self.charge
        proforma = self.modified_sequence + f"/{self.charge}"
        rs_peptide = rustyms.LinearPeptide(proforma)
        iso_dist = [0.001] + list(peptide_formula_dist(rs_peptide.formula()[0]))
        return SpeclibElement(
            precursor=PrecursorEntry(
                sequence=proforma,
                charge=self.charge,
                decoy=self.decoy,
            ),
            elution_group=ElutionGroup(
                id=int_id,
                mobility=self.ion_mobility,
                rt_seconds=0,
                precursor_mzs=[
                    self.precursor_mz + (neutron_mass_frac * isotope)
                    for isotope in [-1, 0, 1, 2]
                ],
                fragment_mzs={k: v.mz for k, v in self.ion_dict.items()},
                fragment_intensities={k: v.intensity for k, v in self.ion_dict.items()},
                precursor_intensities=iso_dist,
            ),
        )


def diann_pep_to_proforma(col: pl.Expr) -> pl.Expr:
    return col.str.replace_all("(UniMod:", "[U:", literal=True).str.replace_all(
        ")", "]", literal=True
    )


def bundle_precursors(df: pl.DataFrame) -> Generator[DiannTransitionGroup, None, None]:
    cdf = df.group_by(["transition_group_id", "FileName"]).agg(
        pl.col("PrecursorMz").unique(),
        pl.col("IonMobility").unique(),
        pl.col("decoy").unique(),
        pl.col("PeptideSequence").unique(),
        diann_pep_to_proforma(pl.col("ModifiedPeptide")).unique(),
        pl.col("PrecursorCharge").unique(),
        (
            pl.col("FragmentType")
            + pl.col("FragmentSeriesNumber").cast(pl.String)
            + "^"
            + pl.col("FragmentCharge").cast(pl.String)
        ).alias("Ions"),
        pl.col("ProductMz"),
        pl.col("LibraryIntensity"),
    )

    for x in cdf.iter_rows(named=True):
        assert len(x["PrecursorMz"]) == 1
        assert len(x["IonMobility"]) == 1
        assert len(x["decoy"]) == 1
        assert len(x["PeptideSequence"]) == 1
        assert len(x["ModifiedPeptide"]) == 1
        assert len(x["PrecursorCharge"]) == 1

        tmp = DiannTransitionGroup(
            transition_group_id=x["transition_group_id"],
            file_name=x["FileName"],
            precursor_mz=x["PrecursorMz"][0],
            ion_mobility=x["IonMobility"][0],
            decoy=x["decoy"][0],
            stripped_sequence=x["PeptideSequence"][0],
            modified_sequence=x["ModifiedPeptide"][0],
            charge=x["PrecursorCharge"][0],
            ion_dict={
                ion: MzIntPair(mz, inten)
                for ion, mz, inten in zip(
                    x["Ions"], x["ProductMz"], x["LibraryIntensity"], strict=True
                )
            },
        )
        yield tmp


def main():
    args = build_parser().parse_args()
    df = pl.read_csv(args.input, separator="\t")
    writer = SpeclibWriter(args.output)
    bundled = list(bundle_precursors(df))
    with writer as f:
        id = 0
        pbar = tqdm(bundled, desc="Converting targets")
        for x in pbar:
            elem = x.as_speclib_element(id)
            f.append(elem)
            id += 1
            if args.add_decoys:
                if elem.precursor.decoy:
                    raise ValueError("Cannot add decoys to decoys")
                decoy_elem = elem.generate_massshift_decoy(
                    id=id, decoy_strategy=DecoyStrategy.REVERSE
                )
                f.append(decoy_elem)
                id += 1

            pbar.set_postfix({"written": id})
