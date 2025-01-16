import argparse
import itertools
import multiprocessing
import os
from dataclasses import dataclass
from typing import Generator

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
    parser.add_argument("--num_workers", type=int, default=multiprocessing.cpu_count())
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


def diann_pep_to_proforma(pep: str) -> str:
    return pep.replace("(UniMod:", "[U:").replace(")", "]")


class DiannTransitionGroupBuilder:
    def __init__(self):
        self.transition_group_id = None
        self.file_name = None
        self.precursor_mz = None
        self.ion_mobility = None
        self.decoy = None
        self.stripped_sequence = None
        self.modified_sequence = None
        self.charge = None
        self.ions = []
        self.mzs = []
        self.intensities = []

    @classmethod
    def new(cls, values: dict) -> "DiannTransitionGroupBuilder":
        builder = cls()
        builder.transition_group_id = values["transition_group_id"]
        builder.file_name = values["FileName"]
        builder.precursor_mz = float(values["PrecursorMz"])
        builder.ion_mobility = float(values["IonMobility"])
        builder.decoy = bool(int(values["decoy"]))
        builder.stripped_sequence = values["PeptideSequence"]
        builder.modified_sequence = diann_pep_to_proforma(values["ModifiedPeptide"])
        builder.charge = int(values["PrecursorCharge"])
        builder.add_ion(values)
        return builder

    def add_ion(self, values: dict) -> None:
        ion = f"{values['FragmentType']}{values['FragmentSeriesNumber']}^{values['FragmentCharge']}"
        self.ions.append(ion)
        self.mzs.append(float(values["ProductMz"]))
        self.intensities.append(float(values["LibraryIntensity"]))

    def build(self) -> DiannTransitionGroup:
        return DiannTransitionGroup(
            transition_group_id=self.transition_group_id,
            file_name=self.file_name,
            precursor_mz=self.precursor_mz,
            ion_mobility=self.ion_mobility,
            decoy=self.decoy,
            stripped_sequence=self.stripped_sequence,
            modified_sequence=self.modified_sequence,
            charge=self.charge,
            ion_dict={
                ion: MzIntPair(mz, inten)
                for ion, mz, inten in zip(
                    self.ions,
                    self.mzs,
                    self.intensities,
                    strict=True,
                )
            },
        )


def bundle_precursors(file) -> Generator[DiannTransitionGroup, None, None]:
    last_id = None
    last_filename = None
    current_builder = None

    with open(file) as f:
        # Get column positions from header
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        for line in f:
            # Split the line and get key fields
            values = line.strip().split("\t")
            row_dict = {name: values[i] for name, i in col_idx.items()}
            curr_id = row_dict["transition_group_id"]
            curr_filename = row_dict["FileName"]

            # If we've hit a new group, yield the previous one
            if last_id is not None and (
                curr_id != last_id or curr_filename != last_filename
            ):
                yield current_builder.build()
                current_builder = DiannTransitionGroupBuilder.new(row_dict)
            elif current_builder is None:
                current_builder = DiannTransitionGroupBuilder.new(row_dict)
            else:
                current_builder.add_ion(row_dict)

            last_id = curr_id
            last_filename = curr_filename

        # Don't forget to yield the last group
        if current_builder is not None:
            yield current_builder.build()


def _aprox_count_lines(file):
    # Estimate number of lines by reading the 2-200th line
    # and counting the bytes in them. THEN divide the total
    # file size by the average line size.
    byte_sizes = []
    tot_file_size = os.path.getsize(file)

    with open(file) as f:
        last = None
        count = 0
        local_char_size = 0
        for line in f:
            local_char_size += len(line)
            curr = tuple(line.split("\t")[:2])
            if curr != last:
                count += 1
                last = curr
                byte_sizes.append(local_char_size)
                local_char_size = 0

            if count > 20000:
                break

    if count < 20000:
        return count

    byte_sizes = byte_sizes[1:]
    mean_byte_size = sum(byte_sizes) / len(byte_sizes)
    count = int(tot_file_size / mean_byte_size)

    return count


def _count_lines(file):
    with open(file) as f:
        last = None
        count = 0
        for line in tqdm(f, desc="Counting lines"):
            curr = tuple(line.split("\t")[:2])
            if curr != last:
                count += 1
                last = curr

        return count


def process_element(args):
    x, start_id, add_decoys = args
    elem = x.as_speclib_element(start_id)

    decoy_elem = None
    if add_decoys:
        if elem.precursor.decoy:
            raise ValueError("Cannot add decoys to decoys")
        decoy_elem = elem.generate_massshift_decoy(
            id=start_id + 1,
            decoy_strategy=DecoyStrategy.REVERSE,
        )

    return (elem, decoy_elem)


def _parallel_process(args, writer, nlines, max_workers):
    with writer as f:
        precursors = list()

        with multiprocessing.Pool(max_workers) as pool:
            # Create arguments for each precursor
            precs = bundle_precursors(args.input)
            ranges = _infinite_range(0, 2 if args.add_decoys else 1)
            process_args = zip(
                precs,
                ranges,
                itertools.repeat(args.add_decoys),
            )
            inner_iter = pool.imap_unordered(process_element, process_args)

            # Process and write results as they come in
            pbar = tqdm(
                inner_iter,
                total=len(precursors),
                desc="Converting targets",
            )
            # The chaining is required bc the "total" in tqdm can be smaller
            # than the number of precursors, which makes it truncate prematurely.
            # The chain makes sure we process the rest.
            for result_group in itertools.chain(pbar, tqdm(inner_iter)):
                for elem in result_group:
                    if elem is not None:
                        f.append(elem)


def _infinite_range(start, step):
    while True:
        yield start
        start += step


def main():
    args = build_parser().parse_args()
    writer = SpeclibWriter(args.output)
    aprox_lines = _aprox_count_lines(args.input)
    print(f"Approximate number of entries: {aprox_lines}")
    if args.num_workers >= 2:
        print(f"Using {args.num_workers} workers")
        _parallel_process(args, writer, aprox_lines, args.num_workers)
    else:
        print(
            "Using serial processing (pass "
            "--num_workers > 2 to use parallel processing)"
        )
        _serial_process(args, writer, aprox_lines)


def _serial_process(args, writer, nlines):
    with writer as f:
        id = 0
        precs = bundle_precursors(args.input)
        pbar = tqdm(precs, desc="Converting targets", total=nlines)
        # The chaining is required bc the "total" in tqdm can be smaller
        # than the number of precursors, which makes it truncate prematurely.
        # The chain makes sure we process the rest.
        for x in itertools.chain(pbar, tqdm(precs)):
            elem = x.as_speclib_element(id)
            f.append(elem)
            id += 1
            if args.add_decoys:
                if elem.precursor.decoy:
                    raise ValueError("Cannot add decoys to decoys")
                decoy_elem = elem.generate_massshift_decoy(
                    id=id,
                    decoy_strategy=DecoyStrategy.REVERSE,
                )
                f.append(decoy_elem)
                id += 1

            pbar.set_postfix({"written": id})


if __name__ == "__main__":
    main()
