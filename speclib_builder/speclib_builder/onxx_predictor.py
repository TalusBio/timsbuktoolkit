from __future__ import annotations

try:
    from elfragmentadonnx.model import OnnxPeptideTransformer
except ImportError:
    print("Could not import the OnnxPeptideTransformer, ML prediction will not work")
    OnnxPeptideTransformer = None

import re
from dataclasses import dataclass
from typing import Generator, Iterable

from .base import MzIntPair, PeptideElement
from .builder import EntryElements, PeptideAnnotator

MOD_REGEX = re.compile(r"(\[[A-Z:0-9]*?\])")


@dataclass(slots=True)
class OnnxPeptideTransformerAnnotator(PeptideAnnotator):
    inference_model: OnnxPeptideTransformer
    min_ordinal: int
    max_ordinal: int
    min_intensity: float
    num_yielded: int = 0
    max_tokens: int = 30

    def model(self, elem: PeptideElement) -> EntryElements:
        return list(self.batched_model([elem]))[0]

    def filter_tokenizable(
        self, elems: Iterable[PeptideElement]
    ) -> list[PeptideElement]:
        out = []
        for elem in elems:
            if (token_count := self.count_tokens(elem.peptide)) > self.max_tokens:
                print(f"Skipping {elem.peptide} because it has {token_count} tokens")
            else:
                out.append(elem)
        return out

    def yielding_adapter(
        self, elems: Iterable[PeptideElement]
    ) -> Generator[tuple[str, int, float] | None, None, None]:
        for elem in elems:
            yield elem.peptide, elem.charge, elem.nce

    @staticmethod
    def count_tokens(peptide: str) -> int:
        # Pretty heuristic for now ... seems to work fine and fast enough
        num_mods = peptide.count("[")
        stripped_peptide = MOD_REGEX.sub("", peptide).strip("-")
        return num_mods + len(stripped_peptide) + 2

    def batched_model(
        self, elements: list[PeptideElement]
    ) -> Generator[EntryElements, None, None]:
        elems = self.filter_tokenizable(elements)

        for pe, elem in zip(
            elems,
            self.inference_model.predict_batched_annotated(
                self.yielding_adapter(elems),
                min_intensity=self.min_intensity,
                min_ordinal=self.min_ordinal,
                max_ordinal=self.max_ordinal,
            ),
            strict=True,
        ):
            if elem is None:
                continue
            pep, ion_dict = elem
            ion_dict = {k: MzIntPair(v[0], v[1]) for k, v in ion_dict.items()}
            yield EntryElements(
                pep,
                ion_dict=ion_dict,
                decoy=pe.decoy,
                id=self.num_yielded,
            )
            self.num_yielded += 1

    @staticmethod
    def get_default() -> "OnnxPeptideTransformerAnnotator":
        if OnnxPeptideTransformer is None:
            raise ImportError(
                "Could not import the OnnxPeptideTransformer, ML prediction will not work"
            )
        return OnnxPeptideTransformerAnnotator(
            inference_model=OnnxPeptideTransformer.default_model(),
            min_ordinal=2,
            max_ordinal=1000,
            min_intensity=0.001,
        )
