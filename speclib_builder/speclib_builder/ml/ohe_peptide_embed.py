import math

import numpy as np
from rustyms import LinearPeptide

# Model pretty much translating the implementation from sage

VALID_AA: str = "ACDEFGHIKLMNPQRSTVWY"
FEATURES: int = len(VALID_AA) * 3 + 2
N_TERMINAL: int = len(VALID_AA)
C_TERMINAL: int = len(VALID_AA) * 2
PEPTIDE_LEN: int = FEATURES - 3
PEPTIDE_MASS_LN: int = FEATURES - 2
INTERCEPT: int = FEATURES - 1


class PeptideOHEEmbedder:
    def embed(self, peptide: LinearPeptide) -> np.ndarray:
        stripped_seq = peptide.stripped_sequence
        mass = peptide.formula()[0].monoisotopic_mass()
        return self.embed_stripped_sequence(stripped_seq, mass)

    def embed_stripped_sequence(self, stripped_seq: str, mass: float) -> np.ndarray:
        embedding = [0.0] * FEATURES
        for aa_idx, residue in enumerate(stripped_seq):
            try:
                idx = VALID_AA.index(residue)
            except ValueError:
                raise ValueError(f"Invalid residue {residue} at position {aa_idx}")
            embedding[idx] += 1.0
            # Embed N- and C-terminal AA's (2 on each end, excluding K/R)
            if aa_idx <= 1:
                embedding[N_TERMINAL + idx] += 1.0
            elif aa_idx >= len(stripped_seq) - 2:
                embedding[C_TERMINAL + idx] += 1.0
        embedding[PEPTIDE_LEN] = len(stripped_seq)
        embedding[PEPTIDE_MASS_LN] = math.log1p(mass)
        embedding[INTERCEPT] = 1.0

        return np.array(embedding)
