"""Shuffled-entrapment peptide generation (Algorithm 1 of Noble et al, FDRBench paper).

`shuffle_keeping_c_terminal` permutes the interior residues of a peptide and
keeps the C-terminal residue fixed. `generate_shuffled_entrapment` produces
the (target, shuffle) pairs per Algorithm 1, with deduplication and
max-attempt fallback that drops targets unable to produce r distinct shuffles.
"""

from __future__ import annotations

import random


def shuffle_keeping_c_terminal(peptide: str, seed: int | None = None) -> str:
    """Shuffle the interior residues, keep the C-terminal residue fixed.

    For peptides of length <= 2 returns the input unchanged (no interior).
    """
    if len(peptide) <= 2:
        return peptide
    rng = random.Random(seed)
    interior = list(peptide[:-1])
    rng.shuffle(interior)
    return "".join(interior) + peptide[-1]


def generate_shuffled_entrapment(
    targets: set[str],
    r: int = 1,
    seed: int = 42,
) -> list[tuple[str, str]]:
    """Algorithm 1: for each target, generate r distinct shuffles.

    Uses `max_attempts = 20 + r` attempts per target. Targets that cannot
    produce r unique shuffles (e.g., homopolymers) are dropped — i.e., they
    contribute zero pairs to the output.

    Returns a list of (target_peptide, shuffle_peptide) pairs. Pairs are
    sorted by target for determinism. RNG is seeded once at the start so
    re-runs with the same `seed` produce the same pairs.
    """
    rng = random.Random(seed)
    max_attempts = 20 + r
    out: list[tuple[str, str]] = []
    for p_target in sorted(targets):  # determinism via sorted iteration
        shuffles: set[str] = set()
        for _ in range(max_attempts):
            if len(shuffles) >= r:
                break
            # New shuffle drawn from the seeded shared RNG
            interior = list(p_target[:-1])
            rng.shuffle(interior)
            cand = "".join(interior) + p_target[-1]
            if cand != p_target and cand not in shuffles and cand not in targets:
                shuffles.add(cand)
        if len(shuffles) >= r:
            for s in sorted(shuffles):
                out.append((p_target, s))
        # else: drop p_target entirely (no entries appended)
    return out
