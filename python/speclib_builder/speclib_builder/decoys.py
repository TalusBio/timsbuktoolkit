import enum
from typing import Generator, Literal

import rustyms
from pyteomics.proforma import ProForma
from pyteomics.proforma import parse as proforma_parse

MUTATE_DICT = {}
for old, new in zip("GAVLIFMPWSCTYHKRQEND", "LLLVVLLLLTSSSSLLNDQE"):
    MUTATE_DICT[old] = new


class DecoyStrategy(enum.Enum):
    REVERSE = "REVERSE"
    MUTATE = "MUTATE"
    EDGE_MUTATE = "EDGE_MUTATE"


MIN_ORD = ord("A")
MAX_ORD = ord("Z")


def as_decoy(x: str, decoy_strategy: DecoyStrategy) -> str:
    for c in x:
        if ord(c) < MIN_ORD or ord(c) > MAX_ORD:
            return as_decoy_proforma(x, decoy_strategy)

    if decoy_strategy == DecoyStrategy.REVERSE:
        out = x[0] + x[-2:0:-1] + x[-1]
    elif decoy_strategy == DecoyStrategy.MUTATE:
        out = "".join([MUTATE_DICT[w] for w in x])
    elif decoy_strategy == DecoyStrategy.EDGE_MUTATE:
        out = "".join([x[0] + MUTATE_DICT[x[1]] + x[2:-2] + MUTATE_DICT[x[-2]] + x[-1]])
    else:
        raise NotImplementedError(f"Unknown decoy strategy {decoy_strategy}")
    return out


def as_decoy_proforma(proforma: str, decoy_strategy: DecoyStrategy) -> str:
    parsed = proforma_parse(proforma)
    inner_seq = parsed[0]

    if decoy_strategy == DecoyStrategy.REVERSE:
        inner_seq = [inner_seq[0]] + inner_seq[-2:0:-1] + [inner_seq[-1]]
    elif decoy_strategy == DecoyStrategy.MUTATE:
        out = []
        for x in inner_seq:
            new_aa = MUTATE_DICT[x[0]]
            out.append((new_aa, x[1] if new_aa == x[0] else None))
        inner_seq = out
    elif decoy_strategy == DecoyStrategy.EDGE_MUTATE:
        new_aa_1 = MUTATE_DICT[inner_seq[1][0]]
        new_aa_neg2 = MUTATE_DICT[inner_seq[-2][0]]
        inner_seq[1] = (
            new_aa_1,
            inner_seq[1][1] if new_aa_1 == inner_seq[1][0] else None,
        )
        inner_seq[-2] = (
            new_aa_neg2,
            inner_seq[-2][1] if new_aa_neg2 == inner_seq[-2][0] else None,
        )

    else:
        raise NotImplementedError(f"Unknown decoy strategy {decoy_strategy}")

    tmp = ProForma(inner_seq, parsed[1])
    return str(tmp).replace("UNIMOD:", "U:")


def yield_with_decoys(
    peptides: list[str], decoy_strategy: DecoyStrategy
) -> Generator[tuple[str, bool, int], None, None]:
    for id, peptide in enumerate(peptides):
        yield peptide, False, id
        try:
            yield as_decoy(peptide, decoy_strategy), True, id
        except KeyError as e:
            print(f"No decoy for {peptide} because KeyError: {e}")


def build_massshift_dict(
    seq: str,
    decoy_strategy: DecoyStrategy,
    *,
    max_charge: int = 3,
    fragmentation_model: rustyms.FragmentationModel = rustyms.FragmentationModel.CidHcd,
) -> tuple[dict[str, float], str]:
    fw_peptide = rustyms.LinearPeptide(seq)
    fw_frags = fw_peptide.generate_theoretical_fragments(
        max_charge, fragmentation_model
    )
    fw_frags = [x for x in fw_frags if x.neutral_loss is None]
    fw_ion_dict = {f"{f.ion}^{f.charge}": f.formula.mass() / f.charge for f in fw_frags}

    rev = as_decoy(seq, decoy_strategy)
    rev_peptide = rustyms.LinearPeptide(rev)
    rev_frags = rev_peptide.generate_theoretical_fragments(
        max_charge, fragmentation_model
    )
    rev_frags = [x for x in rev_frags if x.neutral_loss is None]
    rev_ion_dict = {
        f"{f.ion}^{f.charge}": f.formula.mass() / f.charge for f in rev_frags
    }

    diff_dict = {k: rev_ion_dict[k] - fw_ion_dict[k] for k in fw_ion_dict}
    return diff_dict, rev
