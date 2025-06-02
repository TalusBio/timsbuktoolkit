import json
import time
import numpy as np

import streamlit as st
from matplotlib import pyplot as plt
from pydantic import BaseModel
from pyteomics import parser as pyteomics_parser
from speclib_builder.base import PeptideElement

from .constants import BSA_FASTA_ENTRY
from .models import (
    Extractions,
    MainScoreElements,
    Response,
    ResponseData,
    SearchResults,
)


def infinite_colour_loop():
    options = ["#ff0000", "#00ff00", "#0000ff"]
    i = 0
    while True:
        if i >= len(options):
            i = 0
        yield options[i]
        i += 1


def digest(sequence: str) -> set[str]:
    unique_peptides = set()
    new_peptides = pyteomics_parser.cleave(sequence, "trypsin")
    unique_peptides.update(new_peptides)
    return unique_peptides


def digest_maybe_fasta(sequence: str) -> list[str]:
    splits = sequence.split("\n")
    splits = [x.strip() for x in splits]
    if splits[0].startswith(">"):
        splits = splits[1:]
    digests = digest("".join(splits))
    return [x for x in digests if len(x) > 5 and len(x) < 30]


class TargetDecoyPair(BaseModel):
    target: PeptideElement
    decoy: PeptideElement


def input_component() -> TargetDecoyPair:
    options = ["Sequence", "Examples", "Digest"]
    option = st.selectbox("Input", options)
    if option == "Sequence":
        st.subheader("Sequence")
        cols = st.columns(2)
        target = PeptideElement(
            peptide=cols[0].text_input("Peptide", "TLSDYNIQK"),
            charge=cols[0].slider("Charge", 2, 5, 2),
            nce=cols[0].slider("NCE", 10, 50, 35),
            decoy=cols[0].checkbox("Decoy"),
        )
        decoy = PeptideElement(
            peptide=cols[1].text_input("Peptide", "TQINYDSLK"),
            charge=cols[1].slider("Charge", 2, 5, 2, key="decoy_charge"),
            nce=cols[1].slider("NCE", 10, 50, 35, key="decoy_nce"),
            decoy=cols[1].checkbox("Decoy", key="decoy_decoy"),
        )
        return TargetDecoyPair(target=target, decoy=decoy)
    if option == "Digest":
        st.subheader("Digest")
        fasta = st.text_area("Fasta", BSA_FASTA_ENTRY)
        elems = list(digest_maybe_fasta(fasta))
        target_seq = st.selectbox("Peptide", elems, key="target_peptide_dg")
        decoy_seq = target_seq[0] + target_seq[::-1][1:-1] + target_seq[-1]
        target_charge = st.slider("Charge", 2, 5, 2, key="target_charge_dg")
        target_nce = st.slider("NCE", 10, 50, 35, key="target_nce_dg")
        target = PeptideElement(
            peptide=target_seq,
            charge=target_charge,
            nce=target_nce,
            decoy=False,
        )
        decoy = PeptideElement(
            peptide=decoy_seq,
            charge=target_charge,
            nce=target_nce,
            decoy=True,
        )
        return TargetDecoyPair(target=target, decoy=decoy)

    if option == "Examples":
        st.subheader("Examples")
        examples = [
            ("TLSDYNIQK", 2, "TLSDYNIQK"),
            ("ESTLHLVLR", 2, "ELVLHLTSR"),
            ("DIKPENLLLGSAGELK", 3, "DLEGASGLLLNEPKIK"),
            ("VTEGLTDVILYHQPDDK", 3, "KFEEFQTDMAAHEER"),
            ("IAQDLEMYGVNYFSIK", 2, "STGNFLTLTQAIDK"),
            ("TFEMSDFIVDTR", 2, "MTGLVDEAIDTK"),
            ("VIVDFSSPNIAK", 2, "ELLGQGLLLR"),
        ]
        edict = {f"{k} - {w}": (k, v, w) for k, v, w in examples}
        picked = st.selectbox("Example", list(edict.keys()))

        picked, charge, decoy = edict[picked]
        target = PeptideElement(
            peptide=picked,
            charge=charge,
            nce=35,
            decoy=False,
        )
        decoy = PeptideElement(
            peptide=decoy,
            charge=charge,
            nce=35,
            decoy=True,
        )
        return TargetDecoyPair(target=target, decoy=decoy)

    raise NotImplementedError()


def show_results(data, column, subtitle=None, key_prefix=""):
    if subtitle is not None:
        column.subheader(subtitle)

    if isinstance(data, str):
        column.error(data)
        st.stop()

    if data["status"] != "success":
        column.write(str(data))
        st.stop()

    # # breakpoint()
    # extractions = Extractions(**data["data"]["extractions"])
    # main_score_elements = MainScoreElements(**data["data"]["main_score_elements"])
    # search_results = SearchResults(**data["data"]["search_results"])
    # tmp = ResponseData(**data["data"])
    res: Response = Response(**data)
    main_score = res.data.search_results.main_score
    column.subheader("Main Score: " + str(main_score))

    column.dataframe(res.data.search_results.as_table(), use_container_width=True)

    best_rt = res.data.search_results.obs_rt_seconds / 60
    min_rt = (res.data.extractions.min_rt() / 1000) / 60
    max_rt = (res.data.extractions.max_rt() / 1000) / 60
    default_min = best_rt - 0.5
    default_max = best_rt + 0.5
    min_rt_show = column.slider(
        "Minimum retention time (minutes)",
        min_value=min_rt,
        max_value=max_rt,
        value=default_min,
        key=key_prefix + "min_rt",
    )
    max_rt_show = column.slider(
        "Maximum retention time (minutes)",
        max_value=max_rt,
        min_value=min_rt,
        value=default_max,
        key=key_prefix + "max_rt",
    )

    fig = res.data.plot_main_score(min_rt_show * 1000 * 60, max_rt_show * 1000 * 60)
    plt.axvline(x=best_rt, color="red", alpha=0.5)
    column.pyplot(fig, clear_figure=True, use_container_width=True)

    fig = res.data.main_score_elements.plot(
        min_rt_show * 1000 * 60,
        max_rt_show * 1000 * 60,
        rt_use=np.array(res.data.extractions.precursors.rts_ms),
        vlines_ms=[best_rt * 1000 * 60],
    )
    column.pyplot(fig, clear_figure=True, use_container_width=True)

    fig = res.data.extractions.precursors.plot_transition(
        min_rt_show * 1000 * 60, max_rt_show * 1000 * 60
    )
    plt.axvline(x=best_rt, color="red", alpha=0.5)
    column.pyplot(fig, clear_figure=True, use_container_width=True)

    fig = res.data.extractions.fragments.plot_transition(
        min_rt_show * 1000 * 60, max_rt_show * 1000 * 60
    )
    plt.axvline(x=best_rt, color="red", alpha=0.5)
    column.pyplot(fig, clear_figure=True, use_container_width=True)

    # Button to download data
    data_dict = res.data.model_dump()
    data_json = json.dumps(data_dict)
    column.download_button(
        label="Download data",
        data=data_json,
        file_name="data.json",
        mime="application/json",
        key=key_prefix + "download_button",
    )
