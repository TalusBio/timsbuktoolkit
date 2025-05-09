import dataclasses
import json
import socket
import time
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
import streamlit as st
from matplotlib import pyplot as plt
from pyteomics import parser as pyteomics_parser
from speclib_builder.base import PeptideElement
from speclib_builder.builder import DummyAnnotator, EntryBuilder
from speclib_builder.onxx_predictor import OnnxPeptideTransformerAnnotator

sample_data = "../../sample_out.json"

st.set_page_config(layout="wide")

YOUR_PORT = 3724


def query_server(host="localhost", port: int = YOUR_PORT, query_data: dict = {}):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.connect((host, port))

        # Send query
        query = json.dumps(query_data).encode("utf-8")
        s.sendall(query)
        s.shutdown(socket.SHUT_WR)

        # Read response with timeout
        s.settimeout(2.0)  # 2 second timeout
        chunks = []

        try:
            while True:
                chunk = s.recv(2048)
                if not chunk:
                    break
                chunks.append(chunk)
        except socket.timeout:
            pass  # It's okay if we timeout after receiving data

        # print(chunks)
        response = b"".join(chunks).decode("utf-8")

        try:
            return json.loads(response)
        except json.JSONDecodeError:
            return response


def dataclass_from_dict(klass, d):
    # from: https://stackoverflow.com/a/54769644/4295016
    try:
        fieldtypes = {f.name: f.type for f in dataclasses.fields(klass)}
    except TypeError as e:
        ALLOWED_PASSTHROUGH_TYPES = (int, float, str, bool)
        if type(d) in ALLOWED_PASSTHROUGH_TYPES:
            return d

        if "must be called with a dataclass type or instance" in str(e):
            return d
    except Exception as e:
        raise KeyError(f"Could not convert to {klass} because {e}") from e

    return klass(**{f: dataclass_from_dict(fieldtypes[f], d[f]) for f in d})


def infinite_colour_loop():
    options = ["#ff0000", "#00ff00", "#0000ff"]
    i = 0
    while True:
        if i >= len(options):
            i = 0
        yield options[i]
        i += 1


@dataclass
class ArrayResponse:
    ims_means: dict[str, list]
    intensities: dict[str, list]
    mz_means: dict[str, list]
    retention_time_miliseconds: list[int]
    weighted_ims_mean: list[float]

    def plot_transition(self, min_rt, max_rt):
        # In the same figure/axes show every intensity across retention time
        rt_use = np.array(self.retention_time_miliseconds)
        ranges = np.searchsorted(rt_use, [min_rt, max_rt])
        rt_plot = (rt_use[ranges[0] : ranges[1]] / 1000) / 60

        fig, ax = plt.subplots()
        for k, v in self.intensities.items():
            ax.plot(rt_plot, v[ranges[0] : ranges[1]], label=k)

        ax.set_xlabel("Retention Time (min)")
        ax.set_ylabel("Intensity")
        ax.legend()
        fig.tight_layout()
        return fig


@dataclass
class Extractions:
    id: int
    ms1_arrays: ArrayResponse
    ms2_arrays: ArrayResponse


@dataclass
class MainScoreElements:
    ms1_coelution_score: list[float]
    ms1_cosine_ref_sim: list[float]
    ms2_coelution_score: list[float]
    ms2_cosine_ref_sim: list[float]
    ms2_lazyscore: list[float]
    ms2_lazyscore_vs_baseline: list[float]
    ms2_corr_v_gauss: list[float]
    # TODO: REMOVE
    hyperscore: list[float]
    split_lazyscore: list[float]
    # END

    ref_time_ms: list[int]
    ms2_lazyscore_vs_baseline_std: float

    def plot(self, min_rt_ms, max_rt_ms, vlines_ms: list[int] | None = None):
        # Make a plot grid, where each row is a different score element
        # but all share the same retention time axis
        ncol = 3
        nrow = 3
        fig, ax = plt.subplots(nrows=nrow, ncols=ncol, figsize=(10, 12))

        rt_use = np.array(self.ref_time_ms)
        ranges = np.searchsorted(rt_use, [min_rt_ms, max_rt_ms])
        rt_plot = (rt_use[ranges[0] : ranges[1]] / 1000) / 60

        mask_label = np.zeros(len(self.ref_time_ms), dtype=bool)
        if vlines_ms is not None:
            for vline in vlines_ms:
                local_range = np.searchsorted(rt_use, [vline - 5_000, vline + 5_000])
                mask_label[local_range[0] : local_range[1]] = True

        score_name_pairs = [
            ("MS1 Coelution Score", self.ms1_coelution_score),
            ("MS1 Cosine Ref Sim", self.ms1_cosine_ref_sim),
            ("MS2 Coelution Score", self.ms2_coelution_score),
            ("MS2 Cosine Ref Sim", self.ms2_cosine_ref_sim),
            ("MS2 LazyScore", self.ms2_lazyscore),
            ("MS2 LazyScore Baseline", self.ms2_lazyscore_vs_baseline),
            ("Hyperscore", self.hyperscore),
            ("Split LazyScore", self.split_lazyscore),
            ("MS2 Corr v Gauss", self.ms2_corr_v_gauss),
        ]

        for i, (name, score) in enumerate(score_name_pairs):
            local_score = np.array(score).astype(float)
            score_ir = local_score[mask_label]
            score_or = local_score[~mask_label]

            # Signal to noise ratio would be (max(ir) - max(or))/std(or)
            max_ir = np.nanmax(score_ir) if len(score_ir) > 0 else 0
            max_or = np.nanmax(score_or) if len(score_or) > 0 else 0
            std_or = np.nanstd(score_or) if len(score_or) > 0 else 1e-3
            snr = (max_ir - max_or) / std_or

            ax[i // ncol, i % ncol].plot(rt_plot, score[ranges[0] : ranges[1]])
            ax[i // ncol, i % ncol].set_title(name)

            ax[i // ncol, i % ncol].text(
                0.95,
                0.95,
                f"SNR: {np.round(snr, 2)}",
                transform=ax[i // ncol, i % ncol].transAxes,
                ha="right",
                va="top",
            )

        if vlines_ms is not None:
            vlines_minutes = np.array(vlines_ms) / 1000 / 60
            for li in range(len(vlines_minutes)):
                for r in range(nrow):
                    for c in range(ncol):
                        ax[r, c].axvline(
                            x=vlines_minutes[li], color="k", linestyle="--", alpha=0.5
                        )

        for i in range(nrow):
            for j in range(ncol):
                ax[i, j].set_xlabel("Retention Time (min)")

        fig.tight_layout()
        return fig

    def min_rt(self):
        return min(self.ref_time_ms)

    def max_rt(self):
        return max(self.ref_time_ms)


@dataclass
class SearchResults:
    is_target: bool
    lazyerscore: float
    norm_lazyerscore_vs_baseline: float
    lazyerscore_vs_baseline: float
    ms2_corr_v_gauss: float
    main_score: float
    delta_next: float
    delta_second_next: float
    nqueries: int
    ms1_coelution_score: float
    ms1_cosine_ref_similarity: float
    ms1_mobility_error_0: float
    ms1_mobility_error_1: float
    ms1_mobility_error_2: float
    ms1_mz_error_0: float
    ms1_mz_error_1: float
    ms1_mz_error_2: float
    ms1_summed_precursor_intensity: float
    ms2_coelution_score: float
    ms2_cosine_ref_similarity: float
    ms2_mobility_error_0: float
    ms2_mobility_error_1: float
    ms2_mobility_error_2: float
    ms2_mobility_error_3: float
    ms2_mobility_error_4: float
    ms2_mobility_error_5: float
    ms2_mobility_error_6: float
    ms2_mz_error_0: float
    ms2_mz_error_1: float
    ms2_mz_error_2: float
    ms2_mz_error_3: float
    ms2_mz_error_4: float
    ms2_mz_error_5: float
    ms2_mz_error_6: float

    ms1_inten_ratio_0: float
    ms1_inten_ratio_1: float
    ms1_inten_ratio_2: float
    ms2_inten_ratio_0: float
    ms2_inten_ratio_1: float
    ms2_inten_ratio_2: float
    ms2_inten_ratio_3: float
    ms2_inten_ratio_4: float
    ms2_inten_ratio_5: float
    ms2_inten_ratio_6: float

    delta_ms1_ms2_mobility: float
    sq_delta_ms1_ms2_mobility: float

    delta_theo_rt: float
    sq_delta_theo_rt: float

    ms2_summed_transition_intensity: float
    npeaks: int
    obs_mobility: float
    obs_rt_seconds: float
    precursor_charge: int
    precursor_mobility_query: float
    precursor_mz: float
    precursor_rt_query_seconds: float
    sequence: str

    def as_table(self):
        return pd.DataFrame({
            "key": self.__dict__.keys(),
            "value": [str(x) for x in self.__dict__.values()],
        })


@dataclass
class ResponseData:
    extractions: Extractions
    main_score_elements: MainScoreElements
    longitudinal_main_score: list[float]
    search_results: SearchResults

    def plot_main_score(self, min_rt, max_rt):
        fig, ax = plt.subplots()
        rts = np.array(self.main_score_elements.ref_time_ms)
        ranges = np.searchsorted(rts, [min_rt, max_rt])
        rt_plot = (rts[ranges[0] : ranges[1]] / 1000) / 60
        ax.plot(rt_plot, self.longitudinal_main_score[ranges[0] : ranges[1]])
        ax.set_xlabel("Retention Time (min)")
        ax.set_ylabel("Main Score")
        fig.tight_layout()
        return fig


@dataclass
class Response:
    status: str
    data: ResponseData


BSA_FASTA_ENTRY = """>sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4
MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPF
DEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEP
ERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYY
ANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVA
RLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKE
CCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRR
HPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEK
LGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLIL
NRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLP
DTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVV
STQTALA"""


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


@dataclass
class TargetDecoyPair:
    target: PeptideElement
    decoy: PeptideElement


def input_component() -> PeptideElement:
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
        return TargetDecoyPair(target, decoy)
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
        return TargetDecoyPair(target, decoy)

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
        return TargetDecoyPair(target, decoy)

    raise NotImplementedError()


def main():
    st.title("TimsSeek RTs")
    st.markdown("This is a demo of how to use the timsseek rts server.")

    host = st.text_input("Host", "localhost")
    port = st.text_input("Port", "3724")

    entry_builder = EntryBuilder(
        min_mz=100,
        max_mz=2000,
        max_ions_keep=10,
        min_ion_mz=250,
        max_ion_mz=2000,
        min_ions=3,
    )

    annotator = DummyAnnotator()
    ml_annotator = OnnxPeptideTransformerAnnotator.get_default()
    peptide = input_component()

    if ml_annotator is not None:
        st.info("Using ML annotator")
        annotator = ml_annotator
    else:
        st.warning("Using dummy annotator")

    query_data_target = entry_builder.build_entry(
        annotator.model(peptide.target)
    ).as_rt_entry()
    query_data_decoy = entry_builder.build_entry(
        annotator.model(peptide.decoy)
    ).as_rt_entry()
    with st.expander("Query data - target"):
        st.write(query_data_target)
    with st.expander("Query data - decoy"):
        st.write(query_data_decoy)

    stime = time.monotonic()
    data_target = query_server(host, int(port), query_data_target)
    data_decoy = query_server(host, int(port), query_data_decoy)

    if isinstance(data_target, str):
        st.error(f"Query failed: {data_target}")
        st.write(data_target)

    if isinstance(data_decoy, str):
        st.error(f"Query failed: {data_decoy}")
        st.write(data_decoy)

    etime = time.monotonic()
    # st.write(f"Query took {etime - stime} seconds")

    # with open(sample_data) as f:
    #     stime = time.monotonic()
    #     data = json.load(f)
    #     etime = time.monotonic()
    #     st.write(f"Loaded data in {etime-stime} seconds")

    cols = st.columns(2)
    show_results(
        data_target, subtitle="Target results", key_prefix="target_", column=cols[0]
    )
    show_results(
        data_decoy, subtitle="Decoy results", key_prefix="decoy_", column=cols[1]
    )


def show_results(data, column, subtitle=None, key_prefix=""):
    if subtitle is not None:
        column.subheader(subtitle)

    if isinstance(data, str):
        column.error(data)
        st.stop()

    if data["status"] != "success":
        column.write(str(data))
        st.stop()

    res: Response = dataclass_from_dict(Response, data)
    main_score = res.data.search_results.main_score
    column.subheader("Main Score: " + str(main_score))

    column.dataframe(res.data.search_results.as_table(), use_container_width=True)

    best_rt = res.data.search_results.obs_rt_seconds / 60
    min_rt = (res.data.main_score_elements.min_rt() / 1000) / 60
    max_rt = (res.data.main_score_elements.max_rt() / 1000) / 60
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
        vlines_ms=[best_rt * 1000 * 60],
    )
    column.pyplot(fig, clear_figure=True, use_container_width=True)

    fig = res.data.extractions.ms1_arrays.plot_transition(
        min_rt_show * 1000 * 60, max_rt_show * 1000 * 60
    )
    plt.axvline(x=best_rt, color="red", alpha=0.5)
    column.pyplot(fig, clear_figure=True, use_container_width=True)

    fig = res.data.extractions.ms2_arrays.plot_transition(
        min_rt_show * 1000 * 60, max_rt_show * 1000 * 60
    )
    plt.axvline(x=best_rt, color="red", alpha=0.5)
    column.pyplot(fig, clear_figure=True, use_container_width=True)

    # Button to download data
    data_dict = asdict(res.data)
    data_json = json.dumps(data_dict)
    column.download_button(
        label="Download data",
        data=data_json,
        file_name="data.json",
        mime="application/json",
        key=key_prefix + "download_button",
    )


if __name__ == "__main__":
    main()
