import dataclasses
import json
import socket
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd
import streamlit as st
from matplotlib import pyplot as plt

from speclib_builder.base import PeptideElement
from speclib_builder.builder import DummyAnnotator, EntryBuilder
from speclib_builder.onxx_predictor import OnnxPeptideTransformerAnnotator

sample_data = "../../sample_out.json"


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

        print(chunks)
        response = b"".join(chunks).decode("utf-8")

        try:
            return json.loads(response)
        except json.JSONDecodeError:
            return response


def dataclass_from_dict(klass, d):
    # from: https://stackoverflow.com/a/54769644/4295016
    try:
        fieldtypes = {f.name: f.type for f in dataclasses.fields(klass)}
        return klass(**{f: dataclass_from_dict(fieldtypes[f], d[f]) for f in d})
    except TypeError as e:
        if "must be called with a dataclass type or instance" in str(e):
            return d
    except Exception as e:
        raise KeyError(f"Could not convert to {klass} because {e}") from e


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
    summed_intensity: list[int]
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
    ref_time_ms: list[int]

    def plot(self, min_rt_ms, max_rt_ms, vlines_ms: list[int] | None = None):
        # Make a plot grid, where each row is a different score element
        # but all share the same retention time axis
        fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(10, 8))

        rt_use = np.array(self.ref_time_ms)
        ranges = np.searchsorted(rt_use, [min_rt_ms, max_rt_ms])
        rt_plot = (rt_use[ranges[0] : ranges[1]] / 1000) / 60

        ax[0, 0].plot(rt_plot, self.ms1_coelution_score[ranges[0] : ranges[1]])
        ax[0, 0].set_title("MS1 Coelution Score")
        ax[0, 1].plot(rt_plot, self.ms1_cosine_ref_sim[ranges[0] : ranges[1]])
        ax[0, 1].set_title("MS1 Cosine Ref Sim")
        ax[0, 2].plot(rt_plot, self.ms2_lazyscore_vs_baseline[ranges[0] : ranges[1]])
        ax[0, 2].set_title("MS2 LazyScore Baseline")
        ax[1, 0].plot(rt_plot, self.ms2_coelution_score[ranges[0] : ranges[1]])
        ax[1, 0].set_title("MS2 Coelution Score")
        ax[1, 1].plot(rt_plot, self.ms2_cosine_ref_sim[ranges[0] : ranges[1]])
        ax[1, 1].set_title("MS2 Cosine Ref Sim")
        ax[1, 2].plot(rt_plot, self.ms2_lazyscore[ranges[0] : ranges[1]])
        ax[1, 2].set_title("MS2 LazyScore")

        if vlines_ms is not None:
            vlines_minutes = np.array(vlines_ms) / 1000 / 60
            for i in range(len(vlines_minutes)):
                ax[0, 0].axvline(
                    x=vlines_minutes[i], color="k", linestyle="--", alpha=0.5
                )
                ax[0, 1].axvline(
                    x=vlines_minutes[i], color="k", linestyle="--", alpha=0.5
                )
                ax[0, 2].axvline(
                    x=vlines_minutes[i], color="k", linestyle="--", alpha=0.5
                )
                ax[1, 0].axvline(
                    x=vlines_minutes[i], color="k", linestyle="--", alpha=0.5
                )
                ax[1, 1].axvline(
                    x=vlines_minutes[i], color="k", linestyle="--", alpha=0.5
                )
                ax[1, 2].axvline(
                    x=vlines_minutes[i], color="k", linestyle="--", alpha=0.5
                )

        for i in range(2):
            for j in range(3):
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
    main_score: float
    delta_next: float
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


def input_compoinent() -> PeptideElement:
    options = ["Sequence", "Examples"]
    option = st.selectbox("Input", options)
    if option == "Sequence":
        peptide = PeptideElement(
            peptide=st.text_input("Peptide", "TLSDYNIQK"),
            charge=st.slider("Charge", 2, 5, 2),
            nce=st.slider("NCE", 10, 50, 35),
            decoy=st.checkbox("Decoy"),
        )
        return peptide
    if option == "Examples":
        examples = [("TLSDYNIQK", 2), ("ESTLHLVLR", 2), ("DIKPENLLLGSAGELK", 3)]
        edict = {k: (k, v) for k, v in examples}
        picked = st.selectbox("Example", list(edict.keys()))

        picked, charge = edict[picked]
        peptide = PeptideElement(
            peptide=picked,
            charge=charge,
            nce=35,
            decoy=True,
        )
        return peptide

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
    peptide = input_compoinent()

    if ml_annotator is not None:
        st.info("Using ML annotator")
        annotator = ml_annotator
    else:
        st.warning("Using dummy annotator")

    query_data = entry_builder.build_entry(annotator.model(peptide)).as_rt_entry()
    with st.expander("Query data"):
        st.write(query_data)

    stime = time.monotonic()
    data = query_server(host, int(port), query_data)
    if isinstance(data, str):
        st.error(f"Query failed: {data}")
        st.write(data)

    etime = time.monotonic()
    st.write(f"Query took {etime - stime} seconds")

    # with open(sample_data) as f:
    #     stime = time.monotonic()
    #     data = json.load(f)
    #     etime = time.monotonic()
    #     st.write(f"Loaded data in {etime-stime} seconds")

    if data["status"] != "success":
        st.write(str(data))
        st.stop()

    res: Response = dataclass_from_dict(Response, data)
    main_score = res.data.search_results.main_score
    st.subheader("Main Score: " + str(main_score))

    st.dataframe(res.data.search_results.as_table(), use_container_width=True)

    best_rt = res.data.search_results.obs_rt_seconds / 60
    min_rt = (res.data.main_score_elements.min_rt() / 1000) / 60
    max_rt = (res.data.main_score_elements.max_rt() / 1000) / 60
    min_rt_show = st.slider("Minimum retention time (minutes)", min_rt, max_rt, min_rt)
    max_rt_show = st.slider("Maximum retention time (minutes)", min_rt, max_rt, max_rt)

    fig = res.data.plot_main_score(min_rt_show * 1000 * 60, max_rt_show * 1000 * 60)
    plt.axvline(x=best_rt, color="red", alpha=0.5)
    st.pyplot(fig, clear_figure=True, use_container_width=True)

    fig = res.data.main_score_elements.plot(
        min_rt_show * 1000 * 60,
        max_rt_show * 1000 * 60,
        vlines_ms=[best_rt * 1000 * 60],
    )
    st.pyplot(fig, clear_figure=True, use_container_width=True)

    fig = res.data.extractions.ms1_arrays.plot_transition(
        min_rt_show * 1000 * 60, max_rt_show * 1000 * 60
    )
    plt.axvline(x=best_rt, color="red", alpha=0.5)
    st.pyplot(fig, clear_figure=True, use_container_width=True)

    fig = res.data.extractions.ms2_arrays.plot_transition(
        min_rt_show * 1000 * 60, max_rt_show * 1000 * 60
    )
    plt.axvline(x=best_rt, color="red", alpha=0.5)
    st.pyplot(fig, clear_figure=True, use_container_width=True)


if __name__ == "__main__":
    main()
