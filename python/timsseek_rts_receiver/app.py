import time

import streamlit as st
from speclib_builder.builder import DummyAnnotator, EntryBuilder
from speclib_builder.onxx_predictor import OnnxPeptideTransformerAnnotator
from timsseek_rts_receiver.io import query_server
from timsseek_rts_receiver.receiver import input_component, show_results

st.set_page_config(layout="wide")


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

    query_data_target = entry_builder.build_entry(annotator.model(peptide.target))
    query_data_decoy = entry_builder.build_entry(annotator.model(peptide.decoy))
    with st.expander("Query data - target"):
        st.write(query_data_target)
    with st.expander("Query data - decoy"):
        st.write(query_data_decoy)

    stime = time.monotonic()
    data_target = query_server(host, int(port), query_data_target.model_dump())
    data_decoy = query_server(host, int(port), query_data_decoy.model_dump())

    if isinstance(data_target, str):
        st.error(f"Query failed: {data_target}")
        st.write(data_target)

    if isinstance(data_decoy, str):
        st.error(f"Query failed: {data_decoy}")
        st.write(data_decoy)

    etime = time.monotonic()
    # st.write(f"Query took {etime - stime} seconds")

    cols = st.columns(2)
    show_results(
        data_target,
        subtitle="Target results",
        key_prefix="target_",
        column=cols[0],
    )
    show_results(
        data_decoy,
        subtitle="Decoy results",
        key_prefix="decoy_",
        column=cols[1],
    )


if __name__ == "__main__":
    main()
