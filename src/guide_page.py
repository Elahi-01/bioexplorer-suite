import streamlit as st
from .ui import hero

def render():
    hero("A clearer path from data to discovery.", "Two focused Streamlit tools, bundled real datasets, and charts that explain what you are looking at.", ["GENEDIFF", "AMR EXPLORER", "OPEN SOURCE CODE"])
    a, b = st.columns(2)
    with a, st.container(border=True):
        st.subheader("01 · GeneDiff")
        st.write("Compare normalized log2 gene expression across two independent groups. Start with the real breast-tissue GEO study or upload your own matrix.")
        st.markdown("**Explore:** volcano plot, clustered feature heatmap, individual feature distributions, sample PCA, sortable results and CSV exports.")
    with b, st.container(border=True):
        st.subheader("02 · AMR Explorer")
        st.write("Explore source-reported susceptibility classifications. Start with the real FDA canine surveillance dataset or bring your own public study records.")
        st.markdown("**Explore:** resistance percentages with tested counts, uncertainty intervals, annual trends, classification mix, regional heatmaps and CSV exports.")
    st.subheader("Start in three steps")
    st.markdown("1. Choose a tool in the sidebar and open its bundled real dataset.\n2. Adjust the cohort and thresholds, then hover over charts to inspect the evidence.\n3. Export the results, interactive plots and settings so your comparison can be reproduced.")
    st.subheader("How to interpret the visuals")
    st.markdown("- **Coral** marks higher expression or resistant classifications. **Teal** marks lower expression or susceptible classifications.\n- **Gene heatmap colors** compare each feature with its own mean. They do not compare absolute expression between different genes.\n- **Blank AMR heatmap cells** indicate missing or suppressed data. They are not evidence of zero resistance.\n- **FDR and effect size** answer different questions. Use both, and inspect sample quality before drawing conclusions.")
    with st.expander("Run or deploy your own copy"):
        st.code("python -m venv .venv\n# Activate the environment, then:\npython -m pip install -r requirements.txt\npython -m streamlit run app.py", language="bash")
        st.write("Push the extracted project contents to a GitHub repository. In Streamlit Community Cloud, select that repository and app.py as the entrypoint. Choose Python 3.11 in Advanced settings. GitHub Pages does not execute Streamlit Python apps.")
        st.markdown("[Official Streamlit deployment guide](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app)")
