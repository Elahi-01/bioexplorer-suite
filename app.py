"""Streamlit Community Cloud entrypoint: streamlit run app.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import streamlit as st
from src.ui import style, branding, footer

st.set_page_config(page_title="BioExplorer · GeneDiff & AMR", page_icon="✳", layout="wide", initial_sidebar_state="expanded")
style()
branding()
tool = st.sidebar.radio("Workspace", ["GeneDiff", "AMR Explorer", "Quick guide"], key="workspace")
st.sidebar.divider()
if tool == "GeneDiff":
    from src.gene_page import render
elif tool == "AMR Explorer":
    from src.amr_page import render
else:
    from src.guide_page import render
try:
    render()
except (ValueError, OSError) as exc:
    st.error(f"The data could not be analyzed: {exc}")
    st.info("Check that the bundled data folder is present, or correct the uploaded CSV. See README.md for the input format.")
footer()
