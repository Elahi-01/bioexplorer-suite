"""Optional standalone GeneDiff deployment entrypoint."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import streamlit as st
from src.ui import style, branding, footer
from src.gene_page import render
st.set_page_config(page_title="GeneDiff · Expression Explorer", page_icon="✳", layout="wide")
style()
branding()
render()
footer()
