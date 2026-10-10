"""Optional standalone AMR Explorer deployment entrypoint."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import streamlit as st
from src.ui import style, branding, footer
from src.amr_page import render
st.set_page_config(page_title="AMR Explorer · Surveillance", page_icon="✳", layout="wide")
style()
branding()
render()
footer()
