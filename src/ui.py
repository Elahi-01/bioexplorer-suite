from __future__ import annotations
import html
import streamlit as st

TEAL = "#0F8F87"
CORAL = "#E26B70"
NAVY = "#183C55"
MUTED = "#97A9BC"

def style():
    st.markdown("""<style>
    .block-container{max-width:1440px;padding-top:1.4rem;padding-bottom:3rem}
    h1,h2,h3{letter-spacing:-.035em;color:#142B40} h1{font-size:2.5rem!important}
    [data-testid="stSidebar"]{border-right:1px solid #E2E8F0;background:#FFFFFF}
    [data-testid="stSidebar"] .block-container{padding-top:1.8rem}
    [data-testid="stMetric"]{background:white;border:1px solid #E0E7EE;padding:20px 22px;border-radius:15px;box-shadow:0 3px 15px #17344804}
    [data-testid="stMetricLabel"]{color:#60758A;font-size:.8rem}
    [data-testid="stMetricValue"]{font-size:2rem;letter-spacing:-.035em;color:#183C55}
    .brand{display:flex;align-items:center;gap:10px;font-weight:800;font-size:20px;letter-spacing:-.04em;margin-bottom:6px}
    .brand-icon{display:inline-flex;align-items:center;justify-content:center;background:#0F8F87;color:white;width:35px;height:35px;border-radius:10px;font-size:23px}
    .brand-sub{font-size:11px;color:#778C9E;letter-spacing:.14em;margin-bottom:28px}
    .hero{background:linear-gradient(110deg,#12384E,#1D5560);border-radius:20px;padding:28px 32px;margin:12px 0 22px;position:relative;overflow:hidden;color:white}
    .hero h1{color:white;margin:10px 0 0;font-size:2.7rem!important;padding:0}
    .hero p{color:#C3DCE2;max-width:750px;font-size:15px;line-height:1.6;margin:10px 0 0}
    .eyebrow{font-size:10px;font-weight:700;letter-spacing:.15em;text-transform:uppercase;color:#82DDD1}
    .pill{display:inline-block;padding:5px 10px;border:1px solid #FFFFFF2B;border-radius:20px;background:#FFFFFF0C;color:#D8EEF0;font-size:10px;letter-spacing:.02em;margin:13px 7px 0 0}
    .hero-orbit{position:absolute;width:240px;height:240px;border:1px solid #A6F2DA22;border-radius:50%;right:-55px;top:-70px;pointer-events:none}
    .hero-orbit:after{content:"";position:absolute;inset:28px;border:1px solid #A6F2DA22;border-radius:50%}
    .note{background:#EAF5F2;border:1px solid #D2EAE3;border-left:4px solid #0F8F87;border-radius:11px;padding:14px 18px;color:#255555;font-size:13px;line-height:1.6;margin:14px 0}
    .section-kicker{font-size:10px;color:#8A9BAB;letter-spacing:.15em;font-weight:700;margin-top:10px}
    .footer{color:#8699A9;font-size:11px;border-top:1px solid #DDE6ED;margin-top:32px;padding-top:14px;line-height:1.7}
    [data-baseweb="tab-list"]{gap:18px;border-bottom:1px solid #DCE5ED;margin:15px 0 20px}
    [data-baseweb="tab"]{padding:8px 2px 12px;font-weight:600;font-size:13px}
    [data-testid="stDataFrame"]{border-radius:12px}
    div.stButton>button,div.stDownloadButton>button{border-radius:10px}
    @media(max-width:640px){.hero{padding:22px}.hero h1{font-size:2rem!important}.hero p{font-size:13px}.block-container{padding:1rem}.hero-orbit{display:none}}
    </style>""", unsafe_allow_html=True)

def branding():
    st.sidebar.markdown('<div class="brand"><span class="brand-icon">✳</span> BioExplorer</div><div class="brand-sub">A RESEARCH WORKSPACE</div>', unsafe_allow_html=True)

def hero(name, description, tags):
    pills = "".join(f'<span class="pill">{html.escape(tag)}</span>' for tag in tags)
    st.markdown(f'<div class="hero"><div class="hero-orbit"></div><div class="eyebrow">BIOEXPLORER / DISCOVERY WORKSPACE</div><h1>{html.escape(name)}</h1><p>{html.escape(description)}</p>{pills}</div>', unsafe_allow_html=True)

def note(text):
    st.markdown(f'<div class="note">{html.escape(text)}</div>', unsafe_allow_html=True)

def footer():
    st.markdown('<div class="footer">BioExplorer · Research and education only · Source-linked public datasets<br>GeneDiff explores expression differences. AMR Explorer describes observed susceptibility patterns.</div>', unsafe_allow_html=True)

def show_plot(fig, key):
    st.plotly_chart(fig, use_container_width=True, theme=None, key=key, config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})

def plot_export(fig, filename, label="Download interactive plot"):
    st.download_button(label, fig.to_html(include_plotlyjs=True, full_html=True), filename, "text/html", key=filename)

def safe_csv(df):
    """Neutralize formula-like text in exports opened by spreadsheet programs."""
    df = df.copy()
    for c in df.select_dtypes(include=["object", "string", "category"]).columns:
        s = df[c].astype(str)
        df[c] = s.map(lambda v: "'" + v if v.lstrip().startswith(("=", "+", "-", "@", "\t", "\r")) else v)
    return df.to_csv(index=False).encode("utf-8-sig")
