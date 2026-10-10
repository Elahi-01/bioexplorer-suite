from __future__ import annotations
import json
import streamlit as st
from .data_loader import DATA, load_amr, read_csv, validate_amr
from .amr_analysis import summarize, filter_amr
from .plots import resistance_bars, sir_plot, trends, amr_heatmap
from .ui import hero, note, show_plot, safe_csv, plot_export

@st.cache_data(show_spinner="Loading FDA susceptibility records...", max_entries=1)
def bundled():
    df = load_amr()
    # Repeated labels are categorical to reduce memory on hosted instances.
    for c in ["Bacteria", "Antibiotic", "Result", "State", "Host", "Specimen", "Agency"]:
        df[c] = df[c].astype("category")
    return df

def choices(values, preferred=()):
    values = sorted(values)
    return [v for v in preferred if v in values] or values[:min(4, len(values))]

def render():
    source = st.sidebar.radio("Susceptibility data", ["Real FDA dataset", "Upload my CSV"], key="amr_source")
    if source == "Real FDA dataset":
        hero("AMR Explorer", "Explore antibiotic susceptibility patterns in real veterinary isolates. Compare tested panels, follow trends, and keep the denominator in view.", ["REAL FDA DATA", "VETERINARY · DOGS", "2017–2024 SNAPSHOT"])
        df = bundled()
        manifest = json.loads((DATA / "manifest.json").read_text())["amr"]
        note("These are clinical bacterial isolates from dogs, reported by FDA Vet-LIRN / USDA NAHLN. Region labels include U.S. states and Canada. The source reports susceptible and resistant results only; no intermediate records are present.")
    else:
        hero("AMR Explorer", "Upload source-reported susceptibility classifications and explore the tested observations. Every percentage stays linked to its tested count.", ["YOUR DATA", "S / I / R", "DESCRIPTIVE ANALYSIS"])
        with st.expander("AMR CSV format", expanded=True):
            st.code("Isolate_ID,Bacteria,Antibiotic,Result,Year,State,Host,Specimen,Agency\nISO-001,Escherichia coli,Ciprofloxacin,R,2024,Region A,Human,Blood,Study A", language="text")
            st.caption("Illustrative row only. Required: Isolate_ID, Bacteria, Antibiotic, Result. Optional: Year, State, Host, Specimen, Agency, MIC. Use one row per isolate-organism-antibiotic result and globally unique isolate IDs. Do not upload identifiable patient data.")
        upload = st.file_uploader("Susceptibility records · CSV or CSV.GZ", type=["csv", "gz"], key="amr_upload")
        if not upload:
            st.info("Upload a classified susceptibility dataset or choose the real FDA data in the sidebar.")
            return
        try:
            df = validate_amr(read_csv(upload))
        except ValueError as exc:
            st.error(str(exc))
            return
        note("The app preserves the uploaded S/I/R classifications. It does not interpret MICs, apply breakpoints or choose antibiotics. The meaning of I depends on the original standard and version.")
    st.sidebar.markdown("**Explore a cohort**")
    organisms = sorted(df.Bacteria.unique())
    selected_bacteria = st.sidebar.multiselect("Bacteria", organisms, default=choices(organisms, ["Escherichia coli"]), key="amr_bacteria")
    available = sorted(df.loc[df.Bacteria.isin(selected_bacteria), "Antibiotic"].unique())
    selected_drugs = st.sidebar.multiselect("Antibiotics", available, default=choices(available, ["Ampicillin", "Enrofloxacin", "Cephalexin", "Tetracycline", "Gentamicin", "Amoxicillin/Clavulanic Acid"]), key="amr_drugs")
    years = sorted(df.Year.unique().tolist())
    if years == [0]:
        selected_years = years
        st.sidebar.caption("Year was not supplied. Time trends are unavailable.")
    elif len(years) == 1:
        selected_years = years
        st.sidebar.caption(f"Available year: {years[0]}")
    else:
        low, high = st.sidebar.select_slider("Year range", options=years, value=(years[0], years[-1]), key="amr_years")
        selected_years = [y for y in years if low <= y <= high]
    min_n = st.sidebar.number_input("Minimum tested count for charts", min_value=1, max_value=10000, value=30, step=5, key="amr_min_n", help="Suppress small groups in summary, heatmap and trend charts. The exported table retains all groups.")
    with st.sidebar.expander("Source and location filters"):
        state_values = sorted(df.State.unique())
        selected_states = st.multiselect("State / region", state_values, default=state_values, key="amr_states")
        specimen_values = sorted(df.Specimen.unique())
        selected_specimens = st.multiselect("Specimen / infection source", specimen_values, default=specimen_values, key="amr_specimens")
        agency_values = sorted(df.Agency.unique())
        selected_agencies = st.multiselect("Reporting agency", agency_values, default=agency_values, key="amr_agencies")
    filtered = filter_amr(df, selected_bacteria, selected_drugs, selected_years, selected_states, selected_specimens, selected_agencies)
    if filtered.empty:
        st.info("No records match these filters. Select a bacterium and antibiotic, or widen the source/location filters.")
        return
    all_summary = summarize(filtered, ["Bacteria", "Antibiotic"])
    summary = all_summary[all_summary.Tested.ge(min_n)].copy()
    k1, k2, k3, k4 = st.columns(4)
    k1.metric("AST results in view", f"{len(filtered):,}")
    k2.metric("Unique isolates in view", f"{filtered.Isolate_ID.nunique():,}")
    k3.metric("Pathogens in view", f"{filtered.Bacteria.nunique()}")
    k4.metric("Antibiotics in view", f"{filtered.Antibiotic.nunique()}")
    st.caption("One isolate can contribute results for several antibiotics. AST results are not a count of independent patients. Resistance is calculated separately for each organism–antibiotic pair, never averaged across the panel.")
    tabs = st.tabs(["Resistance overview", "Trends", "Pattern heatmap", "Records & export", "Data & methods"])
    bar = None
    with tabs[0]:
        left, right = st.columns([1.8, 1])
        with left:
            st.subheader("Resistance by tested antibiotic")
            st.caption(f"Each bar is one organism–antibiotic pair with at least {min_n} tests. Whiskers show Wilson 95% intervals for the observed tested sample.")
            if summary.empty:
                st.info("No pair meets the minimum tested count. Lower the display threshold or widen the filters.")
            else:
                # Long panels are easier to read in an exported table.
                displayed = summary.sort_values("Tested", ascending=False).head(24)
                if len(summary) > 24:
                    st.caption("Showing the 24 pairs with the most tests; all pairs remain in the summary export.")
                bar = resistance_bars(displayed)
                show_plot(bar, "amr_bars")
        with right:
            st.subheader("Read the percentage correctly")
            st.markdown("**Resistance % = resistant results ÷ all classified tested results × 100.**")
            st.write("Untested antibiotics and excluded non-interpretable results do not enter the denominator. Susceptible, intermediate and resistant observations stay separate.")
            if not summary.empty:
                top = summary.sort_values(["Resistance_pct", "Tested"], ascending=False).iloc[0]
                st.markdown("**Observed pattern spotlight**")
                st.write(f"**{top.Bacteria} · {top.Antibiotic}** has {int(top.Resistant):,} resistant results out of {int(top.Tested):,} tested in this filtered dataset (**{top.Resistance_pct:.1f}%**).")
                st.caption("This describes the selected surveillance sample. It does not establish a population prevalence or a treatment recommendation.")
        st.markdown("**Classification mix**")
        if not summary.empty:
            show_plot(sir_plot(summary.sort_values("Tested", ascending=False).head(20)), "amr_sir")
        if source == "Real FDA dataset":
            st.caption("Intermediate is 0 because the FDA snapshot supplies no Intermediate-classified rows. This does not prove intermediate phenotypes are absent.")
    with tabs[1]:
        st.subheader("Track the observed pattern over time")
        if filtered.Year.nunique() < 2 or selected_years == [0]:
            st.info("At least two reported years are needed for a time trend.")
        else:
            default = choices(filtered.Antibiotic.unique(), ["Enrofloxacin", "Ampicillin", "Tetracycline"])[:3]
            trend_drugs = st.multiselect("Antibiotics in the trend chart", sorted(filtered.Antibiotic.unique()), default=default, key="amr_trend_drugs")
            annual = summarize(filtered[filtered.Antibiotic.isin(trend_drugs)], ["Bacteria", "Antibiotic", "Year"])
            annual = annual[annual.Tested.ge(min_n)]
            if annual.empty:
                st.info("No annual group meets the current minimum count.")
            else:
                # An absent/suppressed year is a gap, not a line interpolated
                # across unobserved data.
                observed_years = list(range(min(selected_years), max(selected_years) + 1))
                chunks = []
                for (bact, drug), rows in annual.groupby(["Bacteria", "Antibiotic"], observed=True):
                    rows = rows.set_index("Year").reindex(observed_years).reset_index()
                    rows["Bacteria"], rows["Antibiotic"] = str(bact), str(drug)
                    chunks.append(rows)
                import pandas as pd
                trend_figure = trends(pd.concat(chunks, ignore_index=True))
                show_plot(trend_figure, "amr_trends")
                st.caption("Each point has its own tested denominator. Gaps indicate missing or suppressed years. Changes can reflect specimen mix, laboratory participation, sampling and breakpoint changes; this app does not test a causal time trend.")
                st.dataframe(annual, hide_index=True, width="stretch", column_config={"Resistance_pct": st.column_config.NumberColumn(format="%.1f%%")})
    with tabs[2]:
        st.subheader("Where do patterns differ?")
        layout = st.radio("Compare by", ["Pathogen × antibiotic", "State / region × antibiotic"], horizontal=True, key="amr_heat_layout")
        if layout.startswith("Pathogen"):
            heat_data = all_summary
            index = "Bacteria"
        else:
            organism = st.selectbox("Pathogen for the regional comparison", sorted(filtered.Bacteria.unique()), key="amr_heat_organism")
            heat_data = summarize(filtered[filtered.Bacteria.eq(organism)], ["State", "Antibiotic"])
            index = "State"
        if heat_data.empty:
            st.info("No groups are available.")
        else:
            show_plot(amr_heatmap(heat_data, index, min_n), "amr_heat")
            st.caption(f"Color = resistant percentage. Blank cells = no observations or fewer than {min_n} tests, never 0% by assumption. Hover to inspect counts. Regions are compared within a single organism; sampling is not nationally representative.")
    with tabs[3]:
        st.subheader("Keep the records and denominators together")
        st.dataframe(all_summary, hide_index=True, width="stretch", column_config={"Resistance_pct": st.column_config.NumberColumn(format="%.2f%%"), "CI_low": st.column_config.NumberColumn(format="%.2f"), "CI_high": st.column_config.NumberColumn(format="%.2f")})
        a, b = st.columns(2)
        a.download_button("Organism–antibiotic summary · CSV", safe_csv(all_summary), "amr_summary.csv", "text/csv", key="amr_summary_csv", width="stretch")
        b.download_button("Filtered AST records · CSV", safe_csv(filtered), "amr_filtered_records.csv", "text/csv", key="amr_records_csv", width="stretch")
        settings = {"source": "FDA Vet-LIRN/NAHLN canine surveillance" if source == "Real FDA dataset" else "User upload", "bacteria": selected_bacteria, "antibiotics": selected_drugs, "years": selected_years, "states": selected_states, "specimens": selected_specimens, "agencies": selected_agencies, "minimum_chart_tested": int(min_n), "denominator": "S + I + R per organism-antibiotic group", "confidence_interval": "Wilson 95%, descriptive binomial interval", "ast_records": len(filtered), "unique_isolates": int(filtered.Isolate_ID.nunique())}
        st.download_button("Filter settings · JSON", json.dumps(settings, indent=2), "amr_settings.json", "application/json")
        if bar is not None:
            plot_export(bar, "amr_resistance.html")
        search = st.text_input("Search an isolate ID", placeholder="Search public study IDs", key="amr_record_search")
        records = filtered[filtered.Isolate_ID.str.contains(search, case=False, regex=False)] if search else filtered
        st.caption(f"{len(records):,} matching records. Preview limited to 500; CSV above exports the full filtered cohort.")
        st.dataframe(records.head(500), hide_index=True, width="stretch")
        if source == "Real FDA dataset":
            st.download_button("Complete bundled FDA dataset · CSV.GZ", (DATA / "fda_amr.csv.gz").read_bytes(), "fda_amr.csv.gz", "application/gzip")
    with tabs[4]:
        st.subheader("Public source, transparent methods")
        st.markdown("""The denominator for each group is the number of **observed classified AST results**: S + I + R. An isolate can contribute to multiple drug groups. Each isolate–organism–antibiotic combination is counted once.

Wilson 95% intervals describe binomial uncertainty in the tested sample. They do not correct for non-random sampling, clustering, repeat isolates from the same animal, laboratory effects or changes in breakpoints. Globally unique isolate IDs do not establish unique animals.

The app preserves **source-reported classifications**. MIC columns, if present, are informational. No CLSI/EUCAST breakpoint interpretation is performed. In uploaded data, I retains the source's meaning; EUCAST's I category means susceptible with increased exposure, while other standards use different terminology.
""")
        if source == "Real FDA dataset":
            st.markdown("**FDA snapshot processing**")
            st.write(f"The source workbook has {manifest['source_rows']:,} rows. The bundled dataset has {manifest['bundled_ast_rows']:,} classified AST results from {manifest['isolates']:,} isolate IDs, covering all source years 2017–2024 and all four reported organisms.")
            st.write(f"Excluded: {sum(manifest['interpretations_excluded'].values()):,} non-interpretable/missing classification rows; {manifest['ambiguous_metadata_rows_excluded']:,} rows from {manifest['ambiguous_metadata_isolates_excluded']:,} isolate IDs with conflicting source metadata; {manifest['conflicting_result_rows_excluded']:,} conflicting-result rows; {manifest['missing_key_rows_excluded']:,} rows missing key fields. Removed {manifest['same_key_same_result_duplicates_removed']:,} repeated keys with the same classification. Source MICs were not reinterpreted.")
            st.json(manifest)
            st.markdown("**Original data:** [FDA Animal Pathogen AMR data](https://www.fda.gov/media/132928/download?attachment) · [FDA NARMS source and citation guidance](https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/integrated-reportssummaries) · [Vet-LIRN](https://www.fda.gov/animal-veterinary/science-research/veterinary-laboratory-investigation-and-response-network)")
        st.markdown("[EUCAST category definitions](https://www.eucast.org/newsiandr) · Research and education only. This app does not prescribe or rank antibiotics for an individual patient.")
