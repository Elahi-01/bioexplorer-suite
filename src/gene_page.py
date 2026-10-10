from __future__ import annotations
import json
import streamlit as st
import pandas as pd
from .data_loader import DATA, load_gene, read_csv, validate_expression, validate_metadata, dataframe_fingerprint
from .differential import differential, classify, sample_pca
from .plots import volcano, heatmap, gene_box, pca_plot, distribution_plot
from .ui import hero, note, show_plot, safe_csv, plot_export

@st.cache_data(show_spinner="Loading the public GEO dataset...", max_entries=1)
def bundled():
    return load_gene()

@st.cache_data(show_spinner="Computing Welch tests and FDR...", max_entries=3, hash_funcs={pd.DataFrame: dataframe_fingerprint})
def analyze(expression, reference, target, annotation):
    return differential(expression, reference, target, annotation)

@st.cache_data(show_spinner=False, max_entries=3, hash_funcs={pd.DataFrame: dataframe_fingerprint})
def pca(expression):
    return sample_pca(expression)

def render():
    source = st.sidebar.radio("Expression data", ["Real GEO dataset", "Upload my CSV"], key="gene_source")
    st.sidebar.markdown("**Analysis thresholds**")
    fdr = st.sidebar.select_slider("Adjusted p-value below", options=[0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2], value=0.05, key="gene_fdr")
    fold = st.sidebar.slider("Minimum |log₂ fold change|", 0.0, 4.0, 1.0, 0.1, key="gene_fold")
    st.sidebar.caption("Both thresholds must pass. Up/down are relative to the selected reference group.")
    if source == "Real GEO dataset":
        hero("GeneDiff", "Discover expression differences between normal breast tissue and breast cancer. Explore the evidence, then inspect individual features.", ["REAL DATA · GSE42568", "NORMALIZED LOG₂ MICROARRAY", "WELCH + BH FDR"])
        expr, metadata, annotation = bundled()
        cohort = st.sidebar.radio("Cohort", ["Full study · 121 samples", "Quick comparison · 10 + 10"], key="gene_cohort")
        if cohort.startswith("Quick"):
            metadata = metadata.groupby("Group", sort=False).head(10).reset_index(drop=True)
            expr = expr.loc[:, metadata.Sample]
        note("This study uses microarrays, not RNA-seq. Each row is a probe set; several probes can map to the same gene. Values are the original submitter-provided log2 GC-RMA signals.")
    else:
        hero("GeneDiff", "Bring your normalized log2 expression matrix. Map your samples, compare independent groups, and explore every result.", ["YOUR DATA", "VALIDATED INPUT", "EXPLORATORY ANALYSIS"])
        with st.expander("CSV format and upload rules", expanded=True):
            st.code("Feature,H1,H2,H3,D1,D2,D3\nTP53,5.1,5.3,5.0,7.2,7.5,7.1", language="text")
            st.caption("Example values above are illustrative only. Metadata: Sample,Group. Use unique feature IDs and finite numbers; no missing cells.")
        c1, c2 = st.columns(2)
        matrix_file = c1.file_uploader("Expression CSV or CSV.GZ", type=["csv", "gz"], key="gene_matrix")
        metadata_file = c2.file_uploader("Sample metadata CSV (optional)", type=["csv"], key="gene_metadata")
        scale = st.selectbox("Expression scale", ["Choose the scale", "Normalized log2 expression", "Raw RNA-seq counts", "Linear / unknown scale"], key="gene_scale")
        if scale in ["Raw RNA-seq counts", "Linear / unknown scale"]:
            st.info("This analysis supports normalized log2 expression only. For RNA-seq counts, use a count-aware workflow such as DESeq2 or edgeR. Do not just log-transform raw counts and run this test.")
            return
        if matrix_file is None or scale != "Normalized log2 expression":
            st.info("Upload a matrix and identify its scale to start. The sidebar's real GEO dataset is ready to explore.")
            return
        try:
            expr = validate_expression(read_csv(matrix_file))
            annotation = None
            if metadata_file:
                metadata = validate_metadata(read_csv(metadata_file), expr.columns)
            else:
                c1, c2 = st.columns(2)
                reference = c1.multiselect("Reference samples", expr.columns.tolist(), key="manual_reference")
                target = c2.multiselect("Target samples", expr.columns.tolist(), key="manual_target")
                if set(reference) & set(target):
                    st.error("A sample is selected in both groups. Remove the overlap.")
                    return
                if len(reference) < 2 or len(target) < 2:
                    st.info("Select at least two samples per group. Three or more biological replicates per group are preferable.")
                    return
                metadata = pd.DataFrame({"Sample": reference + target, "Group": ["Reference"] * len(reference) + ["Target"] * len(target)})
                excluded = len(expr.columns) - len(metadata)
                if excluded:
                    st.caption(f"{excluded} unassigned samples will be excluded from this comparison.")
                expr = expr.loc[:, metadata.Sample]
        except ValueError as exc:
            st.error(str(exc))
            return
        if not st.checkbox("These are independent biological samples, with no paired design or covariates that need to be modeled.", key="gene_independent"):
            st.info("Confirm the study design before applying an independent-group Welch test. Paired studies and batch-adjusted comparisons need an appropriate model.")
            return
    groups = metadata.Group.drop_duplicates().tolist()
    if len(groups) < 2:
        st.error("The metadata needs at least two different groups.")
        return
    with st.container(border=True):
        a, b, c = st.columns([2, 2, 3], vertical_alignment="center")
        ref_name = a.selectbox("Reference group", groups, index=groups.index("Normal") if "Normal" in groups else 0, key="gene_ref")
        target_options = [g for g in groups if g != ref_name]
        target_name = b.selectbox("Target group", target_options, key="gene_target")
        reference = metadata.loc[metadata.Group.eq(ref_name), "Sample"].tolist()
        target = metadata.loc[metadata.Group.eq(target_name), "Sample"].tolist()
        c.markdown(f"**Comparison:** {target_name} vs {ref_name}")
        c.caption(f"{len(reference)} reference + {len(target)} target samples · two-sided independent tests")
    if min(len(reference), len(target)) < 2:
        st.error("At least two independent samples are required in each selected group.")
        return
    if min(len(reference), len(target)) < 3:
        st.warning("Only two samples in at least one group. Variance estimates and statistical power will be limited.")
    samples = reference + target
    meta = metadata.set_index("Sample").loc[samples].reset_index()
    selected = expr.loc[:, samples]
    result = classify(analyze(selected, reference, target, annotation), fdr, fold)
    sig = result[result.Direction.isin(["Upregulated", "Downregulated"])]
    up = int(result.Direction.eq("Upregulated").sum())
    down = int(result.Direction.eq("Downregulated").sum())
    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Features examined", f"{len(result):,}")
    k2.metric("Significant features", f"{len(sig):,}")
    k3.metric("Upregulated in target", f"{up:,}")
    k4.metric("Downregulated in target", f"{down:,}")
    tested = int(result.Testable.sum())
    st.caption(f"{tested:,} features received a valid test; {len(result) - tested:,} without a usable variance/test result remain untested. FDR is corrected across all {tested:,} valid tests before applying the effect-size threshold.")
    tabs = st.tabs(["Overview", "Expression heatmap", "Sample QC", "Results & export", "Understand the analysis"])
    vol = volcano(result, fdr, fold)
    with tabs[0]:
        left, right = st.columns([1.8, 1])
        with left:
            st.subheader("The differential expression landscape")
            st.caption("Right = higher expression in target. Left = lower expression. Higher = stronger statistical evidence.")
            show_plot(vol, "gene_volcano")
        with right:
            st.subheader("What the comparison shows")
            st.write(f"At **FDR < {fdr:g}** and **|log₂FC| ≥ {fold:g}**, **{len(sig):,} features** pass both criteria: **{up:,} higher** and **{down:,} lower** in {target_name}.")
            st.caption("Statistical significance does not establish a disease mechanism, a biomarker or a diagnosis.")
            st.markdown("**Highest-ranked significant features**")
            if sig.empty:
                st.info("No feature passes the current thresholds. Inspect the full table or change the exploratory thresholds.")
            else:
                st.dataframe(sig[["Feature", "Gene_symbol", "log2FC", "FDR"]].head(8), hide_index=True, width="stretch", column_config={"log2FC": st.column_config.NumberColumn(format="%.2f"), "FDR": st.column_config.NumberColumn(format="%.2e")})
    with tabs[1]:
        st.subheader("Expression patterns across samples")
        n = st.slider("Features in heatmap", 10, 70, 35, 5, key="gene_heat_n")
        heat, uses_sig = heatmap(selected, result, samples, meta.set_index("Sample").Group.to_dict(), min(n, len(result)))
        st.caption(("Top significant features by FDR. " if uses_sig else "No significant features: showing the most variable features for exploration. ") + "Rows are hierarchically ordered; columns keep reference samples before target samples. Teal = lower and coral = higher than that feature's own mean. Row z-scores are clipped to ±3 for color display.")
        show_plot(heat, "gene_heatmap")
        st.markdown("**Inspect one feature**")
        query = st.text_input("Find a gene symbol or feature ID", placeholder="e.g. TP53, ESR1, 1007_s_at", key="gene_query")
        candidates = result
        if query:
            candidates = result[result.Feature.str.contains(query, case=False, regex=False) | result.Gene_symbol.str.contains(query, case=False, regex=False)]
        candidates = candidates.head(100)
        if candidates.empty:
            st.info("No matching feature was found.")
        else:
            mapping = candidates.set_index("Feature").Gene_symbol.to_dict()
            feature = st.selectbox("Feature", candidates.Feature.tolist(), format_func=lambda x: f"{x} · {mapping[x]}" if mapping[x] else x, key="gene_feature")
            row = result.set_index("Feature").loc[feature]
            st.caption(f"{row.Direction} · log2FC {row.log2FC:.3f} · FDR {row.FDR:.3g}. Each dot is one biological sample.")
            show_plot(gene_box(selected, feature, reference, target, ref_name, target_name), "gene_box")
    with tabs[2]:
        st.subheader("Check the samples before interpreting features")
        scores, explained, used = pca(selected)
        st.caption(f"PCA uses the {used:,} most variable features, centered across samples. Group separation can also reflect batch effects or other confounders.")
        show_plot(pca_plot(scores, explained, meta), "gene_pca")
        st.markdown("**Expression distributions**")
        st.caption("Comparable distributions are expected after normalization. Box plots show quartiles and 1.5×IQR whiskers; outlier dots are omitted.")
        show_plot(distribution_plot(selected, meta), "gene_distributions")
        st.dataframe(meta, hide_index=True, width="stretch")
    with tabs[3]:
        st.subheader("Explore and export the evidence")
        a, b = st.columns([2, 1])
        search = a.text_input("Search results", placeholder="Gene symbol or feature ID", key="gene_table_search")
        directions = b.multiselect("Direction filter", ["Upregulated", "Downregulated", "Not significant", "Not testable"], default=["Upregulated", "Downregulated"], key="gene_directions")
        view = result[result.Direction.isin(directions)]
        if search:
            view = view[view.Feature.str.contains(search, case=False, regex=False) | view.Gene_symbol.str.contains(search, case=False, regex=False)]
        columns = ["Feature", "Gene_symbol", "log2FC", "p_value", "FDR", "Direction", "Mean_reference", "Mean_target", "Testable"]
        st.caption(f"{len(view):,} matching rows. The screen shows at most 5,000; downloads contain every matching row.")
        st.dataframe(view[columns].head(5000), hide_index=True, width="stretch", column_config={"FDR": st.column_config.NumberColumn(format="%.3e"), "p_value": st.column_config.NumberColumn(format="%.3e"), "log2FC": st.column_config.NumberColumn(format="%.3f")})
        a, b, c = st.columns(3)
        a.download_button("All analysis results · CSV", safe_csv(result[columns]), "genediff_all_results.csv", "text/csv", key="gene_all_csv", width="stretch")
        b.download_button("Significant features · CSV", safe_csv(sig[columns]), "genediff_significant.csv", "text/csv", key="gene_sig_csv", width="stretch")
        c.download_button("Filtered table · CSV", safe_csv(view[columns]), "genediff_filtered.csv", "text/csv", key="gene_filtered_csv", width="stretch")
        settings = {"source": "GSE42568" if source == "Real GEO dataset" else "User upload", "scale": "normalized log2", "reference_group": ref_name, "target_group": target_name, "reference_samples": reference, "target_samples": target, "method": "two-sided Welch independent t-test", "multiple_testing": "Benjamini-Hochberg across valid tests", "fdr_strictly_below": fdr, "absolute_log2fc_at_least": fold, "features_examined": len(result), "features_tested": tested, "significant": len(sig), "no_covariate_adjustment": True}
        st.download_button("Analysis settings · JSON", json.dumps(settings, indent=2), "genediff_settings.json", "application/json")
        plot_export(vol, "genediff_volcano.html")
        plot_export(heat, "genediff_heatmap.html", "Download interactive heatmap")
        if source == "Real GEO dataset":
            with st.expander("Download the real input data"):
                st.download_button("Expression matrix · CSV.GZ", (DATA / "gse42568_expression.csv.gz").read_bytes(), "gse42568_expression.csv.gz", "application/gzip")
                st.download_button("Sample metadata · CSV", (DATA / "gse42568_metadata.csv").read_bytes(), "gse42568_metadata.csv", "text/csv")
    with tabs[4]:
        st.subheader("From expression values to a comparison")
        st.markdown("""1. **Validate** finite, normalized log2 expression and unique sample/feature IDs.
2. **Compare** target and reference means. `log2FC = mean(log2 target) − mean(log2 reference)`, an estimated log2 ratio of geometric means. A +1 change corresponds to about 2×; −1 to about ½×.
3. **Test** each feature with a two-sided Welch t-test, allowing unequal variances between independent groups. Both groups without measurable within-group variance are not tested.
4. **Correct** all finite test p-values with Benjamini–Hochberg. FDR controls the expected false-discovery proportion under its assumptions; it is not the probability that one particular gene is false.
5. **Explore** features passing both FDR and absolute log2FC thresholds. Heatmaps use per-feature z-scores for display, not for testing.

This prototype does not model batch, age, tumor subtype, other covariates, paired samples or technical replicate dependence. For publishable microarray inference, consider a design-aware limma workflow. For raw RNA-seq counts use a count-aware method such as DESeq2/edgeR; this app does not run those methods.
""")
        st.markdown("**Sources:** [GEO GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568) · [Example sample and log2 scale](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1045202) · [GPL570 annotation](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) · [Welch t-test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html) · [BH FDR](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html)")
