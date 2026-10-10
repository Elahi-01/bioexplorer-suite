from __future__ import annotations
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.cluster.hierarchy import linkage, leaves_list
from .ui import TEAL, CORAL, NAVY, MUTED

COLORS = {"Upregulated": CORAL, "Downregulated": TEAL, "Not significant": "#CBD5E1", "Resistant": CORAL, "Intermediate": "#E7B85D", "Susceptible": TEAL}

def polished(fig, height=440):
    fig.update_layout(template="plotly_white", height=height, paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="white", font={"family": "Arial, sans-serif", "size": 12, "color": NAVY}, margin={"l": 40, "r": 25, "t": 30, "b": 55}, legend={"orientation": "h", "y": 1.12, "x": 0, "title": None}, hoverlabel={"bgcolor": "white", "font_size": 12})
    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(gridcolor="#EEF2F7", zeroline=False)
    return fig

def volcano(results, fdr, fold):
    fig = go.Figure()
    valid = results.loc[results.Testable].copy()
    valid["score"] = -np.log10(valid.FDR.clip(lower=1e-300))
    for direction in ["Not significant", "Downregulated", "Upregulated"]:
        d = valid[valid.Direction.eq(direction)]
        fig.add_trace(go.Scattergl(x=d.log2FC, y=d.score, mode="markers", name=direction, marker={"size": 5 if direction == "Not significant" else 6, "color": COLORS[direction], "opacity": 0.45 if direction == "Not significant" else 0.8}, customdata=d[["Feature", "Gene_symbol", "FDR", "p_value"]].to_numpy(), hovertemplate="<b>%{customdata[0]}</b> · %{customdata[1]}<br>log2FC: %{x:.3f}<br>FDR: %{customdata[2]:.3g}<br>p: %{customdata[3]:.3g}<extra></extra>"))
    fig.add_hline(y=-np.log10(fdr), line_dash="dash", line_color="#9CADBD", annotation_text=f"FDR {fdr:g}", annotation_position="top right")
    for x in [-fold, fold]:
        fig.add_vline(x=x, line_dash="dash", line_color="#9CADBD")
    for _, row in valid[valid.Direction.isin(["Upregulated", "Downregulated"])].head(4).iterrows():
        label = row.Gene_symbol.split(" /// ")[0] if row.Gene_symbol else row.Feature
        fig.add_annotation(x=row.log2FC, y=row.score, text=str(label)[:22], showarrow=True, arrowhead=0, ay=-25, font={"size": 10})
    fig.update_xaxes(title="log₂ fold change · target − reference")
    fig.update_yaxes(title="−log₁₀ adjusted p-value (FDR)")
    return polished(fig, 490)

def heatmap(expression, results, samples, sample_groups, n=35):
    sig = results[results.Direction.isin(["Upregulated", "Downregulated"])]
    significant = not sig.empty
    ids = sig.head(n).Feature.tolist() if significant else expression.loc[:, samples].var(axis=1).nlargest(n).index.tolist()
    values = expression.loc[ids, samples].to_numpy(float)
    std = values.std(axis=1, ddof=1, keepdims=True)
    z = (values - values.mean(axis=1, keepdims=True)) / np.where(std > 1e-12, std, 1)
    if len(ids) > 1:
        order = leaves_list(linkage(z, method="average", metric="euclidean"))
        z, values = z[order], values[order]
        ids = [ids[i] for i in order]
    symbols = results.set_index("Feature").Gene_symbol.to_dict()
    labels = [(symbols.get(i, "").split(" /// ")[0][:20] + " · " + i) if symbols.get(i) else i for i in ids]
    groups = [sample_groups[s] for s in samples]
    custom = np.empty((len(ids), len(samples), 3), dtype=object)
    custom[:, :, 0] = values
    custom[:, :, 1] = np.array(ids)[:, None]
    custom[:, :, 2] = np.array(groups)[None, :]
    fig = go.Figure(go.Heatmap(z=np.clip(z, -3, 3), x=samples, y=labels, customdata=custom, colorscale=[[0, TEAL], [0.5, "#F7F9FC"], [1, CORAL]], zmin=-3, zmax=3, colorbar={"title": "Row z", "thickness": 12}, hovertemplate="%{customdata[1]}<br>%{x} · %{customdata[2]}<br>log2 expression: %{customdata[0]:.3f}<br>Row z-score: %{z:.2f}<extra></extra>"))
    boundary = sum(g == groups[0] for g in groups)
    fig.add_vline(x=boundary - 0.5, line_color="white", line_width=3)
    fig.update_xaxes(tickangle=-55, title="Samples grouped as reference, then target")
    fig.update_yaxes(autorange="reversed", tickfont={"size": 9})
    return polished(fig, max(460, n * 14 + 130)), significant

def gene_box(expression, feature, reference, target, ref_name, target_name):
    samples = reference + target
    d = pd.DataFrame({"Sample": samples, "Group": [ref_name] * len(reference) + [target_name] * len(target), "Expression": expression.loc[feature, samples].values})
    fig = px.box(d, x="Group", y="Expression", color="Group", points="all", hover_data=["Sample"], color_discrete_map={ref_name: TEAL, target_name: CORAL})
    fig.update_yaxes(title="Normalized log₂ expression")
    return polished(fig, 370)

def pca_plot(scores, explained, metadata):
    d = metadata.copy()
    d["PC1"], d["PC2"] = scores[:, 0], scores[:, 1]
    fig = px.scatter(d, x="PC1", y="PC2", color="Group", hover_name="Sample", hover_data=["Title"] if "Title" in d else None, color_discrete_sequence=[TEAL, CORAL], labels={"PC1": f"PC1 · {explained[0]:.1%} variance", "PC2": f"PC2 · {explained[1]:.1%} variance"})
    fig.update_traces(marker_size=10, marker_line_width=1, marker_line_color="white")
    return polished(fig, 430)

def distribution_plot(expression, metadata):
    samples = expression.columns.tolist()
    fig = go.Figure()
    groups = metadata.set_index("Sample").Group.to_dict()
    order = metadata.Group.drop_duplicates().tolist()
    for s in samples:
        # Summary-based box plots display all expression values without sending
        # 6.6 million point coordinates to the browser.
        v = expression[s].to_numpy(float)
        q1, median, q3 = np.percentile(v, [25, 50, 75])
        spread = q3 - q1
        lower = v[v >= q1 - 1.5 * spread].min()
        upper = v[v <= q3 + 1.5 * spread].max()
        color = [TEAL, CORAL][order.index(groups[s]) % 2]
        fig.add_trace(go.Box(name=s, q1=[q1], median=[median], q3=[q3], lowerfence=[lower], upperfence=[upper], marker_color=color, showlegend=False))
    fig.update_xaxes(tickangle=-60, title="Sample")
    fig.update_yaxes(title="Normalized log₂ expression")
    return polished(fig, 380)

def resistance_bars(summary, category="Antibiotic"):
    d = summary.copy()
    d["label"] = d[category].astype(str)
    if "Bacteria" in d and category != "Bacteria" and d.Bacteria.nunique() > 1:
        d["label"] = d.Bacteria.astype(str) + " · " + d[category].astype(str)
    d = d.sort_values("Resistance_pct")
    fig = go.Figure(go.Bar(x=d.Resistance_pct, y=d.label, orientation="h", marker_color=CORAL, error_x={"type": "data", "array": d.CI_high - d.Resistance_pct, "arrayminus": d.Resistance_pct - d.CI_low, "color": NAVY, "thickness": 1}, customdata=d[["Resistant", "Tested", "CI_low", "CI_high"]].to_numpy(), hovertemplate="%{y}<br>Resistance: %{x:.1f}%<br>R / tested: %{customdata[0]} / %{customdata[1]}<br>Wilson 95% CI: %{customdata[2]:.1f}–%{customdata[3]:.1f}%<extra></extra>"))
    fig.update_xaxes(title="Resistant / tested (%)", range=[0, 105])
    return polished(fig, max(400, min(1000, len(d) * 30 + 100)))

def sir_plot(summary):
    d = summary.copy()
    d["label"] = d.Antibiotic.astype(str) if d.Bacteria.nunique() == 1 else d.Bacteria.astype(str) + " · " + d.Antibiotic.astype(str)
    fig = go.Figure()
    for c in ["Susceptible", "Intermediate", "Resistant"]:
        fig.add_trace(go.Bar(x=d.label, y=d[c] / d.Tested * 100, name=c, marker_color=COLORS[c], customdata=d[[c, "Tested"]].to_numpy(), hovertemplate=c + ": %{y:.1f}%<br>Count: %{customdata[0]} / %{customdata[1]}<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_yaxes(title="Observed AST results (%)", range=[0, 100])
    fig.update_xaxes(tickangle=-40, title="Antibiotic")
    return polished(fig, 460)

def trends(summary):
    d = summary.sort_values("Year").copy()
    d["Series"] = d.Bacteria.astype(str) + " · " + d.Antibiotic.astype(str)
    fig = px.line(d, x="Year", y="Resistance_pct", color="Series", markers=True, hover_data={"Resistant": True, "Tested": True, "CI_low": ":.1f", "CI_high": ":.1f"}, color_discrete_sequence=[CORAL, TEAL, NAVY, "#BDA2CB", "#E4B461"])
    fig.update_yaxes(title="Resistant / tested (%)", range=[0, 105])
    fig.update_xaxes(dtick=1, title="Year")
    fig.update_layout(legend={"orientation": "h", "y": -0.25})
    return polished(fig, 500)

def amr_heatmap(summary, index, min_n):
    d = summary.copy()
    d.loc[d.Tested < min_n, "Resistance_pct"] = np.nan
    p = d.pivot(index=index, columns="Antibiotic", values="Resistance_pct")
    counts = d.pivot(index=index, columns="Antibiotic", values="Tested").reindex(index=p.index, columns=p.columns)
    resistant = d.pivot(index=index, columns="Antibiotic", values="Resistant").reindex(index=p.index, columns=p.columns)
    custom = np.stack([counts.to_numpy(), resistant.to_numpy()], axis=-1)
    fig = go.Figure(go.Heatmap(z=p.values, x=p.columns, y=p.index, customdata=custom, zmin=0, zmax=100, colorscale=[[0, "#E7F4EF"], [0.5, "#F8D49B"], [1, CORAL]], colorbar={"title": "R %"}, hoverongaps=False, hovertemplate="%{y} · %{x}<br>Resistance: %{z:.1f}%<br>R / tested: %{customdata[1]} / %{customdata[0]}<extra></extra>"))
    fig.update_xaxes(tickangle=-40)
    fig.update_yaxes(autorange="reversed")
    return polished(fig, max(380, min(900, len(p) * 25 + 140)))
