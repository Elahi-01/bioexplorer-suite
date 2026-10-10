"""Descriptive S/I/R summaries; never pools rates across antibiotic panels."""
from __future__ import annotations
import numpy as np
import pandas as pd

def summarize(df, by):
    if df.empty:
        return pd.DataFrame(columns=list(by) + ["Tested", "Resistant", "Intermediate", "Susceptible", "Resistance_pct", "CI_low", "CI_high"])
    counts = df.groupby(list(by), observed=True, dropna=False).Result.value_counts().unstack(fill_value=0)
    for c in ["Resistant", "Intermediate", "Susceptible"]:
        if c not in counts:
            counts[c] = 0
    counts["Tested"] = counts[["Resistant", "Intermediate", "Susceptible"]].sum(axis=1)
    # Categorical value_counts/unstack may materialize unused combinations.
    # An unobserved group is missing, never a tested zero-resistance group.
    counts = counts.loc[counts.Tested.gt(0)].copy()
    n = counts.Tested.to_numpy(float)
    p = counts.Resistant.to_numpy(float) / n
    z = 1.959963984540054
    center = (p + z * z / (2 * n)) / (1 + z * z / n)
    radius = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    counts["Resistance_pct"] = p * 100
    counts["CI_low"] = np.maximum(0, center - radius) * 100
    counts["CI_high"] = np.minimum(1, center + radius) * 100
    return counts.reset_index()

def filter_amr(df, bacteria, antibiotics, years, states, specimens, agencies):
    return df.loc[df.Bacteria.isin(bacteria) & df.Antibiotic.isin(antibiotics) & df.Year.isin(years) & df.State.isin(states) & df.Specimen.isin(specimens) & df.Agency.isin(agencies)].copy()
