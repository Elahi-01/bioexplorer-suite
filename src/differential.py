"""Independent-group exploratory analysis of normalized log2 expression."""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from scipy import stats

def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    valid = np.isfinite(p)
    if ((p[valid] < 0) | (p[valid] > 1)).any():
        raise ValueError("P-values must be between zero and one.")
    out = np.full_like(p, np.nan)
    if valid.any():
        out[valid] = stats.false_discovery_control(p[valid], method="bh")
    return out

def differential(expression, reference, target, annotation=None):
    if len(reference) < 2 or len(target) < 2:
        raise ValueError("Select at least two independent biological samples per group.")
    if set(reference) & set(target):
        raise ValueError("A sample cannot belong to both comparison groups.")
    if len(reference) != len(set(reference)) or len(target) != len(set(target)):
        raise ValueError("Sample lists cannot contain duplicates.")
    if not set(reference + target).issubset(expression.columns):
        raise ValueError("A selected sample is missing from the expression matrix.")
    a = expression.loc[:, reference].to_numpy(float)
    b = expression.loc[:, target].to_numpy(float)
    # Both groups without measurable within-group variation have no defensible
    # Welch variance estimate; they remain visible but receive no p/FDR value.
    testable = (a.var(axis=1, ddof=1) > 1e-12) | (b.var(axis=1, ddof=1) > 1e-12)
    p = np.full(len(expression), np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        p[testable] = stats.ttest_ind(b[testable], a[testable], axis=1, equal_var=False, alternative="two-sided").pvalue
    testable &= np.isfinite(p)
    p[~testable] = np.nan
    result = pd.DataFrame({"Feature": expression.index, "Mean_reference": a.mean(axis=1), "Mean_target": b.mean(axis=1), "log2FC": b.mean(axis=1) - a.mean(axis=1), "p_value": p, "FDR": bh_adjust(p), "Testable": testable})
    if annotation is not None and not annotation.empty:
        result = result.merge(annotation[["Feature", "Gene_symbol"]], on="Feature", how="left", validate="one_to_one")
    else:
        result["Gene_symbol"] = ""
    result["Gene_symbol"] = result.Gene_symbol.fillna("").replace("---", "")
    return result.sort_values(["FDR", "p_value", "Feature"], na_position="last").reset_index(drop=True)

def classify(result, fdr=0.05, fold=1.0):
    if not 0 < fdr <= 1 or fold < 0:
        raise ValueError("Use 0 < FDR threshold <= 1 and a non-negative fold-change threshold.")
    out = result.copy()
    sig = out.Testable & out.FDR.lt(fdr) & out.log2FC.abs().ge(fold)
    out["Direction"] = np.select([~out.Testable, sig & out.log2FC.gt(0), sig & out.log2FC.lt(0)], ["Not testable", "Upregulated", "Downregulated"], default="Not significant")
    return out

def sample_pca(expression, max_features=2000):
    selected = expression.var(axis=1).nlargest(min(max_features, len(expression))).index
    x = expression.loc[selected].to_numpy(float).T
    x -= x.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :2] * s[:2]
    explained = s[:2] ** 2 / max(float(np.sum(s ** 2)), 1e-15)
    return scores, explained, len(selected)
