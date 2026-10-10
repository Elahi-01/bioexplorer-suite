import io
import json
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy.stats import ttest_ind
from src.data_loader import read_csv, validate_expression, validate_metadata, validate_amr, DATA, load_gene, load_amr, dataframe_fingerprint
from src.differential import differential, classify, bh_adjust
from src.amr_analysis import summarize

def expression_fixture():
    return pd.DataFrame([[1.0, 1.2, 0.8, 3.0, 3.3, 2.7], [5.0, 5.1, 4.9, 2.0, 2.2, 1.8], [4, 4, 4, 4, 4, 4], [2, 2, 2, 8, 8, 8]], index=pd.Index(["up", "down", "constant", "no_variance"], name="Feature"), columns=["H1", "H2", "H3", "D1", "D2", "D3"])

def test_welch_effect_direction_and_non_testable_features():
    x = expression_fixture()
    r = differential(x, ["H1", "H2", "H3"], ["D1", "D2", "D3"]).set_index("Feature")
    assert r.loc["up", "log2FC"] == pytest.approx(2)
    assert r.loc["down", "log2FC"] == pytest.approx(-3)
    assert r.loc["up", "p_value"] == pytest.approx(ttest_ind([3, 3.3, 2.7], [1, 1.2, 0.8], equal_var=False).pvalue)
    assert r.loc["constant", "FDR"] != r.loc["constant", "FDR"]
    assert not r.loc["no_variance", "Testable"]
    classified = classify(r.reset_index(), 0.05, 1).set_index("Feature")
    assert classified.loc["up", "Direction"] == "Upregulated"
    assert classified.loc["down", "Direction"] == "Downregulated"
    assert classified.loc["no_variance", "Direction"] == "Not testable"

def test_bh_against_published_scipy_example_and_nan_exclusion():
    expected = [0.6, 0.6, 0.2, 0.004, np.nan]
    np.testing.assert_allclose(bh_adjust([0.5, 0.6, 0.1, 0.001, np.nan]), expected, equal_nan=True)
    with pytest.raises(ValueError):
        bh_adjust([1.2])

def test_reference_reversal_preserves_p_and_reverses_effect():
    x = expression_fixture()
    a = differential(x, ["H1", "H2", "H3"], ["D1", "D2", "D3"]).set_index("Feature").sort_index()
    b = differential(x, ["D1", "D2", "D3"], ["H1", "H2", "H3"]).set_index("Feature").sort_index()
    np.testing.assert_allclose(a.log2FC, -b.log2FC)
    np.testing.assert_allclose(a.p_value, b.p_value, equal_nan=True)

def test_strict_fdr_and_inclusive_effect_threshold():
    r = pd.DataFrame({"FDR": [0.05, 0.049], "Testable": [True, True], "log2FC": [2., 1.]})
    assert classify(r, 0.05, 1).Direction.tolist() == ["Not significant", "Upregulated"]

@pytest.mark.parametrize("csv_text", ["Gene,H1,H1,D1,D2\ng1,1,2,3,4", "Gene,H1,H2,D1,D2\ng1,1,,3,4", "Gene,H1,H2,D1,D2\ng1,100,200,300,400", "Gene,H1,H2,D1,D2\ng1,1,2,3,4\ng1,2,3,4,5"])
def test_bad_expression_is_rejected(csv_text):
    with pytest.raises(ValueError):
        validate_expression(read_csv(io.BytesIO(csv_text.encode())))

def test_metadata_must_match_all_sample_ids():
    with pytest.raises(ValueError):
        validate_metadata(pd.DataFrame({"Sample": ["H1", "D2"], "Group": ["H", "D"]}), ["H1", "D1"])

def test_analysis_cache_fingerprint_covers_every_input_row():
    d = pd.DataFrame({"sample": np.arange(60000, dtype=float)})
    before = dataframe_fingerprint(d)
    d.loc[59999, "sample"] += 0.001
    assert dataframe_fingerprint(d) != before

def amr_fixture():
    return pd.DataFrame({"Isolate_ID": ["a", "b", "c", "d", "e"], "Bacteria": ["Ecoli"] * 5, "Antibiotic": ["DrugA"] * 3 + ["DrugB"] * 2, "Result": ["R", "I", "S", "S", "S"], "Year": [2020, 2020, 2021, 2020, 2020]})

def test_amr_sir_denominator_and_zero_resistance_interval():
    d = validate_amr(amr_fixture())
    s = summarize(d, ["Bacteria", "Antibiotic"]).set_index("Antibiotic")
    assert s.loc["DrugA", "Tested"] == 3
    assert s.loc["DrugA", "Resistance_pct"] == pytest.approx(100 / 3)
    assert s.loc["DrugA", "Intermediate"] == 1
    assert s.loc["DrugB", "Resistance_pct"] == 0
    assert s.loc["DrugB", "CI_high"] > 0

def test_amr_counts_are_weighted_when_years_are_combined():
    d = validate_amr(amr_fixture()).query("Antibiotic == 'DrugA'")
    annual = summarize(d, ["Year"])
    combined = summarize(d, ["Bacteria", "Antibiotic"])
    assert combined.Resistance_pct.iloc[0] == pytest.approx(100 / 3)
    assert combined.Resistance_pct.iloc[0] != annual.Resistance_pct.mean()

def test_unused_categorical_combinations_are_not_observations():
    d = validate_amr(amr_fixture())
    for c in ["Bacteria", "Antibiotic", "Result"]:
        d[c] = d[c].astype("category")
    d["Bacteria"] = d.Bacteria.cat.add_categories(["Unobserved organism"])
    d["Antibiotic"] = d.Antibiotic.cat.add_categories(["Untested drug"])
    s = summarize(d, ["Bacteria", "Antibiotic", "Year"])
    assert s.Tested.gt(0).all()
    assert np.isfinite(s[["Resistance_pct", "CI_low", "CI_high"]]).all().all()
    assert s.Tested.sum() == len(d)

def test_duplicate_ast_and_unknown_classifications_are_rejected():
    d = amr_fixture()
    with pytest.raises(ValueError):
        validate_amr(pd.concat([d, d.iloc[[0]]]))
    d.loc[0, "Result"] = "Unknown"
    with pytest.raises(ValueError):
        validate_amr(d)

def test_inconsistent_isolate_metadata_is_rejected():
    d = amr_fixture()
    d.loc[3, "Isolate_ID"] = "a"
    d.loc[3, "Year"] = 2022
    with pytest.raises(ValueError):
        validate_amr(d)

def test_real_bundles_match_manifest_and_have_no_duplicate_ast():
    manifest = json.loads((DATA / "manifest.json").read_text())
    for name, info in manifest["files"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == info["sha256"]
    gene, meta, annotation = load_gene()
    assert gene.shape == (54675, 121)
    assert meta.Group.value_counts().to_dict() == {"Breast cancer": 104, "Normal": 17}
    assert annotation.Feature.is_unique
    amr = load_amr()
    assert len(amr) == manifest["amr"]["bundled_ast_rows"]
    assert amr.Isolate_ID.nunique() == manifest["amr"]["isolates"]
    assert not amr.duplicated(["Isolate_ID", "Bacteria", "Antibiotic"]).any()
