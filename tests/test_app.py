from pathlib import Path
from streamlit.testing.v1 import AppTest

APP = Path(__file__).resolve().parents[1] / "app.py"

def test_gene_real_data_thresholds_and_no_significant_state():
    at = AppTest.from_file(str(APP), default_timeout=90).run()
    assert not at.exception
    assert not at.error
    assert at.metric[0].value == "54,675"
    initial = int(at.metric[1].value.replace(",", ""))
    assert initial > 0
    at.slider(key="gene_fold").set_value(4.0).run()
    assert not at.exception
    assert int(at.metric[1].value.replace(",", "")) <= initial
    at.radio(key="gene_cohort").set_value("Quick comparison · 10 + 10").run()
    assert not at.exception
    at.radio(key="gene_source").set_value("Upload my CSV").run()
    assert not at.exception
    at.selectbox(key="gene_scale").set_value("Raw RNA-seq counts").run()
    assert not at.exception
    assert any("count-aware" in x.value for x in at.info)

def test_amr_filters_multi_organism_regional_and_empty_state():
    at = AppTest.from_file(str(APP), default_timeout=90)
    at.session_state["workspace"] = "AMR Explorer"
    at.run()
    assert not at.exception
    assert not at.error
    assert int(at.metric[0].value.replace(",", "")) > 0
    at.multiselect(key="amr_bacteria").set_value(["Escherichia coli", "Klebsiella spp"]).run()
    assert not at.exception
    at.radio(key="amr_heat_layout").set_value("State / region × antibiotic").run()
    assert not at.exception
    at.multiselect(key="amr_drugs").set_value([]).run()
    assert not at.exception
    assert any("No records match" in x.value for x in at.info)

def test_guide_and_standalone_entrypoints():
    at = AppTest.from_file(str(APP), default_timeout=90)
    at.session_state["workspace"] = "Quick guide"
    at.run()
    assert not at.exception
    for entry in ["genediff_app.py", "amr_app.py"]:
        at = AppTest.from_file(str(APP.parent / entry), default_timeout=90).run()
        assert not at.exception
        assert not at.error
