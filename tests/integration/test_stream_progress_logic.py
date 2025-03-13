import json
import pytest
from typing import Optional
from flask import Flask, current_app
from blueprints.prioritisation_by_nomaly_scores import stream_progress, prioritisation_bp

# Dummy service classes to simulate app dependencies.
class DummyService:
    def __init__(self, hdf):
        self._hdf = hdf

class DummyServices:
    def __init__(self):
        self.phenotype: Optional[DummyService] = None
        self.genotype: Optional[DummyService] = None
        self.nomaly_score: Optional[DummyService] = None
        self.stats: Optional[DummyService] = None

# Dummy variant processor that immediately executes the submitted function.
class DummyVariantProcessor:
    def submit(self, func, *args, **kwargs):
        func(*args, **kwargs)
        class DummyFuture:
            def result(self):
                return None
        return DummyFuture()

@pytest.fixture
def app():
    app = Flask(__name__)
    app.config["TESTING"] = True
    # Register the blueprint that defines the /stream_progress route.
    app.register_blueprint(prioritisation_bp)
    # Set up dummy nomaly services.
    services = DummyServices()
    services.phenotype = DummyService("dummy_pheno")
    services.genotype = DummyService("dummy_geno")
    services.nomaly_score = DummyService("dummy_nomaly")
    services.stats = DummyService("dummy_stats")
    app.extensions["nomaly_services"] = services
    return app

@pytest.fixture
def client(app):
    return app.test_client()

# Dummy data that mimics what get_top_variants would normally return.
dummy_data = {
    "metric1_top_variants": ["var1", "var2"],
    "metric1_top_gene_set": ["gene1"],
    "num_rp": 10,
    "num_rn": 5,
    "metric1_pvalue": 0.05,
    "metric1_tpr": 0.8,
    "metric1_fpr": 0.2,
    "metric1_lrp": 1.2,
    "metric1_tp": 2,
    "metric1_tn": 3,
    "metric1_fp": 1,
    "metric1_fn": 0,
    "metric1_threshold": 0.5,

    "roc_stats_mcc_top_variants": ["var3"],
    "roc_stats_mcc_top_gene_set": ["gene2"],
    "roc_stats_mcc_pvalue": 0.01,
    "roc_stats_mcc_tpr": 0.7,
    "roc_stats_mcc_fpr": 0.3,
    "roc_stats_mcc_lrp": 1.1,
    "roc_stats_mcc_tp": 1,
    "roc_stats_mcc_tn": 4,
    "roc_stats_mcc_fp": 0,
    "roc_stats_mcc_fn": 1,
    "roc_stats_mcc_threshold": 0.6,

    "roc_stats_yjs_top_variants": ["var4"],
    "roc_stats_yjs_top_gene_set": ["gene3"],
    "roc_stats_yjs_pvalue": 0.02,
    "roc_stats_yjs_tpr": 0.9,
    "roc_stats_yjs_fpr": 0.1,
    "roc_stats_yjs_lrp": 1.3,
    "roc_stats_yjs_tp": 3,
    "roc_stats_yjs_tn": 2,
    "roc_stats_yjs_fp": 1,
    "roc_stats_yjs_fn": 1,
    "roc_stats_yjs_threshold": 0.7,

    "roc_stats_lrp_top_variants": ["var5"],
    "roc_stats_lrp_top_gene_set": ["gene4"],
    "roc_stats_lrp_pvalue": 0.03,
    "roc_stats_lrp_tpr": 0.6,
    "roc_stats_lrp_fpr": 0.4,
    "roc_stats_lrp_lrp": 1.4,
    "roc_stats_lrp_tp": 1,
    "roc_stats_lrp_tn": 5,
    "roc_stats_lrp_fp": 1,
    "roc_stats_lrp_fn": 0,
    "roc_stats_lrp_threshold": 0.8,
}

# For protective requests we add extra expected fields.
dummy_data_protective = dummy_data.copy()
dummy_data_protective.update({
    "roc_stats_lrn_protective_top_variants": ["var6"],
    "roc_stats_lrn_protective_top_gene_set": ["gene5"],
    "roc_stats_lrn_protective_pvalue": 0.04,
    "roc_stats_lrn_protective_tpr": 0.85,
    "roc_stats_lrn_protective_fpr": 0.15,
    "roc_stats_lrn_protective_lrp": 1.5,
    "roc_stats_lrn_protective_tp": 2,
    "roc_stats_lrn_protective_tn": 3,
    "roc_stats_lrn_protective_fp": 0,
    "roc_stats_lrn_protective_fn": 1,
    "roc_stats_lrn_protective_threshold": 0.9,
})

# Dummy get_top_variants functions for testing.
def dummy_get_top_variants(*args, **kwargs):
    return dummy_data

def dummy_get_top_variants_protective(*args, **kwargs):
    return dummy_data_protective

def dummy_get_top_variants_error(*args, **kwargs):
    raise Exception("Test error")

# Helper function to extract JSON messages from the streaming response.
def extract_messages(response):
    msgs = []
    for chunk in response.response:
        if chunk.startswith(b"data: "):
            payload = chunk.decode("utf-8").strip()
            # Remove the "data: " prefix.
            if payload.startswith("data: "):
                payload = payload[6:].strip()
            try:
                msg = json.loads(payload)
                msgs.append(msg)
            except json.JSONDecodeError:
                continue
    return msgs

def test_stream_progress_success(app, client, monkeypatch):
    # Test the normal /stream_progress execution (non-protective path).
    monkeypatch.setattr(
        "blueprints.prioritisation_by_nomaly_scores_refactored.get_top_variants",
        dummy_get_top_variants,
    )
    monkeypatch.setattr("blueprints.prioritisation_by_nomaly_scores.get_top_variants", dummy_get_top_variants)
    monkeypatch.setattr("blueprints.prioritisation_by_nomaly_scores.variant_processor", DummyVariantProcessor())
    # Call the endpoint without setting the protective flag.
    response = client.get("/stream_progress/test/dummy")
    msgs = extract_messages(response)
    types = [m.get("type") for m in msgs]
    # Expect a results message and a done message.
    assert "results" in types
    assert types[-1] == "done"
    results_msg = next(m for m in msgs if m.get("type") == "results")
    # Check that core metric keys are present in the results message.
    assert "metric1" in results_msg["data"]
    assert "roc_stats_mcc" in results_msg["data"]
    assert "roc_stats_yjs" in results_msg["data"]
    assert "roc_stats_lrp" in results_msg["data"]

def test_stream_progress_error(app, client, monkeypatch):
    # Test the error handling path where get_top_variants raises an exception.
    monkeypatch.setattr(
        "blueprints.prioritisation_by_nomaly_scores_refactored.get_top_variants",
        dummy_get_top_variants_error,
    )
    monkeypatch.setattr("blueprints.prioritisation_by_nomaly_scores.get_top_variants", dummy_get_top_variants_error)
    monkeypatch.setattr("blueprints.prioritisation_by_nomaly_scores.variant_processor", DummyVariantProcessor())
    response = client.get("/stream_progress/test/dummy")
    msgs = extract_messages(response)
    types = [m.get("type") for m in msgs]
    # Expect an error message and then a done message.
    assert "error" in types
    assert types[-1] == "done"

def test_stream_progress_protective(app, client, monkeypatch):
    # Test the protective flag behavior.
    monkeypatch.setattr(
        "blueprints.prioritisation_by_nomaly_scores_refactored.get_top_variants",
        dummy_get_top_variants_protective,
    )
    monkeypatch.setattr("blueprints.prioritisation_by_nomaly_scores.get_top_variants", dummy_get_top_variants_protective)
    monkeypatch.setattr("blueprints.prioritisation_by_nomaly_scores.variant_processor", DummyVariantProcessor())
    # Call the endpoint with the protective parameter set to true.
    response = client.get("/stream_progress/test/dummy?protective=true")
    msgs = extract_messages(response)
    types = [m.get("type") for m in msgs]
    # There should be a results message.
    assert "results" in types
    results_msg = next(m for m in msgs if m.get("type") == "results")
    # Verify that the protective data is present.
    assert "roc_stats_lrn_protective" in results_msg["data"]
    assert results_msg["data"]["roc_stats_lrn_protective"]["top_variants"] == dummy_data_protective["roc_stats_lrn_protective_top_variants"]
    # Final message should be done.
    assert types[-1] == "done"
