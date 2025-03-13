from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from blueprints.prioritisation_by_nomaly_scores import get_top_variants

# Expected output values for our test case
EXPECTED_OUTPUT = {
    "metric1_top_variants": 3,  # Number of variants above the threshold
    "metric1_top_gene_set": 3,  # Number of genes for those variants
    "num_rp": 497.0,
    "num_rn": 247437.0,
    "metric1_pvalue": 0.0026,
    "metric1_tpr": 0.012,
    "metric1_fpr": 1.03,
    "metric1_lrp": 0.011,
    "metric1_tp": 6,
    "metric1_tn": 685,
    "metric1_fp": 257045,
    "metric1_fn": 532,
    "metric1_threshold": 0.022,
    # Similar expectations for other metrics
    "roc_stats_mcc_top_variants": 3,
    "roc_stats_mcc_pvalue": 0.0015,
    "roc_stats_yjs_top_variants": 15,
    "roc_stats_yjs_pvalue": 0.045,
    "roc_stats_lrp_top_variants": 6,
    "roc_stats_lrp_pvalue": 0.16,
}

# Additional expected values for the protective case
EXPECTED_PROTECTIVE_OUTPUT = {
    **EXPECTED_OUTPUT,
    "roc_stats_lrn_protective_top_variants": 0,
    "roc_stats_lrn_protective_pvalue": 1.0,
    "roc_stats_lrn_protective_tp": 0.0,
    "roc_stats_lrn_protective_tn": 247353,
    "roc_stats_lrn_protective_fp": 84,
    "roc_stats_lrn_protective_fn": 497,
    "roc_stats_lrn_protective_threshold": 0.0303,
}


@pytest.fixture
def mock_cache_path(monkeypatch):
    """Mock the cache path function to use a temp dir"""

    def mock_get_cache_path(*args):
        return Path("/tmp/nonexistent_path")  # Ensure cache misses

    monkeypatch.setattr(
        "blueprints.prioritisation_by_nomaly_scores.get_cache_path", mock_get_cache_path
    )


@pytest.fixture
def mock_logger():
    """Create a mock stream logger"""
    logger = MagicMock()
    logger.info = MagicMock()
    return logger


@pytest.fixture
def mock_term_variants():
    """Create a realistic set of term variants for testing"""
    return pd.DataFrame(
        {
            "variant_id": ["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"],
            "gene": ["GENE1", "GENE2", "GENE3"],
            "hmm_score": [0.8, 0.7, 0.6],
            "aa": ["p.Ala123Val", "p.Ser456Thr", "p.Gly789Ala"],
        }
    )


@pytest.fixture
def mock_phenotype_data():
    """Create mock phenotype data with cases and controls"""
    # Create 10 individuals: 3 cases, 5 controls, 2 excluded
    eids = np.array([1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010])
    phenotypes = np.array(
        [1, 1, 1, 0, 0, 0, 0, 0, 9, 9]
    )  # 1=case, 0=control, 9=excluded
    return eids, phenotypes


@pytest.fixture
def mock_nomaly_scores():
    """Create mock Nomaly scores for cases and controls"""
    # Create scores for 3 cases and 5 controls
    # Higher scores (>0.022) indicate higher likelihood of being a case
    case_scores = np.array([0.03, 0.025, 0.01])  # 2 out of 3 cases above threshold
    control_scores = np.array(
        [0.005, 0.015, 0.03, 0.01, 0.005]
    )  # 1 out of 5 controls above threshold
    return case_scores, control_scores


@pytest.fixture
def mock_stats_data():
    """Create mock statistics data"""
    return {
        "num_rp": 3.0,  # 3 cases
        "num_rn": 5.0,  # 5 controls
        "metric1_pvalue": 0.04,
        "metric1_tp": 2,  # 2 true positives (cases above threshold)
        "metric1_fp": 4,  # 4 false positives (controls below threshold)
        "metric1_fn": 1,  # 1 false negative (case below threshold)
        "metric1_tn": 1,  # 1 true negative (control above threshold)
        "roc_stats_mcc_pvalue": 0.03,
        "roc_stats_mcc_or": 2.0,
        "roc_stats_mcc_threshold": 0.02,
        "roc_stats_mcc_tp": 2,
        "roc_stats_mcc_fp": 1,
        "roc_stats_mcc_fn": 1,
        "roc_stats_mcc_tn": 4,
        "roc_stats_yjs_pvalue": 0.02,
        "roc_stats_yjs_or": 3.0,
        "roc_stats_yjs_threshold": 0.01,
        "roc_stats_yjs_tp": 3,
        "roc_stats_yjs_fp": 2,
        "roc_stats_yjs_fn": 0,
        "roc_stats_yjs_tn": 3,
        "roc_stats_lrp_pvalue": 0.05,
        "roc_stats_lrp_or": 4.0,
        "roc_stats_lrp_threshold": 0.03,
        "roc_stats_lrp_tp": 1,
        "roc_stats_lrp_fp": 1,
        "roc_stats_lrp_fn": 2,
        "roc_stats_lrp_tn": 4,
        "roc_stats_lrn_protective_pvalue": 0.06,
        "roc_stats_lrn_protective_or": 0.5,
        "roc_stats_lrn_protective_threshold": 0.025,
        "roc_stats_lrn_protective_tp": 1,
        "roc_stats_lrn_protective_fp": 2,
        "roc_stats_lrn_protective_fn": 2,
        "roc_stats_lrn_protective_tn": 3,
    }


def test_get_top_variants_basic_functionality(
    unit_test_app,
    mock_cache_path,
    mock_logger,
    mock_term_variants,
    phenotype_service,  # Using the fixture instead of mock_phenotype_data
    stats_service,  # Using the fixture instead of mock_stats_data
    nomaly_scores_service,
    genotype_service,  # Using the fixture instead of local mock
):
    """Test the basic functionality of get_top_variants with controlled test data"""

    # Mock get_term_variants to return our test variants
    with patch(
        "blueprints.prioritisation_by_nomaly_scores.get_term_variants"
    ) as mock_get_term_variants:
        mock_get_term_variants.return_value = mock_term_variants

        # Call the function we're testing
        result = get_top_variants(
            phecode="250.2",  # Test phecode
            term="GO:0030800",  # Test term
            phenotype_service=phenotype_service,
            genotype_service=genotype_service,
            nomaly_scores_service=nomaly_scores_service,
            stats_service=stats_service,
            stream_logger=mock_logger,
            no_cache=True,  # Skip cache lookups
        )

        # Verify the function processed the data correctly

        # 1. Check that basic statistics were preserved
        stats_data = stats_service._hdf.get_stats_by_term_phecode("GO:0030800", "250.2")
        assert result["num_rp"] == stats_data["num_rp"]
        assert result["num_rn"] == stats_data["num_rn"]

        # 2. Check that metric1 values were calculated correctly
        assert result["metric1_threshold"] == 0.022  # This is hardcoded in the function
        assert result["metric1_tp"] == 2  # 2 cases above threshold
        assert result["metric1_fp"] == 4  # 4 controls below threshold

        # 3. Check that derived metrics were calculated
        assert "metric1_tpr" in result  # True positive rate
        assert "metric1_fpr" in result  # False positive rate
        assert "metric1_lrp" in result  # Likelihood ratio positive

        # 4. Check that top variants were processed
        assert "metric1_top_variants" in result
        assert "metric1_top_gene_set" in result

        # 5. Check that all expected statistics were processed
        for stat in ["metric1", "roc_stats_mcc", "roc_stats_yjs", "roc_stats_lrp"]:
            assert f"{stat}_top_variants" in result
            assert f"{stat}_top_gene_set" in result
            assert f"{stat}_tpr" in result
            assert f"{stat}_fpr" in result
            assert f"{stat}_lrp" in result


def test_get_top_variants_protective(
    unit_test_app,
    mock_cache_path,
    mock_logger,
    mock_term_variants,
    phenotype_service,  # Using the fixture instead of mock_phenotype_data
    stats_service,  # Using the fixture instead of mock_stats_data
    nomaly_scores_service,
    genotype_service,  # Using the fixture instead of local mock
):
    """Test the protective mode of get_top_variants"""

    # Mock get_term_variants
    with patch(
        "blueprints.prioritisation_by_nomaly_scores.get_term_variants"
    ) as mock_get_term_variants:
        mock_get_term_variants.return_value = mock_term_variants

        # Call the function with protective=True
        result = get_top_variants(
            phecode="250.2",
            term="GO:0030800",
            phenotype_service=phenotype_service,
            genotype_service=genotype_service,
            nomaly_scores_service=nomaly_scores_service,
            stats_service=stats_service,
            stream_logger=mock_logger,
            protective=True,
            no_cache=True,
        )

        # Verify protective mode was processed
        assert "roc_stats_lrn_protective_top_variants" in result
        assert "roc_stats_lrn_protective_top_gene_set" in result
        assert "roc_stats_lrn_protective_tpr" in result
        assert "roc_stats_lrn_protective_fpr" in result
        assert "roc_stats_lrn_protective_lrp" in result


def test_cache_handling(unit_test_app, monkeypatch, mock_logger):
    """Test that cache is correctly used when available"""

    # Create a simple mock result
    mock_result = {
        "num_rp": 3.0,
        "num_rn": 5.0,
        "metric1_pvalue": 0.04,
        "metric1_top_variants": [{"variant_id": "test"}],
        "metric1_top_gene_set": [{"gene": "test"}],
    }

    # Mock cache functions
    mock_load = MagicMock(return_value=mock_result)
    mock_save = MagicMock()

    monkeypatch.setattr(
        "blueprints.prioritisation_by_nomaly_scores.load_cached_results", mock_load
    )
    monkeypatch.setattr(
        "blueprints.prioritisation_by_nomaly_scores.save_results_to_cache", mock_save
    )

    # Create minimal mock services
    phenotype_service = MagicMock()
    genotype_service = MagicMock()
    nomaly_scores_service = MagicMock()
    stats_service = MagicMock()

    # Call function without no_cache flag
    result = get_top_variants(
        phecode="250.2",
        term="GO:0030800",
        phenotype_service=phenotype_service,
        genotype_service=genotype_service,
        nomaly_scores_service=nomaly_scores_service,
        stats_service=stats_service,
        stream_logger=mock_logger,
    )

    # Verify cache was loaded
    mock_load.assert_called_once()
    # Save should not be called since we returned cached data
    mock_save.assert_not_called()

    # Result should match our cached data
    assert result == mock_result
