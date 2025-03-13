import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch
from queue import Queue

from blueprints.prioritisation_by_nomaly_scores import (
    StreamLogger,
)

from blueprints.prioritisation_by_nomaly_scores_refactored import (
    read_nomaly_filtered_genotypes_new,
    fetch_phenotype_data,
    fetch_nomaly_scores,
    fetch_stats_data,
    compute_derived_stats,
    process_statistic,
    term_variant_prioritisation,
    process_gene_level_stats,
)


@pytest.fixture
def mock_variant_scores():
    """Test variant scores matching our test HDF5 variants."""
    return pd.DataFrame(
        {"VS00": [0.1, 0.2, 0.3], "VS01": [0.4, 0.5, 0.6], "VS11": [0.7, 0.8, 0.9]},
        index=["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"],
    )


@pytest.fixture
def mock_term_variants():
    """Test term variants matching our test HDF5 variants."""
    return pd.DataFrame(
        {
            "variant_id": ["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"],
            "gene": ["GENE1", "GENE2", "GENE3"],
            "hmm_score": [0.5, 0.6, 0.7],
            "aa": ["p.Ala123Val", "p.Ser456Thr", "p.Gly789Ala"],
        }
    )


def test_read_nomaly_filtered_genotypes_new(genotype_service):
    """Test the new genotype reading function."""
    sorted_eids = np.array([1001, 1002, 1003])
    variants = ["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"]

    result = read_nomaly_filtered_genotypes_new(sorted_eids, variants, genotype_service)

    assert result is not None
    assert np.array_equal(result["row_eids"], sorted_eids)
    assert result["col_variants"] == variants
    assert len(result["error_variants"]) == 0
    # The actual data shape will depend on the genotype_service fixture


def test_fetch_phenotype_data(phenotype_service):
    """Test fetching phenotype data."""
    phecode = "571.5"

    # Mock the phenotype_service.get_cases_for_phecode method
    mock_eids = np.array([1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010])
    mock_phenotypes = np.array(
        [1, 1, 1, 0, 0, 0, 0, 0, 9, 9]
    )  # 1=case, 0=control, 9=excluded

    phenotype_service.get_cases_for_phecode.return_value = (mock_eids, mock_phenotypes)

    case_eids, control_eids, exclude_eids, all_eids, phenotypes = fetch_phenotype_data(
        phecode, phenotype_service
    )

    # Basic validation
    assert len(case_eids) == 3  # 3 cases (phenotype=1)
    assert len(control_eids) == 5  # 5 controls (phenotype=0)
    assert len(exclude_eids) == 2  # 2 excluded (phenotype=9)
    assert len(all_eids) == len(phenotypes)
    assert np.all(np.diff(case_eids) > 0), "case_eids should be sorted"


def test_fetch_nomaly_scores(nomaly_scores_service):
    """Test fetching Nomaly scores."""
    case_eids = np.array([1001, 1002, 1003])
    control_eids = np.array([2001, 2002, 2003, 2004, 2005])
    term = "TEST:001"

    # Mock the nomaly_scores_service.get_scores_by_eids_unsorted method
    case_scores = np.array([0.03, 0.025, 0.01])
    control_scores = np.array([0.005, 0.015, 0.03, 0.01, 0.005])

    nomaly_scores_service.get_scores_by_eids_unsorted.side_effect = [
        case_scores,
        control_scores,
    ]

    result_case_scores, result_control_scores = fetch_nomaly_scores(
        case_eids, control_eids, nomaly_scores_service, term
    )

    assert np.array_equal(result_case_scores, case_scores)
    assert np.array_equal(result_control_scores, control_scores)
    assert len(result_case_scores) == len(case_eids)
    assert len(result_control_scores) == len(control_eids)


def test_fetch_stats_data(stats_service):
    """Test fetching statistics data."""
    term = "TEST:001"
    phecode = "571.5"

    with patch.object(
        stats_service,
        "get_stats_by_term_phecode",
        return_value={
            "num_rp": 3.0,
            "num_rn": 5.0,
            "metric1_pvalue": 0.04,
            "roc_stats_mcc_pvalue": 0.03,
            "roc_stats_mcc_or": 2.0,
            "roc_stats_mcc_threshold": 0.02,
            "roc_stats_mcc_tp": 2,
            "roc_stats_mcc_fp": 1,
            "roc_stats_mcc_fn": 1,
            "roc_stats_mcc_tn": 4,
        },
    ):
        stats = fetch_stats_data(term, phecode, stats_service)

        assert stats["num_rp"] == 3.0
        assert stats["num_rn"] == 5.0
        assert stats["metric1_pvalue"] == 0.04
        assert stats["roc_stats_mcc_pvalue"] == 0.03


def test_compute_derived_stats():
    """Test computing derived statistics."""
    stats = {
        "num_rp": 3.0,
        "num_rn": 5.0,
        "metric1_pvalue": 0.04,
    }

    case_eids = np.array([1001, 1002, 1003])
    control_eids = np.array([2001, 2002, 2003, 2004, 2005])
    case_scores = np.array([0.03, 0.025, 0.01])  # 2 above threshold
    control_scores = np.array([0.005, 0.015, 0.03, 0.01, 0.005])  # 1 above threshold
    phecode = "571.5"

    result = compute_derived_stats(
        stats, case_eids, control_eids, case_scores, control_scores, phecode
    )

    assert result["metric1_threshold"] == 0.022
    assert result["metric1_tp"] == 2  # 2 cases above threshold
    assert result["metric1_fp"] == 1  # 1 control above threshold
    assert result["metric1_fn"] == 1  # 1 case below threshold
    assert result["metric1_tn"] == 4  # 4 controls below threshold


def test_term_variant_prioritisation(
    genotype_service, mock_variant_scores, mock_term_variants
):
    """Test the term_variant_prioritisation function."""
    # In the refactored version, we pass the term_variants DataFrame directly
    # instead of the term string
    sorted_eids = np.array([1001, 1002, 1003])

    # Mock the genotype_service.get_genotypes method
    mock_genotypes = np.array(
        [
            [2, 1, 0],  # Individual 1001
            [0, 2, 1],  # Individual 1002
            [1, 0, 2],  # Individual 1003
        ]
    )
    genotype_service.get_genotypes.return_value = (
        mock_genotypes.T
    )  # Function expects transposed matrix

    # Rename the columns to match what the function expects
    mock_variant_scores.columns = ["vs00", "vs01", "vs11"]

    # Mock the term_variant_prioritisation function to use our mock_variant_scores
    with patch(
        "blueprints.prioritisation_by_nomaly_scores_refactored.term_variant_prioritisation",
        side_effect=lambda term_variants,
        eids,
        genotype_service,
        stream_logger=None: pd.DataFrame(
            {
                "variant_id": term_variants["variant_id"],
                "gene": term_variants["gene"],
                "hmm_score": term_variants["hmm_score"],
                "vs": [0.3, 0.4, 0.5],
                "num_individuals": [2, 1, 3],
            }
        ),
    ):
        result = term_variant_prioritisation(
            mock_term_variants, sorted_eids, genotype_service
        )

        # Test the actual business logic:
        assert isinstance(result, pd.DataFrame)

        # Verify required columns are present
        assert all(col in result.columns for col in ["gene", "variant_id", "hmm_score"])

        # Check that gene is a string (after our fix)
        if len(result) > 0:
            assert isinstance(result["gene"].iloc[0], str)


def test_process_gene_level_stats():
    """Test processing gene-level statistics."""
    # Create a sample DataFrame of top variants
    top_variants = pd.DataFrame(
        {
            "variant_id": ["1_123_A/G", "1_456_G/A", "2_789_C/T"],
            "gene": ["GENE1", "GENE1, GENE2", "GENE3"],
            "hmm_score": [0.8, 0.7, 0.6],
            "vs": [0.3, 0.4, 0.5],
            "num_individuals": [10, 15, 20],
        }
    )

    result = process_gene_level_stats(top_variants)

    assert isinstance(result, list)
    assert len(result) > 0
    assert "gene" in result[0]
    assert "variant_id" in result[0]
    assert "hmm_score" in result[0]
    assert "total_vs" in result[0]
    assert "variant_num" in result[0]
    assert "num_individuals" in result[0]


def test_stream_logger():
    """Test the StreamLogger class."""
    queue = Queue()
    logger = StreamLogger(queue)

    test_message = "Test progress message"
    logger.info(test_message)

    message = queue.get()
    assert message["type"] == "progress"
    assert message["data"] == test_message
