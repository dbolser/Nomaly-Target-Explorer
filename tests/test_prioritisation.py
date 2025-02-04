import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch
from queue import Queue

from blueprints.prioritisation_by_nomaly_scores import (
    read_cases_for_disease_code,
    read_nomaly_filtered_genotypes,
    individual_variant_prioritisation,
    term_variant_prioritisation,
    get_top_variants,
    StreamLogger,
)


@pytest.fixture
def mock_variant_scores():
    """Test variant scores matching our test HDF5 variants."""
    return pd.DataFrame(
        {"VS00": [0.1, 0.2, 0.3], "VS01": [0.4, 0.5, 0.6], "VS11": [0.7, 0.8, 0.9]},
        index=["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"],
    )


@pytest.fixture
def mock_cases_info():
    """Test cases matching our test HDF5 individuals."""
    return {
        "phecode": "571.5",
        "Sex": "both",
        "phecode_exclude_range": [],
        "cases": [1001, 1002, 1003],
        "exclude": [],
        "controls_num": 1000,
    }


@pytest.fixture
def mock_term_variants():
    """Test term variants matching our test HDF5 variants."""
    return pd.DataFrame(
        {
            "variant_id": ["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"],
            "gene": ["GENE1", "GENE2", "GENE3"],
            "hmm_score": [0.5, 0.6, 0.7],
            "term": ["TEST:001"] * 3,
        }
    )


def test_read_nomaly_filtered_genotypes(unit_test_app):
    sorted_eids = np.array([1001, 1002, 1003])
    variants = ["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"]

    result = read_nomaly_filtered_genotypes(sorted_eids, variants)

    assert result is not None
    assert np.array_equal(result["row_eids"], sorted_eids)
    assert result["col_variants"] == variants
    assert np.array_equal(result["data"], np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]]))
    assert len(result["error_variants"]) == 0


def test_read_cases_for_disease_code(mock_cases_info):
    with patch(
        "blueprints.prioritisation_by_nomaly_scores.pickle.load",
        return_value=mock_cases_info,
    ):
        cases = read_cases_for_disease_code("571.5")
        assert cases["phecode"] == "571.5"
        assert len(cases["cases"]) == 3
        assert cases["Sex"] == "both"


def test_individual_variant_prioritisation(mock_variant_scores):
    row = np.array([0, 1, 2])  # Example genotypes for one individual
    term_variant_scores = mock_variant_scores

    top_variants = individual_variant_prioritisation(row, term_variant_scores)

    assert isinstance(top_variants, pd.Index)
    assert len(top_variants) > 0
    scores = [term_variant_scores.loc[v].sum() for v in top_variants]
    assert all(scores[i] >= scores[i + 1] for i in range(len(scores) - 1))


def test_term_variant_prioritisation(
    unit_test_app, mock_variant_scores, mock_term_variants
):
    sorted_eids = np.array([1001, 1002, 1003])
    term = "TEST:001"

    with patch(
        "blueprints.prioritisation_by_nomaly_scores.get_term_variants",
        return_value=mock_term_variants,
    ):
        result = term_variant_prioritisation(sorted_eids, mock_variant_scores, term)

        # Test the actual business logic:
        assert isinstance(result, pd.DataFrame)

        # 1. Test that variants are properly filtered based on individual prioritisation
        # The function should select variants where VS scores exceed threshold
        for variant_id in result["variant_id"]:
            scores = mock_variant_scores.loc[variant_id]
            assert any(scores > 1.0), (
                f"Variant {variant_id} included but has no high scores"
            )

        # 2. Test that variants are properly sorted by hmm_score
        assert result["hmm_score"].is_monotonic_decreasing, (
            "Results not sorted by hmm_score"
        )

        # 3. Test that the genotype data was actually used
        # Each returned variant should have appeared in at least one individual's top variants
        genotype_data = unit_test_app.genotype._hdf.genotype_matrix[:]
        assert all(
            any(genotype_data[:, i] > 0)
            for i, variant in enumerate(mock_term_variants["variant_id"])
            if variant in result["variant_id"]
        ), "Returned variants don't match genotype data"

        # 4. Verify required columns are present with correct types
        assert all(col in result.columns for col in ["gene", "variant_id", "hmm_score"])
        assert result["hmm_score"].dtype in (np.float64, float)


def test_get_top_variants(unit_test_app, mock_cases_info, mock_term_variants):
    disease_code = "571.5"
    term = "TEST:001"

    with (
        patch(
            "blueprints.prioritisation_by_nomaly_scores.read_cases_for_disease_code",
            return_value=mock_cases_info,
        ),
        patch(
            "blueprints.prioritisation_by_nomaly_scores.get_term_variants",
            return_value=mock_term_variants,
        ),
    ):
        top_variants, top_gene_set = get_top_variants(disease_code, term)

        assert isinstance(top_variants, pd.DataFrame)
        assert isinstance(top_gene_set, pd.DataFrame)
        assert "gene" in top_variants.columns
        assert "variant_id" in top_variants.columns
        assert "hmm_score" in top_variants.columns
        assert "variant_num" in top_gene_set.columns


def test_stream_logger():
    queue = Queue()
    logger = StreamLogger(queue)

    test_message = "Test progress message"
    logger.info(test_message)

    message = queue.get()
    assert message["type"] == "progress"
    assert message["data"] == test_message
