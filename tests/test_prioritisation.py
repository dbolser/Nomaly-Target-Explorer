import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch, MagicMock

from blueprints.nomaly import GenotypeHDF5
from blueprints.prioritisation_by_nomaly_scores import (
    read_cases_for_disease_code,
    read_nomaly_filtered_genotypes,
    individual_variant_prioritisation,
    term_variant_prioritisation,
    get_top_variants,
)


# Mock data fixtures
@pytest.fixture
def mock_variant_scores():
    return pd.DataFrame(
        {"VS00": [0.1, 0.2, 0.3], "VS01": [0.4, 0.5, 0.6], "VS11": [0.7, 0.8, 0.9]},
        index=["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"],
    )


@pytest.fixture
def mock_cases_info():
    return {
        "phecode": "571.5",
        "Sex": "both",
        "phecode_exclude_range": [],
        "cases": [1001, 1002, 1003],
        "exclude": [],
        "controls_num": 1000,
    }


@pytest.fixture
def mock_genotypes_data():
    return {
        "row_eids": np.array([1001, 1002, 1003]),
        "col_variants": ["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"],
        "data": np.array(
            [
                [0, 1, 2],  # Individual 1
                [1, 0, 1],  # Individual 2
                [2, 1, 0],  # Individual 3
            ]
        ),
        "error_variants": [],
    }


@pytest.fixture
def mock_genotype_hdf5():
    """Create a mock GenotypeHDF5 instance."""
    mock = MagicMock(spec=GenotypeHDF5)
    mock.individual = np.array([1001, 1002, 1003])
    mock.query_variantID_genotypes.side_effect = lambda x: (
        np.array([1001, 1002, 1003]),
        np.array([0, 1, 2]),
    )
    return mock


@pytest.fixture(autouse=True)
def patch_nomaly_genotype(mock_genotype_hdf5):
    """Patch the nomaly_genotype instance in prioritisation_by_nomaly_scores."""
    with patch(
        "blueprints.prioritisation_by_nomaly_scores.nomaly_genotype", mock_genotype_hdf5
    ):
        yield mock_genotype_hdf5


def test_read_cases_for_disease_code(mock_cases_info):
    with patch("pickle.load", return_value=mock_cases_info):
        cases = read_cases_for_disease_code("571.5")
        assert cases["phecode"] == "571.5"
        assert len(cases["cases"]) == 3
        assert cases["Sex"] == "both"


def test_read_nomaly_filtered_genotypes(mock_genotypes_data, mock_genotype_hdf5):
    sorted_eids = np.array([1001, 1002, 1003])
    variants = ["1_186977737_A/G", "1_186977780_G/A", "1_46813503_C/T"]

    result = read_nomaly_filtered_genotypes(sorted_eids, variants)

    assert (result["row_eids"] == mock_genotypes_data["row_eids"]).all()
    assert len(result["col_variants"]) == len(variants)
    assert result["data"].shape == (len(sorted_eids), len(variants))
    assert len(result["error_variants"]) == 0


def test_individual_variant_prioritisation(mock_variant_scores):
    row = np.array([0, 1, 2])  # Example genotypes for one individual
    term_variant_scores = mock_variant_scores

    top_variants = individual_variant_prioritisation(row, term_variant_scores)

    assert isinstance(top_variants, pd.Index)
    assert len(top_variants) > 0
    # Verify the variants are sorted by their scores
    scores = [term_variant_scores.loc[v].sum() for v in top_variants]
    assert all(scores[i] >= scores[i + 1] for i in range(len(scores) - 1))


def test_term_variant_prioritisation(mock_variant_scores, mock_genotypes_data):
    sorted_eids = mock_genotypes_data["row_eids"]
    term = "UP:UPA00240"

    # Create mock DB response first
    mock_db_response = pd.DataFrame(
        {
            "variant_id": [
                "1_186977737_A/G",
                "1_186977780_G/A",
                "1_46813503_C/T",
                "1_46815113_G/A",
                "1_46815158_A/G",
                "1_46815212_C/T",
                "1_46817100_C/T",
                "1_46930161_C/T",
                "1_46930245_C/T",
                "1_46932824_A/G",
            ],
            "gene": [
                "GENE1",
                "GENE2",
                "GENE3",
                "GENE4",
                "GENE5",
                "GENE6",
                "GENE7",
                "GENE8",
                "GENE9",
                "GENE10",
            ],
            "hmm_score": [0.5, 0.6, 0.7, 0.5, 0.6, 0.7, 0.5, 0.6, 0.7, 0.5],
            "term": ["UP:UPA00240"] * 10,
        }
    )

    # Create mock scores for all variants in the DB response
    mock_scores = pd.DataFrame(
        {
            "VS00": np.random.uniform(0.1, 0.3, len(mock_db_response)),
            "VS01": np.random.uniform(0.4, 0.6, len(mock_db_response)),
            "VS11": np.random.uniform(0.7, 0.9, len(mock_db_response)),
        },
        index=mock_db_response["variant_id"],
    )

    with (
        patch("db.get_term_variants", return_value=mock_db_response),
        patch(
            "blueprints.prioritisation_by_nomaly_scores.read_nomaly_filtered_genotypes",
            return_value={
                "row_eids": sorted_eids,
                "col_variants": mock_db_response["variant_id"].tolist(),
                "data": np.random.randint(
                    0, 3, (len(sorted_eids), len(mock_db_response))
                ),
                "error_variants": [],
            },
        ),
    ):
        result = term_variant_prioritisation(sorted_eids, mock_scores, term)

        assert isinstance(result, pd.DataFrame)
        assert "gene" in result.columns
        assert "variant_id" in result.columns
        assert "hmm_score" in result.columns
        assert len(result) > 0
        # Check that the variants in the result are from our mock data
        assert all(
            v in mock_db_response["variant_id"].tolist() for v in result["variant_id"]
        )


def test_get_top_variants(mock_cases_info, mock_variant_scores):
    disease_code = "571.5"
    term = "UP:UPA00240"

    with (
        patch(
            "blueprints.prioritisation_by_nomaly_scores.read_cases_for_disease_code",
            return_value=mock_cases_info,
        ),
        patch(
            "blueprints.prioritisation_by_nomaly_scores.term_variant_prioritisation",
            return_value=pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2"],
                    "variant_id": ["1_186977737_A/G", "1_186977780_G/A"],
                    "hmm_score": [0.5, 0.6],
                }
            ),
        ),
    ):
        top_variants, top_gene_set = get_top_variants(disease_code, term)

        assert isinstance(top_variants, pd.DataFrame)
        assert isinstance(top_gene_set, pd.DataFrame)
        assert "gene" in top_variants.columns
        assert "variant_id" in top_variants.columns
        assert "hmm_score" in top_variants.columns
        assert "variant_num" in top_gene_set.columns


def test_error_handling_missing_variants(mock_genotype_hdf5):
    sorted_eids = np.array([1001, 1002, 1003])
    variants = ["nonexistent1", "nonexistent2"]

    # Configure mock to return None for nonexistent variants
    mock_genotype_hdf5.query_variantID_genotypes.side_effect = lambda x: None

    result = read_nomaly_filtered_genotypes(sorted_eids, variants)
    assert len(result["error_variants"]) == len(variants)


def test_empty_case_handling():
    disease_code = "nonexistent"
    term = "UP:UPA00240"

    with patch(
        "blueprints.prioritisation_by_nomaly_scores.read_cases_for_disease_code",
        return_value={"cases": []},
    ):
        top_variants, top_gene_set = get_top_variants(disease_code, term)
        assert len(top_variants) == 0
        assert len(top_gene_set) == 0


def test_large_dataset_performance(mock_genotype_hdf5):
    # Test with a larger dataset to ensure performance is acceptable
    num_individuals = 1000
    num_variants = 100

    sorted_eids = np.arange(1001, 1001 + num_individuals)
    variants = [f"1_{i}_A/G" for i in range(num_variants)]

    # Configure mock for larger dataset
    mock_genotype_hdf5.individual = sorted_eids
    mock_genotype_hdf5.query_variantID_genotypes.side_effect = lambda x: (
        sorted_eids,
        np.random.randint(0, 3, num_individuals),
    )

    result = read_nomaly_filtered_genotypes(sorted_eids, variants)
    assert result["data"].shape == (num_individuals, num_variants)
