import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch, MagicMock
from queue import Queue
import json
from pathlib import Path

from blueprints.nomaly import GenotypeHDF5
from blueprints.prioritisation_by_nomaly_scores import (
    read_cases_for_disease_code,
    read_nomaly_filtered_genotypes,
    individual_variant_prioritisation,
    term_variant_prioritisation,
    get_top_variants,
    save_results_to_cache,
    load_cached_results,
    StreamLogger,
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
    with patch("blueprints.nomaly_services.services.genotype", mock_genotype_hdf5):
        yield mock_genotype_hdf5


@pytest.fixture
def mock_config(tmp_path):
    """Create a temporary directory for cache testing"""
    with patch("config.Config") as mock_config:
        mock_config.VARIANT_SCORES_DIR = str(tmp_path)
        yield mock_config


def test_read_cases_for_disease_code(mock_cases_info):
    with patch(
        "blueprints.prioritisation_by_nomaly_scores.pickle.load",
        return_value=mock_cases_info,
    ):
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


def test_individual_variant_prioritisation_edge_cases(mock_variant_scores):
    # Test with -1 values (missing/error data)
    row = np.array([-1, 0, 1])
    result = individual_variant_prioritisation(row, mock_variant_scores)
    assert isinstance(result, pd.Index)

    # Test with all zeros
    row = np.array([0, 0, 0])
    result = individual_variant_prioritisation(row, mock_variant_scores)
    assert len(result) <= 5  # Should only return top 5 when scores are low

    # Test with all valid genotypes
    row = np.array([0, 1, 2])
    result = individual_variant_prioritisation(row, mock_variant_scores)
    assert len(result) > 0


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


def test_term_variant_prioritisation_error_handling(mock_variant_scores):
    sorted_eids = np.array([1001, 1002])
    term = "UP:UPA00240"

    # Test with empty variant list
    mock_db_response = pd.DataFrame(columns=["variant_id", "gene", "hmm_score", "term"])

    with patch(
        "blueprints.prioritisation_by_nomaly_scores.get_term_variants",
        return_value=mock_db_response,
    ):
        result = term_variant_prioritisation(sorted_eids, mock_variant_scores, term)
        assert len(result) == 0


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


def test_get_top_variants_with_cache(mock_config, mock_cases_info, mock_variant_scores):
    disease_code = "571.5"
    term = "UP:UPA00240"

    # First call should compute and cache
    with patch(
        "blueprints.prioritisation_by_nomaly_scores.read_cases_for_disease_code",
        return_value=mock_cases_info,
    ):
        variants1, genes1 = get_top_variants(disease_code, term)

        # Second call should use cache
        variants2, genes2 = get_top_variants(disease_code, term)

        # Results should be identical
        pd.testing.assert_frame_equal(variants1, variants2)
        pd.testing.assert_frame_equal(genes1, genes2)


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


def test_cache_operations(mock_config):
    disease_code = "123"
    term = "TEST:001"

    # Test data
    test_variants = pd.DataFrame(
        {"variant_id": ["1_100_A/G"], "gene": ["GENE1"], "hmm_score": [0.5]}
    )
    test_gene_set = pd.DataFrame(
        {"variant_id": ["1_100_A/G"], "hmm_score": [0.5], "variant_num": [1]},
        index=["GENE1"],
    )

    with patch("blueprints.prioritisation_by_nomaly_scores.Config") as patched_config:
        patched_config.VARIANT_SCORES_DIR = mock_config.VARIANT_SCORES_DIR

        cache_path = (
            Path(patched_config.VARIANT_SCORES_DIR)
            / f"variant_prioritization_{disease_code}_{term}.json"
        )

        # Test saving to cache
        save_results_to_cache(disease_code, term, test_variants, test_gene_set)

        # Debug: Print original data
        print(f"\nOriginal variants:\n{test_variants}")
        print(f"\nOriginal gene_set:\n{test_gene_set}")

        # Verify file was created and print contents
        assert cache_path.exists(), f"Cache file not created at {cache_path}"
        with open(cache_path) as f:
            cached_data = json.load(f)
            print(f"\nCache file contents:\n{json.dumps(cached_data, indent=2)}")
        assert "top_variants" in cached_data
        assert "top_gene_set" in cached_data

        # Test loading from cache (still within the same patch context)
        results = load_cached_results(disease_code, term)
        if results is None:
            print("\nload_cached_results returned None!")
            # Try manual load to debug
            try:
                with open(cache_path) as f:
                    data = json.load(f)
                print("\nManual load of cache file:")
                print(f"Keys in data: {data.keys()}")
                print(f"top_variants type: {type(data['top_variants'])}")
                print(f"top_gene_set type: {type(data['top_gene_set'])}")
                variants_df = pd.DataFrame.from_records(data["top_variants"])
                genes_df = pd.DataFrame.from_records(data["top_gene_set"]).set_index(
                    "gene"
                )
                print(f"\nManually created variants_df:\n{variants_df}")
                print(f"\nManually created genes_df:\n{genes_df}")
            except Exception as e:
                print(f"\nError in manual load: {e}")
        else:
            loaded_variants, loaded_gene_set = results
            print(f"\nLoaded variants:\n{loaded_variants}")
            print(f"\nLoaded gene_set:\n{loaded_gene_set}")

        assert results is not None, "Cache load returned None"
        loaded_variants, loaded_gene_set = results

        # Verify loaded data
        assert not loaded_variants.empty
        assert not loaded_gene_set.empty
        assert loaded_variants["variant_id"].iloc[0] == "1_100_A/G"


def test_stream_logger():
    queue = Queue()
    logger = StreamLogger(queue)

    test_message = "Test progress message"
    logger.info(test_message)

    message = queue.get()
    assert message["type"] == "progress"
    assert message["data"] == test_message
