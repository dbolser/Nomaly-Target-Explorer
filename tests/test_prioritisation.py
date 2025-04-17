import pytest
import numpy as np
import pandas as pd
from queue import Queue
from pandas.testing import assert_frame_equal

from blueprints.prioritisation_by_nomaly_scores import (
    StreamLogger,
    process_individual_variants,
    process_individual_variants_vectorized,
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


def test_fetch_nomaly_scores(nomaly_scores_service):
    """Test fetching Nomaly scores."""
    case_eids = np.array([1001, 1002, 1003])
    control_eids = np.array([1004, 1005, 1006, 1007, 1008])
    term = "TEST:001"

    result_case_scores = nomaly_scores_service.get_scores_by_eids_unsorted(
        case_eids, term
    )
    result_control_scores = nomaly_scores_service.get_scores_by_eids_unsorted(
        control_eids, term
    )

    assert len(result_case_scores) == len(case_eids)
    assert len(result_control_scores) == len(control_eids)
    assert np.array_equal(result_case_scores, np.array([0.001, 0.002, 0.003])), (
        "Case scores should be [0.001, 0.002, 0.003]",
        "instead we got",
        result_case_scores,
    )
    assert np.array_equal(
        result_control_scores, np.array([0.004, 0.005, 0.006, 0.007, 0.008])
    ), "Control scores should be [0.004, 0.005, 0.006, 0.007, 0.008]"


def test_fetch_stats_data(stats_service):
    """Test fetching statistics data."""
    term = "TEST:001"
    phecode = "571.5"

    stats = stats_service.get_term_stats(term, phecode=phecode)

    assert len(stats) == 1
    stats = stats.iloc[0]

    assert stats["num_rp"] == 3.0
    assert stats["num_rn"] == 5.0
    assert stats["metric1_pvalue"] == 0.0026
    assert stats["roc_stats_mcc_pvalue"] == 0.0015


def test_stream_logger():
    """Test the StreamLogger class."""
    queue = Queue()
    logger = StreamLogger(queue)

    test_message = "Test progress message"
    logger.log(test_message)

    message = queue.get()
    assert message["type"] == "progress"
    assert message["data"] == test_message


def test_process_individual_variants_equivalence():
    """Verify that vectorized and non-vectorized functions produce identical results."""
    # 1. Create sample data
    variants = [f"v{i}" for i in range(1, 6)]  # v1, v2, v3, v4, v5
    term_variant_scores = pd.DataFrame(
        {
            "variant_id": variants,
            "vs00": [0.1, 0.2, 0.0, 0.5, 0.1],
            "vs01": [1.1, 1.5, 0.5, 2.0, 0.6],
            "vs11": [2.5, 3.0, 1.2, 4.0, 1.3],
        }
    )

    # Genotypes: (num_variants, num_individuals)
    # Individuals: Ind1, Ind2, Ind3, Ind4
    genotypes = np.array(
        [
            [0, 1, 2, 0],  # v1 genotypes for Ind1, Ind2, Ind3, Ind4
            [1, 1, 0, 2],  # v2
            [2, 0, -9, 1],  # v3 (-9 is missing)
            [0, 1, 1, 0],  # v4
            [1, 2, 0, 1],  # v5
        ]
    )

    # 2. Run both functions
    expected_df = process_individual_variants(
        genotypes.copy(), term_variant_scores.copy()
    )
    actual_df = process_individual_variants_vectorized(
        genotypes.copy(), term_variant_scores.copy()
    )

    # 3. Compare results
    # Both functions now sort internally by 'vs' descending.
    # Reset index for consistent comparison as order might be the only difference otherwise.
    expected_df = expected_df.reset_index(drop=True)
    actual_df = actual_df.reset_index(drop=True)

    # Use pandas testing function for robust comparison
    assert_frame_equal(actual_df, expected_df, check_dtype=True)


# You might want to add more tests covering edge cases:
# - No individuals
# - No variants
# - All genotypes missing
# - Scores are all zero or negative
