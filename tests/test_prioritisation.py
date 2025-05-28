import pytest
import numpy as np
import pandas as pd
from queue import Queue

from blueprints.prioritisation_by_nomaly_scores import (
    StreamLogger,
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


def test_fetch_nomaly_scores(mock_nomaly_scores_service):
    """Test fetching Nomaly scores."""
    case_eids = np.array([1001, 1002, 1003])
    control_eids = np.array([1004, 1005, 1006, 1007, 1008])
    term = "TEST:001"

    result_case_scores = mock_nomaly_scores_service.get_scores_by_eids_unsorted(
        case_eids, term
    )
    result_control_scores = mock_nomaly_scores_service.get_scores_by_eids_unsorted(
        control_eids, term
    )

    assert len(result_case_scores) == len(case_eids)
    assert len(result_control_scores) == len(control_eids)
    assert np.array_equal(result_case_scores, np.array([[0.001], [0.002], [0.003]]))
    assert np.array_equal(
        result_control_scores, np.array([[0.004], [0.005], [0.006], [0.007], [0.008]])
    )


def test_fetch_stats_data(mock_stats_service):
    """Test fetching statistics data."""
    term = "TEST:001"
    phecode = "571.5"

    stats = mock_stats_service.get_term_stats(term, phecode=phecode)

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


# You might want to add more tests covering edge cases:
# - No individuals
# - No variants
# - All genotypes missing
# - Scores are all zero or negative
