import pytest
from pathlib import Path
import numpy as np
import os

from config import Config
from data_services.nomaly_score import NomalyScoreService

# === Configuration & Constants ===
# Path to the test HDF5 file - this should be properly configured
TEST_HDF5_FILE = Path("/data/analysis/UKBB/Run-v1/DatabaseInputs/float32_scores.h5")

# Test data constants
TEST_EID = 2103010
TEST_SECOND_EID = 5243997
TEST_TERMS = np.array(["GO:0043551", "GO:0018146", "GO:0030176", "GO:0031109"])

# Expected scores for test data - UPDATE WITH ACTUAL VALUES
EXPECTED_SCORES_2103010 = {
    "GO:0043551": 1.3262272e-02,
    "GO:0018146": 4.2215530e-07,
    "GO:0030176": 6.9814450e-03,
    "GO:0031109": 2.1157183e-08,
}

EXPECTED_SCORES_5243997 = {
    "GO:0043551": 0.0115555042459608,
    "GO:0018146": 5.86770890964509e-17,
    "GO:0030176": 0.0027654172916,  # Update this with actual value
    "GO:0031109": 1.3377684e-10,  # Update this with actual value
}


def check_file_exists():
    """Check if the test HDF5 file exists, skip module if not."""
    if not TEST_HDF5_FILE.exists():
        pytest.skip(
            f"Test HDF5 file not found: {TEST_HDF5_FILE}", allow_module_level=True
        )


# Check at module load time
check_file_exists()


@pytest.fixture(scope="module")
def nomaly_service():
    """Fixture to provide a configured NomalyScoreService instance."""
    service = NomalyScoreService(TEST_HDF5_FILE)
    assert service.initialized, "Service failed to initialize"
    return service


# === Basic Service Tests ===


def test_service_initialization():
    """Test that the service initializes correctly with a valid file."""
    service = NomalyScoreService(TEST_HDF5_FILE)
    assert service.initialized
    assert hasattr(service, "_hdf")
    assert service._hdf.f is not None


def test_service_initialization_none():
    """Test that the service initializes correctly with None."""
    service = NomalyScoreService(None)
    assert not service.initialized
    with pytest.raises(ValueError, match="Service not properly initialized"):
        service.get_scores_by_eids_unsorted(np.array([TEST_EID]))


def test_service_properties(nomaly_service):
    """Test basic properties of the service."""
    # Check that the HDF5 file is properly loaded
    assert nomaly_service._hdf.f is not None

    # Check that we have the expected arrays
    assert hasattr(nomaly_service._hdf, "eids")
    assert hasattr(nomaly_service._hdf, "terms")
    assert hasattr(nomaly_service._hdf, "data_matrix")

    # Check arrays are non-empty
    assert len(nomaly_service._hdf.eids) > 0
    assert len(nomaly_service._hdf.terms) > 0


# === Individual Score Tests ===


def test_get_scores_single_eid_all_terms(nomaly_service):
    """Test retrieving scores for a single EID and all terms."""
    eids = np.array([TEST_EID])
    scores = nomaly_service.get_scores_by_eids_unsorted(eids)

    assert scores is not None
    assert isinstance(scores, np.ndarray)
    assert scores.ndim == 2, f"Expected 2D array, got shape {scores.shape}"
    assert scores.shape[0] == 1, f"Expected 1 row, got {scores.shape[0]}"
    assert scores.shape[1] == len(nomaly_service._hdf.terms), (
        f"Expected {len(nomaly_service._hdf.terms)} columns"
    )


def test_get_scores_single_eid_specific_terms(nomaly_service):
    """Test retrieving scores for a single EID and specific terms."""
    eids = np.array([TEST_EID])
    scores = nomaly_service.get_scores_by_eids_unsorted(eids, TEST_TERMS)

    assert scores is not None
    assert isinstance(scores, np.ndarray)
    assert scores.ndim == 2, f"Expected 2D array, got shape {scores.shape}"
    assert scores.shape == (1, len(TEST_TERMS)), (
        f"Expected shape (1, {len(TEST_TERMS)}), got {scores.shape}"
    )

    # Check specific known values
    for i, term in enumerate(TEST_TERMS):
        if term in EXPECTED_SCORES_2103010:
            expected = EXPECTED_SCORES_2103010[term]
            actual = scores[0, i]
            assert np.isclose(actual, expected, rtol=1e-6, atol=1e-9), (
                f"Score mismatch for {term}: expected {expected}, got {actual}"
            )


def test_get_scores_multiple_eids_specific_terms(nomaly_service):
    """Test retrieving scores for multiple EIDs and specific terms."""
    eids = np.array([TEST_EID, TEST_SECOND_EID])
    terms = np.array(["GO:0043551", "GO:0018146"])  # Subset of TEST_TERMS

    scores = nomaly_service.get_scores_by_eids_unsorted(eids, terms)

    assert scores is not None
    assert isinstance(scores, np.ndarray)
    assert scores.ndim == 2, f"Expected 2D array, got shape {scores.shape}"
    assert scores.shape == (len(eids), len(terms)), (
        f"Expected shape {(len(eids), len(terms))}, got {scores.shape}"
    )

    # Create expected array based on our constants
    expected = np.array(
        [
            [
                EXPECTED_SCORES_2103010["GO:0043551"],
                EXPECTED_SCORES_2103010["GO:0018146"],
            ],
            [
                EXPECTED_SCORES_5243997["GO:0043551"],
                EXPECTED_SCORES_5243997["GO:0018146"],
            ],
        ],
        dtype=np.float32,
    )

    # Check the entire array matches expected values
    assert np.allclose(scores, expected, rtol=1e-6, atol=1e-9), (
        f"Scores mismatch.\nExpected:\n{expected}\nActual:\n{scores}"
    )


# === Edge Cases and Error Handling ===


def test_nonexistent_eid(nomaly_service):
    """Test behavior when requesting a non-existent EID."""
    # This test might need adjustment based on actual implementation
    # If the code should raise an error for non-existent EIDs, use pytest.raises
    # If it should return zeros or NaNs, check for that pattern
    non_existent_eid = 99999999  # Choose an EID that's definitely not in the data

    # The current implementation in NomalyScoreHDF5.get_scores_by_eids_unsorted
    # has an assertion that all EIDs must exist, so we expect it to raise
    with pytest.raises(AssertionError):
        nomaly_service.get_scores_by_eids_unsorted(np.array([non_existent_eid]))


def test_empty_eids_array(nomaly_service):
    """Test behavior with an empty EIDs array."""
    # This might raise an error or return an empty array depending on implementation
    try:
        result = nomaly_service.get_scores_by_eids_unsorted(np.array([]))
        assert result.shape[0] == 0, (
            f"Expected empty first dimension, got shape {result.shape}"
        )
    except (ValueError, AssertionError, IndexError):
        # If implementation raises an error for empty input, that's also acceptable
        pass


def test_nonexistent_terms(nomaly_service):
    """Test behavior when requesting non-existent terms."""
    non_existent_terms = np.array(["FAKE:9999", "NONEXISTENT:1234"])

    # This might either filter out non-existent terms or raise an error
    try:
        result = nomaly_service.get_scores_by_eids_unsorted(
            np.array([TEST_EID]), non_existent_terms
        )
        # If it returns a result, it should have the right shape but might have zeros
        assert result.shape == (1, 0), f"Expected shape (1, 0), got {result.shape}"
    except (ValueError, IndexError):
        # If implementation raises an error for missing terms, that's also acceptable
        pass


# === Performance and Data Consistency Tests ===


def test_data_consistency_multiple_eids(nomaly_service):
    """Test that data for the same EID is consistent across different query patterns."""
    # Get scores for a single EID
    single_result = nomaly_service.get_scores_by_eids_unsorted(
        np.array([TEST_EID]), TEST_TERMS
    )

    # Get scores for multiple EIDs including the test EID
    multi_result = nomaly_service.get_scores_by_eids_unsorted(
        np.array([TEST_EID, TEST_SECOND_EID]), TEST_TERMS
    )

    # The first row of multi_result should exactly match single_result
    assert np.array_equal(single_result[0], multi_result[0]), (
        "Inconsistent results for same EID in different queries"
    )


def test_scores_data_ranges(nomaly_service):
    """Test that scores are within expected ranges."""
    # Nomaly scores should typically be between 0 and 1
    scores = nomaly_service.get_scores_by_eids_unsorted(
        np.array([TEST_EID, TEST_SECOND_EID])
    )

    # Scores should be finite (not NaN or Inf)
    assert np.all(np.isfinite(scores)), "Found non-finite scores (NaN or Inf)"

    # Most scores should be between 0 and 1, but very small values close to 0 are common
    # This is just a sanity check - adjust based on your actual data characteristics
    assert np.all((scores >= 0) & (scores <= 1)), "Found scores outside range [0, 1]"


# === Original Test (Fixed) ===

def test_get_scores_by_eids_unsorted_specific(nomaly_service):
    """
    Tests get_scores_by_eids_unsorted with a specific EID and set of terms
    using known correct values.
    """
    eids_to_test = np.array([2103010])  # Use numpy array as input type hint suggests
    terms_to_test = np.array(
        ["GO:0043551", "GO:0018146", "GO:0046173", "GO:0030176", "GO:0031109"]
    )

    # Expected values based on your manual verification
    expected_scores = np.array(
        [[1.3262272e-02, 4.2215530e-07, 7.5833222e-22, 6.9814450e-03, 2.1157183e-08]],
        dtype=np.float32,
    )  # Match the expected dtype, and make sure it's 2D with shape (1, 5)

    # Call the method under test
    actual_scores = nomaly_service.get_scores_by_eids_unsorted(
        eids_to_test, terms_to_test
    )

    # Assertions
    assert actual_scores is not None, "Method returned None"
    assert isinstance(actual_scores, np.ndarray), (
        f"Expected numpy array, got {type(actual_scores)}"
    )

    # Check shape - should be (num_eids, num_terms) -> (1, 5)
    assert actual_scores.shape == (len(eids_to_test), len(terms_to_test)), (
        f"Expected shape {(len(eids_to_test), len(terms_to_test))}, got {actual_scores.shape}"
    )

    # Check the actual values using np.allclose for floating-point comparisons
    assert np.allclose(actual_scores, expected_scores, rtol=1e-6, atol=1e-25), (
        f"Scores mismatch.\nExpected: {expected_scores}\nActual: {actual_scores}"
    )


def test_get_scores_by_eids_unsorted_multiple_eids(nomaly_service):
    """Test with multiple EIDs"""
    eids_to_test = np.array([2103010, 5243997])
    terms_to_test = np.array(["GO:0043551", "GO:0018146"])

    # Expected values - these should match what's in your constants
    expected_scores = np.array(
        [
            [
                EXPECTED_SCORES_2103010["GO:0043551"],
                EXPECTED_SCORES_2103010["GO:0018146"],
            ],  # First EID
            [
                EXPECTED_SCORES_5243997["GO:0043551"],
                EXPECTED_SCORES_5243997["GO:0018146"],
            ],  # Second EID
        ],
        dtype=np.float32,
    )

    actual_scores = nomaly_service.get_scores_by_eids_unsorted(
        eids_to_test, terms_to_test
    )

    # Check shape
    assert actual_scores.shape == (len(eids_to_test), len(terms_to_test)), (
        f"Shape mismatch: expected {(len(eids_to_test), len(terms_to_test))}, got {actual_scores.shape}"
    )

    # Check values
    assert np.allclose(actual_scores, expected_scores, rtol=1e-6, atol=1e-9), (
        f"Scores mismatch.\nExpected:\n{expected_scores}\nActual:\n{actual_scores}"
    )
