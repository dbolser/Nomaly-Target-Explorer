import pytest
from pathlib import Path
import numpy as np

# Assuming your classes are in data_services/nomaly_score.py
# Adjust the import path as needed
from data_services.nomaly_score import NomalyScoreService

# Define the path to your test HDF5 file
# IMPORTANT: Replace this with the actual path to your HDF5 file for testing
TEST_HDF5_FILE = Path("/data/general/UKBB/Run-v1/DatabaseInputs/float32_scores.h5")

# Check if the test file exists, skip if not
if not TEST_HDF5_FILE.exists():
    pytest.skip(f"Test HDF5 file not found: {TEST_HDF5_FILE}", allow_module_level=True)


@pytest.fixture(scope="module")
def nomaly_service() -> NomalyScoreService:
    """Fixture to initialize the NomalyScoreService once per module."""
    print(f"\nInitializing NomalyScoreService with {TEST_HDF5_FILE}")
    service = NomalyScoreService(TEST_HDF5_FILE)
    assert service.initialized, "Service failed to initialize"
    print("NomalyScoreService initialized.")
    # You might want to add cleanup here if needed, e.g., closing the HDF5 file
    # though h5py might handle it if the service instance goes out of scope.
    # yield service
    # print("Cleaning up NomalyScoreService")
    # if service._hdf and service._hdf.f:
    #     service._hdf.f.close()
    return service  # Simpler for now if explicit close isn't strictly needed


def test_get_scores_by_eids_unsorted_specific(nomaly_service: NomalyScoreService):
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
        [1.3262272e-02, 4.2215530e-07, 7.5833222e-22, 6.9814450e-03, 2.1157183e-08],
        dtype=np.float32,
    )  # Match the expected dtype

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
    # assert actual_scores.shape == (len(eids_to_test), len(terms_to_test)), (
    #     f"Expected shape {(len(eids_to_test), len(terms_to_test))}, got {actual_scores.shape}"
    # )

    # Check the actual values using np.allclose for floating-point comparisons
    # Adjust tolerance (atol, rtol) if needed based on data precision
    assert np.allclose(
        actual_scores.flatten(), expected_scores, rtol=1e-6, atol=1e-25
    ), (
        f"Scores mismatch.\nExpected: {expected_scores}\nActual:   {actual_scores.flatten()}"
    )


def test_get_scores_by_eids_unsorted_multiple_eids(nomaly_service: NomalyScoreService):
    """Test with multiple EIDs"""
    # You'll need to find known correct values for these EIDs/Terms first
    eids_to_test = np.array(
        [2103010, 5243997]
    )  # Example: first two eids from your h5ls output
    terms_to_test = np.array(["GO:0043551", "GO:0018146"])  # Example: first two terms

    # --- IMPORTANT: Replace these with actual expected values ---
    # You would get these by querying the HDF5 file or source data directly
    # For EID 2103010: [0.0132623, 4.22155e-07]
    # For EID 5243997: You need to find these values! Let's use placeholders.
    # Example: Suppose the values for 5243997 are 0.1 and 0.2
    expected_scores = np.array(
        [
            [1.3262272e-02, 4.2215530e-07],  # Row for EID 2103010
            [
                0.0115555042459608,
                5.86770890964509e-17,
            ],  #
        ],
        dtype=np.float32,
    )
    # -----------------------------------------------------------

    actual_scores = nomaly_service.get_scores_by_eids_unsorted(
        eids_to_test, terms_to_test
    )

    # assert actual_scores.shape == (len(eids_to_test), len(terms_to_test)), (
    #     "Shape mismatch"
    # )

    # Replace the placeholder assertion once you have the real expected values
    assert np.allclose(actual_scores, expected_scores, rtol=1e-6, atol=1e-9), (
        f"Scores mismatch.\nExpected:\n{expected_scores}\nActual:\n{actual_scores}"
    )


# Add more tests:
# - Test with terms=None (all terms)
# - Test with EIDs not found (should raise ValueError based on updated code)
# - Test with terms not found (should print warning or raise, depending on implementation)
# - Test edge cases (empty EID list, empty term list?)
