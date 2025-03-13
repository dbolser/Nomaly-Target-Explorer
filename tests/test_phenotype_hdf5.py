"""This is an DATA TEST, where we look at the real production data."""

import pytest
import numpy as np
import h5py
import tempfile

from data_services.phenotype import PhenotypeService

from config import Config


# By setting the scope to "module", we only read the HDF5 file once for all the
# following tests.. I think.
@pytest.fixture(scope="module")
def production_phenotype_service():
    """Fixture to create a PhenotypeService instance."""
    return PhenotypeService(hdf5_file=Config.PHENOTYPES_H5)


def test_get_cases_basic(production_phenotype_service):
    """Test basic case retrieval for a phecode."""
    eids, case_status = production_phenotype_service.get_cases_for_phecode("250.2")

    assert isinstance(eids, np.ndarray)
    assert isinstance(case_status, np.ndarray)
    assert len(eids) == len(case_status)

    # Just a sample...
    expected_eids = np.array([1000015, 1000027, 1000032, 6022427, 6022432, 6022446])

    np.testing.assert_array_equal(eids[:3], expected_eids[:3])
    np.testing.assert_array_equal(eids[-3:], expected_eids[-3:])

    assert np.sum(case_status == 1) == 42973


def test_get_cases_nonexistent_phecode(production_phenotype_service):
    """Test handling of a nonexistent phecode."""
    with pytest.raises(KeyError):
        production_phenotype_service.get_cases_for_phecode("999.999")


# TODO: Move this to a test for the HDF5 file itself
def test_get_cases_invalid_data():
    """Test handling of invalid data in the HDF5 file."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create mismatched data
            f.create_dataset(
                "phenotype_data",
                data=np.array([[0], [0]], dtype=np.int8),  # Only 2 individuals
            )
            f.create_dataset("eids", data=np.array([101001, 101002, 101003]))
            f.create_dataset("phecodes", data=np.array([b"777.7"]))
            f.create_dataset("sex", data=np.array([b"M", b"F", b"M"]))
            f.create_dataset("populations", data=np.array([b"EUR", b"EUR", b"SAS"]))

            with pytest.raises(KeyError):
                _ = PhenotypeService(tmp.name)


def test_get_cases_with_population(production_phenotype_service):
    """Test case retrieval with population filter."""
    # This is a placeholder test - implement once population filtering is added
    eids, case_status = production_phenotype_service.get_cases_for_phecode(
        "250.2", population="SAS"
    )
    # For now, population parameter is ignored, so results should be the same
    expected_eids = np.array([1001013, 1002039, 1002389, 6021272, 6021752, 6021889])

    np.testing.assert_array_equal(eids[:3], expected_eids[:3])
    np.testing.assert_array_equal(eids[-3:], expected_eids[-3:])
    assert np.sum(case_status == 1) == 2261


def test_get_cases_with_sex(production_phenotype_service):
    """Test case retrieval with sex filter."""
    # This is a placeholder test - implement once sex filtering is added
    eids_with_sex_f, case_status_with_sex_f = (
        production_phenotype_service.get_cases_for_phecode("250.2", biological_sex="F")
    )
    eids_with_sex_m, case_status_with_sex_m = (
        production_phenotype_service.get_cases_for_phecode("250.2", biological_sex="M")
    )

    assert np.sum(case_status_with_sex_f == 1) == 17434
    assert np.sum(case_status_with_sex_m == 1) == 25539

    # Test that it's a strict of the cases without sex filter
    eids, case_status = production_phenotype_service.get_cases_for_phecode("250.2")

    assert np.all(np.union1d(eids_with_sex_f, eids_with_sex_m) == eids)
