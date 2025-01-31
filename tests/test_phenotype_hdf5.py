import pytest
import numpy as np
import h5py
import tempfile

from data_services.implementations.hdf5.phenotype import HDF5PhenotypeService
from data_services.interfaces.phenotype import PhenotypeService


@pytest.fixture
def mock_phenotype_file():
    """Create a temporary mock HDF5 file with test phenotype data."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create required datasets
            # Phenotype matrix: rows are phecodes, columns are individuals
            # 1 = case, 0 = control, -1 = excluded
            f.create_dataset(
                "phenotype_matrix",
                data=np.array(
                    [
                        [-1, 0, 1, 0],  # First phecode
                        [1, 0, -1, 1],  # Second phecode
                    ],
                    dtype=np.int8,
                ),
            )
            # Individual IDs
            f.create_dataset("individual", data=np.array([1001, 1002, 1003, 1004]))
            # Phecode labels
            f.create_dataset(
                "phecode",
                data=np.array([b"250.2", b"401.1"]),  # Diabetes and Hypertension
            )
        yield tmp.name


@pytest.fixture
def phenotype_service(mock_phenotype_file):
    """Create a PhenotypeService instance with the mock data."""
    return HDF5PhenotypeService(mock_phenotype_file)


def test_service_initialization(phenotype_service):
    """Test that the service initializes correctly."""
    assert isinstance(phenotype_service, PhenotypeService)
    assert hasattr(phenotype_service, "_hdf")


def test_get_cases_basic(phenotype_service):
    """Test basic case retrieval for a phecode."""
    eids, case_status = phenotype_service.get_cases_for_phecode("250.2")

    assert isinstance(eids, np.ndarray)
    assert isinstance(case_status, np.ndarray)
    assert len(eids) == len(case_status)

    # For phecode 250.2 in our mock data:
    # Individual 1001: excluded (-1)
    # Individual 1002: control (0)
    # Individual 1003: case (1)
    # Individual 1004: control (0)
    expected_eids = np.array([1002, 1003, 1004])
    expected_status = np.array([0, 1, 0])

    np.testing.assert_array_equal(eids, expected_eids)
    np.testing.assert_array_equal(case_status, expected_status)


def test_get_cases_all_excluded(phenotype_service):
    """Test handling of a theoretical phecode where all individuals are excluded."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            f.create_dataset(
                "phenotype_matrix",
                data=np.array([[-1, -1, -1, -1]], dtype=np.int8),
            )
            f.create_dataset("individual", data=np.array([1001, 1002, 1003, 1004]))
            f.create_dataset("phecode", data=np.array([b"999.9"]))

        service = HDF5PhenotypeService(tmp.name)
        eids, case_status = service.get_cases_for_phecode("999.9")

        assert len(eids) == 0
        assert len(case_status) == 0


def test_get_cases_all_controls(phenotype_service):
    """Test handling of a phecode with only controls (no cases)."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            f.create_dataset(
                "phenotype_matrix",
                data=np.array([[0, 0, 0, 0]], dtype=np.int8),
            )
            f.create_dataset("individual", data=np.array([1001, 1002, 1003, 1004]))
            f.create_dataset("phecode", data=np.array([b"888.8"]))

        service = HDF5PhenotypeService(tmp.name)
        eids, case_status = service.get_cases_for_phecode("888.8")

        assert len(eids) == 4
        assert len(case_status) == 4
        assert all(status == 0 for status in case_status)


def test_get_cases_nonexistent_phecode(phenotype_service):
    """Test handling of a nonexistent phecode."""
    with pytest.raises(ValueError):
        phenotype_service.get_cases_for_phecode("999.999")


def test_get_cases_invalid_data(phenotype_service):
    """Test handling of invalid data in the HDF5 file."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create mismatched data
            f.create_dataset(
                "phenotype_matrix",
                data=np.array([[0, 0]], dtype=np.int8),  # Only 2 individuals
            )
            f.create_dataset(
                "individual",
                data=np.array([1001, 1002, 1003]),  # 3 individuals
            )
            f.create_dataset("phecode", data=np.array([b"777.7"]))

        service = HDF5PhenotypeService(tmp.name)
        with pytest.raises(AssertionError):
            service.get_cases_for_phecode("777.7")


def test_get_cases_with_population(phenotype_service):
    """Test case retrieval with population filter."""
    # This is a placeholder test - implement once population filtering is added
    eids, case_status = phenotype_service.get_cases_for_phecode(
        "250.2", population="EUR"
    )
    # For now, population parameter is ignored, so results should be the same
    expected_eids = np.array([1002, 1003, 1004])
    expected_status = np.array([0, 1, 0])

    np.testing.assert_array_equal(eids, expected_eids)
    np.testing.assert_array_equal(case_status, expected_status)


def test_get_cases_with_sex(phenotype_service):
    """Test case retrieval with sex filter."""
    # This is a placeholder test - implement once sex filtering is added
    eids, case_status = phenotype_service.get_cases_for_phecode("250.2", sex="female")
    # For now, sex parameter is ignored, so results should be the same
    expected_eids = np.array([1002, 1003, 1004])
    expected_status = np.array([0, 1, 0])

    np.testing.assert_array_equal(eids, expected_eids)
    np.testing.assert_array_equal(case_status, expected_status)
