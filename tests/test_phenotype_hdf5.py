import pytest
import numpy as np
import h5py
import tempfile

from data_services.phenotype import PhenotypeService


@pytest.fixture
def mock_phenotype_file():
    """Create a temporary mock HDF5 file with test phenotype data."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create required datasets

            # Phenotype matrix: columns are eids rows are phecodes,
            # 1 = case, 0 = control, -1 = excluded
            f.create_dataset(
                "phenotype_data",
                data=np.array(
                    [
                        [9, 1, 0, 1],  # M
                        [0, 0, 1, 0],  # F
                        [1, 9, 0, 1],  # M
                        [0, 1, 0, 0],  # F
                    ],
                    dtype=np.int8,
                ),
            )

            # Individual IDs
            f.create_dataset("eids", data=np.array([101001, 101002, 101003, 101004]))

            # Biological sex labels
            f.create_dataset("affected_sex", data=np.array([b"M", b"F", b"M", b"F"]))

            # Population labels
            f.create_dataset(
                "populations", data=np.array([b"EUR", b"EUR", b"EUR", b"SAS"])
            )

            # Phecode labels
            f.create_dataset(
                "phecodes",
                # Diabetes, Hypertension, Female-specific, Male-specific
                data=np.array([b"250.2", b"401.1", b"635.2", b"601.1"]),
            )

            # Affected sex labels
            f.create_dataset("phecode_sex", data=np.array([b"B", b"B", b"F", b"M"]))

        yield tmp.name


@pytest.fixture
def phenotype_service(mock_phenotype_file):
    """Create a PhenotypeService instance with the mock data."""
    return PhenotypeService(mock_phenotype_file)


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

    expected_eids = np.array([101001, 101002, 101003, 101004])
    expected_status = np.array([9, 0, 1, 0])

    np.testing.assert_array_equal(eids, expected_eids)
    np.testing.assert_array_equal(case_status, expected_status)


def test_get_cases_nonexistent_phecode(phenotype_service):
    """Test handling of a nonexistent phecode."""
    with pytest.raises(KeyError):
        phenotype_service.get_cases_for_phecode("999.999")


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


def test_get_cases_with_population(phenotype_service):
    """Test case retrieval with population filter."""
    # This is a placeholder test - implement once population filtering is added
    eids, case_status = phenotype_service.get_cases_for_phecode(
        "250.2", population="EUR"
    )
    # For now, population parameter is ignored, so results should be the same
    expected_eids = np.array([101001, 101002, 101003])
    expected_status = np.array([9, 0, 1])

    np.testing.assert_array_equal(eids, expected_eids)
    np.testing.assert_array_equal(case_status, expected_status)


def test_get_cases_with_sex(phenotype_service):
    """Test case retrieval with sex filter."""
    # This is a placeholder test - implement once sex filtering is added
    eids, case_status = phenotype_service.get_cases_for_phecode(
        "250.2", biological_sex="F"
    )
    # For now, sex parameter is ignored, so results should be the same
    expected_eids = np.array([101002, 101004])
    expected_status = np.array([0, 0])

    np.testing.assert_array_equal(eids, expected_eids)
    np.testing.assert_array_equal(case_status, expected_status)
