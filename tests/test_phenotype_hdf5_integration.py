import numpy as np
import pytest

from data_services.implementations.hdf5.phenotype import HDF5PhenotypeService
from config import Config

config = Config()
phenotype_service = HDF5PhenotypeService(config.PHENOTYPES_H5)

# Known test cases - these values will need to be updated with actual production data
KNOWN_PHECODE = "250.2"  # Type 2 diabetes
EXPECTED_CASE_COUNT = 42973  # Placeholder - update with actual count
EXPECTED_CONTROL_COUNT = 439119  # Placeholder - update with actual count

FEMALE_SPECIFIC_PHECODE = "620.2"  # Female genital disease
MALE_SPECIFIC_PHECODE = "600"  # Male genital disease


def test_production_file_exists():
    """Verify the production phenotype file exists and is readable."""
    assert phenotype_service is not None
    assert hasattr(phenotype_service, "_hdf")
    assert hasattr(phenotype_service._hdf, "phenotype_data")


def test_known_phecode_query():
    """Test querying a known phecode from production data."""
    eids, case_status = phenotype_service.get_cases_for_phecode(KNOWN_PHECODE)

    assert eids is not None
    assert case_status is not None
    assert isinstance(eids, np.ndarray)
    assert isinstance(case_status, np.ndarray)
    assert len(eids) == len(case_status)

    # Verify case/control counts
    case_count = np.sum(case_status == 1)
    control_count = np.sum(case_status == 0)

    # These assertions will need to be updated with actual values
    assert case_count == EXPECTED_CASE_COUNT, (
        f"Expected {EXPECTED_CASE_COUNT} cases, got {case_count}"
    )
    assert control_count == EXPECTED_CONTROL_COUNT, (
        f"Expected {EXPECTED_CONTROL_COUNT} controls, got {control_count}"
    )


def test_sex_specific_phecodes():
    """Test sex-specific phecodes have appropriate distributions."""
    # Female-specific phecode
    eids_f, status_f = phenotype_service.get_cases_for_phecode(FEMALE_SPECIFIC_PHECODE)
    case_count_f = np.sum(status_f == 1)
    assert case_count_f > 0, (
        "Expected non-zero female cases for female-specific phecode"
    )

    # Male-specific phecode
    eids_m, status_m = phenotype_service.get_cases_for_phecode(MALE_SPECIFIC_PHECODE)
    case_count_m = np.sum(status_m == 1)
    assert case_count_m > 0, "Expected non-zero male cases for male-specific phecode"


def test_individual_consistency():
    """Test that individual IDs are consistent across different phecode queries."""
    # Get cases for two different phecodes
    eids1, _ = phenotype_service.get_cases_for_phecode("250.2")  # Type 2 diabetes
    eids2, _ = phenotype_service.get_cases_for_phecode("401.1")  # Hypertension

    # Check that the individual IDs are in the same range
    assert eids1.min() > 1_000_000
    assert eids2.min() > 1_000_000
    assert eids1.max() < 10_000_000
    assert eids2.max() < 10_000_000

    # Check that there's significant overlap in the individuals
    common_ids = np.intersect1d(eids1, eids2)
    assert len(common_ids) > 0, "Expected overlap in individuals between phecodes"


def test_case_control_ratios():
    """Test that case/control ratios are within expected ranges for common diseases."""
    common_phecodes = [
        "250.2",  # Type 2 diabetes
        "401.1",  # Hypertension
        "272.1",  # Hyperlipidemia
    ]

    for phecode in common_phecodes:
        _, status = phenotype_service.get_cases_for_phecode(phecode)
        case_count = np.sum(status == 1)
        control_count = np.sum(status == 0)

        # Common diseases should have reasonable case counts
        assert case_count > 100, f"Expected >100 cases for {phecode}, got {case_count}"

        # Case/control ratio should be reasonable (e.g., between 1:100 and 1:1)
        ratio = control_count / case_count
        assert 1 <= ratio <= 100, (
            f"Unexpected case/control ratio for {phecode}: {ratio}"
        )


def test_excluded_individuals():
    """Test that excluded individuals are handled correctly."""
    # Get cases for a phecode
    eids, status = phenotype_service.get_cases_for_phecode(KNOWN_PHECODE)

    assert eids.shape == status.shape

    # All values should be either 0 (control), 1 (case) or 9 (exclude)
    assert np.all(np.isin(status, [0, 1, 9])), "Invalid status values in results"


def test_population_stratification():
    """Test population-specific queries (when implemented)."""
    # This is a placeholder test - implement once population filtering is added
    populations = ["EUR", "AFR", "EAS", "SAS"]

    for pop in populations:
        eids, status = phenotype_service.get_cases_for_phecode(
            KNOWN_PHECODE, population=pop
        )
        # For now, population parameter is ignored, so just verify basic sanity
        assert len(eids) > 0, f"No data returned for population {pop}"
        assert eids.shape == status.shape

        assert np.all(np.isin(status, [0, 1, 9]))


def test_sex_stratification():
    """Test sex-specific queries (when implemented)."""
    # This is a placeholder test - implement once sex filtering is added
    for sex in ["male", "female"]:
        eids, status = phenotype_service.get_cases_for_phecode(KNOWN_PHECODE, sex=sex)
        # For now, sex parameter is ignored, so just verify basic sanity
        assert len(eids) > 0, f"No data returned for sex {sex}"
        assert eids.shape == status.shape
        assert np.all(np.isin(status, [0, 1, 9]))
