import numpy as np
import pandas as pd
import pytest

from config import Config
from data_services import PhenotypeService

# Known test cases - these values will need to be updated with actual production data
KNOWN_PHECODE = "250.2"  # Type 2 diabetes
NUM_INDIVIDUALS = 449423
EXPECTED_CASE_COUNT = 37073  # Placeholder - update with actual count
EXPECTED_CONTROL_COUNT = 408655  # Placeholder - update with actual count

FEMALE_SPECIFIC_PHECODE = "635.2"
MALE_SPECIFIC_PHECODE = "601.1"


@pytest.fixture(scope="session")
def phenotype_service():
    return PhenotypeService(Config.PHENOTYPES_HDF)


def test_production_file_exists(phenotype_service):
    """Verify the production phenotype file exists and is readable."""
    assert phenotype_service is not None
    assert hasattr(phenotype_service, "_hdf")
    assert hasattr(phenotype_service._hdf, "phenotype_data")


def test_known_phecode_query(phenotype_service):
    """Test querying a known phecode from production data."""
    phenotypes_df = phenotype_service.get_cases_for_phecode(KNOWN_PHECODE)

    assert phenotypes_df is not None
    assert isinstance(phenotypes_df, pd.DataFrame)
    assert len(phenotypes_df) == NUM_INDIVIDUALS

    # Verify case/control counts
    case_count = np.sum(phenotypes_df["phenotype"] == 1)
    control_count = np.sum(phenotypes_df["phenotype"] == 0)

    # These assertions will need to be updated with actual values
    assert case_count == EXPECTED_CASE_COUNT
    assert control_count == EXPECTED_CONTROL_COUNT


def test_get_cases_nonexistent_phecode(phenotype_service):
    """Test handling of a nonexistent phecode."""
    with pytest.raises(KeyError):
        phenotype_service.get_cases_for_phecode("999.999")


# TODO: Fetch some individual_sex data
def test_sex_specific_phecodes(phenotype_service):
    """Test sex-specific phecodes have appropriate distributions."""
    phenotypes_df = phenotype_service.get_cases_for_phecode(FEMALE_SPECIFIC_PHECODE)

    # Female-specific phecode
    case_count_f = np.sum(phenotypes_df["phenotype"] == 1)
    assert case_count_f > 0

    # Male-specific phecode
    phenotypes_df = phenotype_service.get_cases_for_phecode(MALE_SPECIFIC_PHECODE)
    case_count_m = np.sum(phenotypes_df["phenotype"] == 1)
    assert case_count_m > 0

@pytest.mark.parametrize("phecode", ["250.2", "401.1"])
def test_individual_consistency(phenotype_service, phecode):
    """Test that individual IDs are consistent across different phecode queries."""
    phenotypes_df = phenotype_service.get_cases_for_phecode(phecode)

    # Get cases for two different phecodes
    eids = phenotypes_df["eid"].values

    # Check that the individual IDs are in the same range
    assert eids.min() > 1_000_000
    assert eids.max() < 10_000_000


@pytest.mark.parametrize("phecode", ["250.2", "401.1", "272.1"])
def test_case_control_ratios(phenotype_service, phecode):
    """Test that case/control ratios are within expected ranges for common diseases."""
    phenotypes_df = phenotype_service.get_cases_for_phecode(phecode)

    case_count = np.sum(phenotypes_df["phenotype"] == 1)
    control_count = np.sum(phenotypes_df["phenotype"] == 0)

    # Common diseases should have reasonable case counts
    assert case_count > 100, f"Expected >100 cases for {phecode}, got {case_count}"

    # Case/control ratio should be reasonable (e.g., between 1:100 and 1:1)
    ratio = control_count / case_count
    assert 1 <= ratio <= 100, f"Unexpected case/control ratio for {phecode}: {ratio}"


def test_excluded_individuals(phenotype_service):
    """Test that excluded individuals are handled correctly."""
    phenotypes_df = phenotype_service.get_cases_for_phecode(KNOWN_PHECODE)

    # Get cases for a phecode
    assert np.all(np.isin(phenotypes_df["phenotype"], [0, 1, 9]))


@pytest.mark.parametrize(
    "ancestry_counts",
    [("EUR", 37073), ("AFR", 1431), ("EAS", 219), ("SAS", 2261)],
)
def test_population_stratification(phenotype_service, ancestry_counts):
    """Test population-specific queries (when implemented)."""

    ancestry, count = ancestry_counts
    phenotypes_df = phenotype_service.get_cases_for_phecode(KNOWN_PHECODE, ancestry)
    assert np.all(np.isin(phenotypes_df["phenotype"], [0, 1, 9]))

    case_count = np.sum(phenotypes_df["phenotype"] == 1)
    assert case_count == count
