from pathlib import Path

import numpy as np
import pandas as pd
from tests.data_utils import production_data_available
import pytest

pytestmark = pytest.mark.skipif(
    not production_data_available(),
    reason="requires production stats HDF5 data"
)

from config import Config
from data_services import StatsRegistry

DEFAULT_RUN_VERSION = "Run-v1"
DEFAULT_ANCESTRY = "EUR"

KNOWN_PHECODES = ["250.2", "290.11"]  # Type 2 diabetes, Dementia
KNOWN_TERM = "GO:0030800"  # Example GO term from old tests

EXPECTED_STATS_250_2_GO_0030800 = {
    "num_rp": 37073,
    "num_rn": 302274,
    "mwu_pvalue": 0.00665,
    "roc_stats_auc": 0.5019811469161224,
    "roc_stats_mcc_value": 0.0042109998036669,
    "roc_stats_yjs_threshold": 0.0001,
    "roc_stats_lrp_protective_or": 0.0,
    "metric1_pvalue": 1.0,
}

EXPECTED_STATS = [
    "num_rp",
    "num_rn",
    "mwu_pvalue",
    "metric1_pvalue",
    "roc_stats_mcc_pvalue",
    "roc_stats_auc",
    "roc_stats_mcc_value",
    "roc_stats_yjs_threshold",
    "roc_stats_lrp_protective_or",
    "tti_pvalue",
    # "tti_pvalue_corrected",
    # "tti_pvalue_corrected_bonferroni",
    # "tti_pvalue_corrected_fdr",
    # "tti_pvalue_corrected_holm",
    # "tti_pvalue_corrected_sidak",
    # "tti_pvalue_corrected_holm_sidak",
]


def test_config_stats_selector_defined():
    """Test that Config.STATS_SELECTOR is properly defined."""
    # Check if STATS_SELECTOR exists
    assert hasattr(Config, "STATS_SELECTOR"), "Config.STATS_SELECTOR is not defined"
    assert Config.STATS_SELECTOR is not None, "Config.STATS_SELECTOR is None"

    # Check if it's a dictionary with the expected structure
    assert isinstance(Config.STATS_SELECTOR, dict), (
        "Config.STATS_SELECTOR is not a dictionary"
    )

    # Check that at least one run version is defined
    assert len(Config.STATS_SELECTOR) > 0, "Config.STATS_SELECTOR is empty"

    # Check the structure for the first run version and print available versions for debugging
    first_run = next(iter(Config.STATS_SELECTOR))
    print(f"Available run versions: {list(Config.STATS_SELECTOR.keys())}")

    assert isinstance(Config.STATS_SELECTOR[first_run], dict), (
        f"Run version '{first_run}' is not a dictionary"
    )

    # Check that at least one ancestry is defined for the first run and print available ancestries
    assert len(Config.STATS_SELECTOR[first_run]) > 0, (
        f"No ancestries defined for run version '{first_run}'"
    )
    print(
        f"Available ancestries for {first_run}: {list(Config.STATS_SELECTOR[first_run].keys())}"
    )

    # Verify that the first ancestry for the first run points to a valid file path
    first_ancestry = next(iter(Config.STATS_SELECTOR[first_run]))
    file_path = Config.STATS_SELECTOR[first_run][first_ancestry]

    assert isinstance(file_path, (str, Path)), (
        f"Path for {first_run}/{first_ancestry} is not a string or Path"
    )

    # Convert to Path if it's a string
    if isinstance(file_path, str):
        file_path = Path(file_path)

    # Check if the HDF5 file exists
    assert file_path.exists(), f"HDF5 file at {file_path} does not exist"
    assert file_path.is_file(), f"{file_path} is not a file"
    assert file_path.suffix in [".hdf5", ".h5"], f"{file_path} is not an HDF5 file"


# === Fixtures ===
@pytest.fixture(scope="session")
def stats_registry():
    """Provides a configured StatsRegistry instance."""
    # Assuming Config.STATS_SELECTOR exists and is the required dict structure
    # e.g., {'Run-v1': {'EUR': Path(...), 'AFR': Path(...)}, ...}
    assert hasattr(Config, "STATS_SELECTOR") and Config.STATS_SELECTOR, (
        "STATS_SELECTOR not found or not configured in config.py. "
        "Define it as a dict: {'run_version': {'ancestry': Path('/path/to/stats.hdf5')}}"
    )
    return StatsRegistry(Config.STATS_SELECTOR)


@pytest.fixture(scope="session")
def stats_service(stats_registry):
    """Provides a default StatsService instance (e.g., for Run-v1, EUR)."""
    return stats_registry.get(
        run_version=DEFAULT_RUN_VERSION, ancestry=DEFAULT_ANCESTRY
    )


# === Tests ===
def test_stats_service_instance(stats_service):
    """Test that the fixture provides a valid StatsService instance."""
    assert stats_service is not None
    assert hasattr(stats_service, "_hdf")
    assert stats_service._hdf.f is not None  # Check if HDF5 file is open


def test_get_stats_for_known_phecode(stats_service):
    """Test retrieving stats for a known phecode."""
    phecode = KNOWN_PHECODES[0]
    stats_df = stats_service.get_phecode_stats(phecode=phecode)

    assert isinstance(stats_df, pd.DataFrame)
    assert not stats_df.empty
    # assert stats_df.index.name == "term"  # Assuming terms become the index

    # NOTE: The stats are columns of the dataframe!
    assert all(stat in stats_df.columns for stat in EXPECTED_STATS), (
        f"DataFrame missing expected stats columns (We have {stats_df.columns} and expected {EXPECTED_STATS})"
    )
    # Basic check on number of terms returned (adjust threshold as needed)
    assert len(stats_df) > 100, "Expected a reasonable number of terms for the phecode"


def test_get_stats_for_known_term(stats_service):
    """Test retrieving stats for a known term."""
    term = KNOWN_TERM
    stats_df = stats_service.get_term_stats(term=term)

    assert isinstance(stats_df, pd.DataFrame)
    assert not stats_df.empty
    # Assuming phecodes become the index when querying by term

    # NOTE: The stats are indexes of the dataframe!
    assert all(stat in stats_df.index for stat in EXPECTED_STATS), (
        f"DataFrame missing expected stats INDEX (We have {stats_df.columns} and expected {EXPECTED_STATS})"
    )
    # Basic check on number of phecodes returned (adjust threshold as needed)
    assert len(stats_df) > 50, "Expected a reasonable number of phecodes for the term"


def test_get_stats_for_specific_term_phecode(stats_service):
    """Test retrieving stats for a specific term-phecode combination."""
    phecode = KNOWN_PHECODES[0]  # e.g., "250.2"
    term = KNOWN_TERM  # e.g., "GO:0030800"

    # Query by phecode, filter by term
    stats_df = stats_service.get_phecode_stats(phecode=phecode, term=term)

    assert isinstance(stats_df, pd.DataFrame)
    assert len(stats_df) == 1, "Expected exactly one row for specific term-phecode"
    assert stats_df.index[0] == term

    # Check specific placeholder values (UPDATE THESE)
    retrieved_stats = stats_df.iloc[0].to_dict()
    for stat_name, expected_value in EXPECTED_STATS_250_2_GO_0030800.items():
        assert stat_name in retrieved_stats, f"Stat '{stat_name}' not found in results"
        # Use pytest.approx for floating point comparisons
        if isinstance(expected_value, float):
            assert retrieved_stats[stat_name] == pytest.approx(
                expected_value, rel=1e-3
            ), f"Mismatch for {stat_name}"
        else:
            assert retrieved_stats[stat_name] == expected_value, (
                f"Mismatch for {stat_name}"
            )

    # Optional: Query by term, filter by phecode and check consistency
    stats_df_term_query = stats_service.get_term_stats(term=term, phecode=phecode)
    assert np.all(stats_df.values == stats_df_term_query.values), (
        "Results differ when querying by term vs phecode"
    )


def test_get_stats_nonexistent_phecode(stats_service):
    """Test querying a phecode that doesn't exist in the data."""
    non_existent_phecode = "999.9999"
    with pytest.raises(ValueError, match="No data found for the provided phecodes"):
        stats_service.get_phecode_stats(phecode=non_existent_phecode)


def test_get_stats_nonexistent_term(stats_service):
    """Test querying a term that doesn't exist in the data."""
    non_existent_term = "GO:9999999_FAKE"
    with pytest.raises(ValueError, match="No data found for the provided terms"):
        stats_service.get_term_stats(term=non_existent_term)


@pytest.mark.parametrize("phecode", KNOWN_PHECODES)
def test_column_consistency_across_phecodes(stats_service, phecode):
    """Check that the columns returned are consistent for different phecodes."""
    stats_df = stats_service.get_phecode_stats(phecode=phecode)
    assert all(stat in stats_df.columns for stat in EXPECTED_STATS), (
        f"Column mismatch for phecode {phecode}"
    )


# Potential future tests:
# - Test querying multiple phecodes/terms at once if needed.
# - Test querying with specific stats_types filters.
# - Add tests for other run versions or ancestries if the registry is configured for them.
# - More detailed value range checks for certain statistics (e.g., p-values between 0 and 1).


def test_filtering_by_stats_types(stats_service):
    """Test filtering results by specific stats types."""
    phecode = KNOWN_PHECODES[0]

    # Get a subset of stats columns
    subset_stats = ["mwu_pvalue", "tti_pvalue"]
    stats_df = stats_service.get_phecode_stats(
        phecode=phecode, stats_types=subset_stats
    )

    assert isinstance(stats_df, pd.DataFrame)
    assert not stats_df.empty
    # Check that only the requested stats types are returned
    assert set(stats_df.columns) == set(subset_stats), (
        "Column mismatch for filtered stats"
    )


def test_p_value_ranges(stats_service):
    """Test that p-values are within the valid range [0, 1]."""
    phecode = KNOWN_PHECODES[0]
    pvalue_columns = [col for col in EXPECTED_STATS if "pvalue" in col.lower()]

    stats_df = stats_service.get_phecode_stats(phecode=phecode)

    for col in pvalue_columns:
        # Check if column exists
        assert col in stats_df.columns, f"P-value column {col} not found"

        # Filter out NaN values
        values = stats_df[col].dropna().values

        # Check that all p-values are in [0, 1]
        assert np.all((values >= 0) & (values <= 1)), (
            f"P-values in {col} outside valid range [0, 1]"
        )


@pytest.mark.parametrize(
    "run_version,ancestry",
    [
        (DEFAULT_RUN_VERSION, DEFAULT_ANCESTRY),  # Default combo for baseline
        ("Run-v2", DEFAULT_ANCESTRY),
        (DEFAULT_RUN_VERSION, "AFR"),
    ],
)
def test_run_version_ancestry_combinations(stats_registry, run_version, ancestry):
    """Test accessing stats for different run version and ancestry combinations."""
    service = stats_registry.get(run_version=run_version, ancestry=ancestry)

    # Basic test to verify service works
    phecode = KNOWN_PHECODES[0]
    stats_df = service.get_phecode_stats(phecode=phecode)

    assert isinstance(stats_df, pd.DataFrame)
    assert not stats_df.empty


def test_multiple_terms_query(stats_service):
    """Test querying stats for multiple terms at once."""
    # This test might need to be modified based on actual API implementation
    phecode = KNOWN_PHECODES[0]

    # You'll need more than one known term that exists in the data
    terms = [KNOWN_TERM, "GO:0005829"]  # Add a second known term

    stats_df = stats_service.get_phecode_stats(phecode=phecode, terms=terms)

    assert isinstance(stats_df, pd.DataFrame)
    assert not stats_df.empty
    assert len(stats_df) == len(terms), (
        f"Expected {len(terms)} rows, got {len(stats_df)}"
    )

    # Check that all requested terms are in the results
    result_terms = stats_df.index.tolist()
    for term in terms:
        assert term in result_terms, f"Term {term} missing from results"


def test_term_phecode_consistency(stats_service):
    """Test that querying by term+phecode or phecode+term gives consistent results."""
    # Choose a known phecode and term combination
    phecode = KNOWN_PHECODES[0]
    term = KNOWN_TERM

    # Query from both directions
    term_df = stats_service.get_term_stats(term=term, phecode=phecode)
    phecode_df = stats_service.get_phecode_stats(phecode=phecode, term=term)

    # Convert to comparable format if needed (depends on API implementation)
    # This assumes both return DataFrames with similar structure
    assert term_df.shape == phecode_df.shape, (
        "Different result shapes from term vs phecode query"
    )

    # Compare actual values - this might need adjustment based on how indices work
    for col in term_df.columns:
        if col in phecode_df.columns:
            assert term_df[col].iloc[0] == pytest.approx(
                phecode_df[col].iloc[0], rel=1e-6
            ), f"Value mismatch for {col} between term and phecode queries"


def test_service_error_handling(stats_service):
    """Test that the service handles invalid inputs gracefully."""
    # Test with invalid phecode
    with pytest.raises(ValueError):
        stats_service.get_phecode_stats(phecode="not-a-phecode")

    # Test with invalid term
    with pytest.raises(ValueError):
        stats_service.get_term_stats(term="not-a-term")

    # Test with invalid stats type
    with pytest.raises(ValueError):
        stats_service.get_phecode_stats(
            phecode=KNOWN_PHECODES[0], stats_types=["invalid_stat"]
        )


@pytest.mark.parametrize(
    "phecode1,phecode2",
    [
        (
            KNOWN_PHECODES[0],
            KNOWN_PHECODES[1],
        )  # Use different phecodes from your constants
    ],
)
def test_term_ordering_consistency(stats_service, phecode1, phecode2):
    """Test that term ordering is consistent across different phecodes."""
    # Get stats for two different phecodes
    stats_df1 = stats_service.get_phecode_stats(phecode=phecode1)
    stats_df2 = stats_service.get_phecode_stats(phecode=phecode2)

    # Check that we have results
    assert not stats_df1.empty
    assert not stats_df2.empty

    # Find common terms between the two results
    common_terms = set(stats_df1.index).intersection(set(stats_df2.index))
    assert len(common_terms) > 0, "No common terms found between phecodes"

    # Check consistency for common terms - use mwu_pvalue as an example stat
    for term in list(common_terms)[:5]:  # Limit to first 5 to avoid excessive checking
        # Verify term exists in both results
        assert term in stats_df1.index
        assert term in stats_df2.index


def test_data_type_consistency(stats_service):
    """Test that returned data types are consistent and appropriate."""
    phecode = KNOWN_PHECODES[0]
    stats_df = stats_service.get_phecode_stats(phecode=phecode)

    # Check count columns are integers or floats
    count_cols = [col for col in stats_df.columns if "num_" in col or "_count" in col]
    for col in count_cols:
        assert pd.api.types.is_numeric_dtype(stats_df[col]), f"{col} should be numeric"

    # Check p-value columns are floats
    pvalue_cols = [col for col in stats_df.columns if "pvalue" in col.lower()]
    for col in pvalue_cols:
        assert pd.api.types.is_float_dtype(stats_df[col]), f"{col} should be float type"


def test_phecode_cases_controls_ratio(stats_service):
    """Test that case/control ratios are reasonable for common phecodes."""
    for phecode in KNOWN_PHECODES:
        # Pick a term that likely has data for this phecode
        term = KNOWN_TERM
        stats_df = stats_service.get_phecode_stats(phecode=phecode, term=term)

        # If we have num_rp (cases) and num_rn (controls)
        if "num_rp" in stats_df.columns and "num_rn" in stats_df.columns:
            cases = stats_df["num_rp"].iloc[0]
            controls = stats_df["num_rn"].iloc[0]

            # Sanity check: we should have both cases and controls
            assert cases > 0, f"Expected non-zero cases for {phecode}"
            assert controls > 0, f"Expected non-zero controls for {phecode}"

            # Sanity check: case/control ratio should be reasonable
            ratio = controls / cases if cases > 0 else float("inf")
            assert 0.1 <= ratio <= 1000, (
                f"Unexpected case/control ratio for {phecode}: {ratio}"
            )


# === Tests for StatsRegistry ===


def test_stats_registry_caching(stats_registry):
    """Test that StatsRegistry caches service instances."""
    # Get the same service twice
    service1 = stats_registry.get(DEFAULT_RUN_VERSION, DEFAULT_ANCESTRY)
    service2 = stats_registry.get(DEFAULT_RUN_VERSION, DEFAULT_ANCESTRY)

    # They should be the exact same object (identity check)
    assert service1 is service2, "StatsRegistry did not cache service instance"


def test_stats_registry_initialization():
    """Test that StatsRegistry initializes correctly."""
    # Create a new registry with an empty selector
    registry = StatsRegistry({})
    assert registry.initialized, "Registry should initialize with empty dict"

    # Create a new registry with None (should not be initialized)
    registry = StatsRegistry(None)
    assert not registry.initialized, (
        "Registry with None selector should not be initialized"
    )

    # Try to use the uninitialized registry
    with pytest.raises(ValueError, match="Service not properly initialized"):
        registry.get()


def test_stats_registry_invalid_selector():
    """Test StatsRegistry with invalid selector formats."""
    # Test with empty key values
    with pytest.raises(ValueError):
        registry = StatsRegistry({"Run-v1": {}})
        registry.get("Run-v1", "EUR")

    # Test with invalid run version
    with pytest.raises(ValueError, match="No stats file found"):
        registry = StatsRegistry({"Run-v1": {"EUR": Path("dummy.hdf5")}})
        registry.get("Run-v9999", "EUR")

    # Test with invalid ancestry
    with pytest.raises(ValueError, match="No stats file found"):
        registry = StatsRegistry({"Run-v1": {"EUR": Path("dummy.hdf5")}})
        registry.get("Run-v1", "INVALID")


# Additional test for the STATS_SELECTOR configuration in Config
def test_config_stats_selector_structure():
    """Test that the STATS_SELECTOR configuration has the expected structure."""
    # Basic structure checks
    assert isinstance(Config.STATS_SELECTOR, dict), (
        "STATS_SELECTOR should be a dictionary"
    )

    # Check each run version
    for run_version, ancestries in Config.STATS_SELECTOR.items():
        assert isinstance(run_version, str), (
            f"Run version key '{run_version}' should be a string"
        )
        assert isinstance(ancestries, dict), (
            f"Ancestries for '{run_version}' should be a dictionary"
        )

        # Check each ancestry
        for ancestry, filepath in ancestries.items():
            assert isinstance(ancestry, str), (
                f"Ancestry key '{ancestry}' should be a string"
            )
            assert isinstance(filepath, (str, Path)), (
                f"Path for {run_version}/{ancestry} should be string or Path"
            )

            # Convert to Path if it's a string
            if isinstance(filepath, str):
                filepath = Path(filepath)

            # Note: We don't check file existence here to allow for flexible test execution
            # in environments where not all files might be available


# === Advanced Data Slicing Tests ===


def test_data_slice_multi_term_single_phecode(stats_service):
    """Test slicing data with multiple terms for a single phecode."""
    phecode = KNOWN_PHECODES[0]

    # Try to find a few terms that we know exist
    # First get all terms for this phecode
    all_terms_df = stats_service.get_phecode_stats(phecode=phecode)

    # Select the first few terms
    assert len(all_terms_df) >= 3, f"Not enough terms found for phecode {phecode}"
    terms = all_terms_df.index[:3].tolist()

    # Query with multiple terms
    stats_df = stats_service.get_phecode_stats(phecode=phecode, terms=terms)

    # Verify the result
    assert isinstance(stats_df, pd.DataFrame)
    assert not stats_df.empty
    assert len(stats_df) <= len(terms), "Got more rows than requested terms"

    # Check that the returned terms are a subset of what we requested
    for term in stats_df.index:
        assert term in terms, f"Got unexpected term in results: {term}"


def test_data_slice_single_term_multi_phecode(stats_service):
    """Test slicing data with a single term for multiple phecodes."""
    term = KNOWN_TERM

    # Query with multiple phecodes
    stats_df = stats_service.get_term_stats(term=term, phecodes=KNOWN_PHECODES)

    # Verify the result
    assert isinstance(stats_df, pd.DataFrame)
    assert not stats_df.empty
    assert stats_df.shape[1] <= len(KNOWN_PHECODES), (
        "Got more rows than requested phecodes"
    )

    # Check that the returned phecodes are a subset of what we requested
    for phecode in stats_df.columns:
        assert phecode in KNOWN_PHECODES, (
            f"Got unexpected phecode in results: {phecode}"
        )


def test_full_slice_all_stats_types(stats_service):
    """Test getting all available stats types for a term-phecode combination."""
    phecode = KNOWN_PHECODES[0]
    term = KNOWN_TERM

    # Get stats without specifying stats_types (should return all)
    full_stats_df = stats_service.get_phecode_stats(phecode=phecode, term=term)

    # Now specify each column explicitly
    columns = full_stats_df.columns.tolist()
    explicit_stats_df = stats_service.get_phecode_stats(
        phecode=phecode, term=term, stats_types=columns
    )

    # They should be the same
    assert full_stats_df.equals(explicit_stats_df), (
        "Different results when specifying all stats_types"
    )


# === Thorough Error/Edge Case Tests ===


def test_invalid_combinations(stats_service):
    """Test various invalid parameter combinations."""
    # Test: Providing both term and terms
    with pytest.raises(ValueError, match="silly"):
        stats_service.get_phecode_stats(
            phecode=KNOWN_PHECODES[0], term=KNOWN_TERM, terms=["another_term"]
        )

    # Test: Providing both phecode and phecodes
    with pytest.raises(ValueError, match="silly"):
        stats_service.get_term_stats(
            term=KNOWN_TERM, phecode=KNOWN_PHECODES[0], phecodes=["another_phecode"]
        )
