from unittest.mock import patch, MagicMock
import pytest
from werkzeug.security import generate_password_hash
import numpy as np
import h5py
import tempfile
from pathlib import Path
import os
import pandas as pd

from app import create_app
from data_services.registry import ServiceRegistry
from db import get_db_connection


""" Our lovely fixtures:

test_app (session scope)
└── unit_test_app
    ├── mock_genotype_hdf5_file (session scope)
    │   └── mock_genotype_hdf5_file_with_npy
    ├── mock_phenotype_file
    └── unit_test_app_client

integration_app (session scope)
├── integration_app_client
└── auth_integration_app_client
    └── test_admin


TODOs for later cleanup:
1. The genotype fixture's matrix shape looks wrong (5x4 instead of 4x5 - variants should be columns)
2. Some fixture names are inconsistent (`mock_genotype_hdf5_file` vs `mock_phenotype_file`)
3. The test hierarchy comment could live in a proper docstring at the top of `conftest.py`
4. `test_tests.py` has some commented-out assertions that should be cleaned up
5. Integration tests could use better documentation about what real data they expect
6. Could use more explicit tests for error cases in the streaming endpoints

But for now, let's add the integration test and call it a day! Would you like me to add the test?


"""


@pytest.fixture(scope="session")
def test_app():
    """Base Flask app for all tests with common configuration."""
    # Override the environment to use TestingConfig
    os.environ["FLASK_ENV"] = "testing"
    
    # Create app with testing config
    _app = create_app("testing")
    
    # Additional test-specific configuration
    _app.config.update({
        "SERVER_NAME": "localhost.localdomain",
        "PREFERRED_URL_SCHEME": "http",
        "APPLICATION_ROOT": "/",
    })
    
    return _app


@pytest.fixture(scope="session")
def mock_genotype_hdf5_file():
    """Create a test HDF5 file with merged test data."""
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as f:
        with h5py.File(f.name, "w") as hdf:
            # Add test individuals (using all 4 to cover both cases)
            eids = np.array([1001, 1002, 1003, 1004])
            hdf.create_dataset("fam", data=eids)

            biological_sex = np.array(["M", "F", "M", "F"], dtype=np.string_)
            hdf.create_dataset("sex", data=biological_sex)

            ancestry = np.array(["EUR", "EUR", "EUR", "SAS"], dtype=np.string_)
            hdf.create_dataset("ancestry", data=ancestry)

            # Add test variants (combining both sets)
            bim = np.array(
                [
                    b"1:100:A:T",
                    b"2:200:C:G",
                    b"1:186977737:A:G",
                    b"1:186977780:A:G",
                    b"1:46813503:C:T",
                ]
            )
            hdf.create_dataset("bim", data=bim)
            hdf.create_dataset("genotype_variant_id", data=bim)

            plink_variant_id = np.array(
                [
                    b"1:100_A/T",
                    b"2:200_C/G",
                    b"1:186977737_A/G",
                    b"1:186977780_A/G",
                    b"1:46813503_C/T",
                ]
            )
            hdf.create_dataset("plink_variant_id", data=plink_variant_id)

            # Use underscore format for nomaly_variant_ids
            nomaly_variant_id = np.array(
                [
                    b"1_100_A/T",
                    b"2_200_C/G",
                    b"1_186977737_A/G",
                    b"1_186977780_G/A",  # Note that it's flipped
                    b"Missing",
                ]
            )
            hdf.create_dataset("nomaly_variant_id", data=nomaly_variant_id)

            # Add test genotypes (4 individuals x 5 variants)
            # Combining both matrices and ensuring consistent patterns
            genotypes = np.array(
                [
                    [0, 1, 1, 2],  # Variant 1
                    [1, 0, 0, 1],  # Variant 2
                    [2, 2, 2, 0],  # Variant 3
                    [-1, -1, 2, 1],  # Variant 4
                    [0, 1, 0, 2],  # Variant 5
                ]
            )
            hdf.create_dataset("genotype_matrix", data=genotypes)

    yield f.name
    Path(f.name).unlink()  # Clean up after tests


@pytest.fixture
def mock_genotype_hdf5_file_with_npy(mock_genotype_hdf5_file):
    """The current implementation of GenotypesHDF5 requires a .npy file to be
    present 'next to' the HDF5 file. This fixture creates that .npy file."""
    with h5py.File(mock_genotype_hdf5_file, "r") as f:
        matrix = f["genotype_matrix"]
        assert isinstance(matrix, h5py.Dataset)
        np_matrix = matrix[:]
        np.save(f"{mock_genotype_hdf5_file}.npy", np_matrix)
        yield mock_genotype_hdf5_file
        os.unlink(f"{mock_genotype_hdf5_file}.npy")


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
def unit_test_app(
    test_app, phenotype_service, genotype_service, stats_service, nomaly_scores_service
):
    """App configured for unit tests with mocked services."""

    # Override config for testing
    test_app.config.update(
        {
            "LOGIN_DISABLED": True,
        }
    )

    with test_app.app_context():
        from data_services.registry import ServiceRegistry

        # Create test services registry with our mock services
        services = ServiceRegistry()

        # Use our mock services instead of creating new ones
        services.genotype = genotype_service
        services.phenotype = phenotype_service
        services.stats = stats_service
        services.nomaly_score = nomaly_scores_service

        # Register the services with the app
        test_app.extensions["nomaly_services"] = services

        yield test_app


@pytest.fixture
def phenotype_service():
    """Mock phenotype service for unit tests."""
    from unittest.mock import MagicMock
    import numpy as np

    mock_service = MagicMock()
    mock_hdf = MagicMock()

    # Define our standard test cohort
    # 10 individuals: 3 cases, 5 controls, 2 excluded
    test_eids = np.array([1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010])
    test_phenotypes = {
        "250.2": np.array(
            [1, 1, 1, 0, 0, 0, 0, 0, 9, 9]
        ),  # 1=case, 0=control, 9=excluded
        "290.11": np.array(
            [0, 1, 1, 0, 0, 1, 0, 0, 9, 9]
        ),  # Different pattern for another disease
        "401.1": np.array(
            [1, 0, 1, 1, 0, 0, 0, 0, 9, 9]
        ),  # Different pattern for another disease
    }

    def mock_get_cases_for_phecode(phecode, population=None, biological_sex=None):
        """Return test phenotype data for the given phecode."""
        if phecode not in test_phenotypes:
            # Return default pattern for unknown phecodes
            return test_eids, test_phenotypes["250.2"]
        return test_eids, test_phenotypes[phecode]

    mock_hdf.get_cases_for_phecode = mock_get_cases_for_phecode
    mock_service._hdf = mock_hdf

    return mock_service


@pytest.fixture
def genotype_service():
    """Mock genotype service for unit tests."""
    from unittest.mock import MagicMock
    import numpy as np

    mock_service = MagicMock()
    mock_hdf = MagicMock()

    # Define our standard test cohort - must match phenotype_service
    test_eids = np.array([1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010])

    # Define test genotypes for specific variants
    # 0=homozygous ref, 1=heterozygous, 2=homozygous alt, -1=missing
    test_genotypes = {
        "1_186977737_A/G": np.array([2, 1, 0, 1, 0, 0, 1, 0, 0, 0]),  # Common variant
        "1_186977780_G/A": np.array([1, 2, 0, 0, 1, 0, 0, 0, 0, 0]),  # Less common
        "1_46813503_C/T": np.array([0, 1, 2, 0, 0, 1, 0, 0, 0, 0]),  # Rare variant
        # Add more variants as needed for different test scenarios
    }

    def mock_query_variantID_genotypes(variant_id):
        """Return genotypes for specific variants."""
        if variant_id in test_genotypes:
            return test_eids, test_genotypes[variant_id]
        # Return None for unknown variants (matches real behavior)
        return None

    def mock_get_genotypes(eids=None, vids=None, nomaly_ids=True):
        """Return a genotype matrix for the given eids and variants."""
        if eids is None:
            eids = test_eids
        if vids is None:
            return np.zeros((len(eids), 0))

        # Create matrix with genotypes for requested variants
        genotype_matrix = np.zeros((len(eids), len(vids)))
        for i, vid in enumerate(vids):
            if vid in test_genotypes:
                # Find indices of requested eids in our test data
                eid_indices = np.searchsorted(test_eids, eids)
                genotype_matrix[:, i] = test_genotypes[vid][eid_indices]

        return genotype_matrix

    # Set up the mock
    mock_hdf.query_variantID_genotypes = mock_query_variantID_genotypes
    mock_hdf.get_genotypes = mock_get_genotypes
    mock_hdf.individual = test_eids

    mock_service._hdf = mock_hdf

    return mock_service


@pytest.fixture
def stats_service():
    """Mock stats service for unit tests."""
    from unittest.mock import MagicMock
    import numpy as np
    import pandas as pd

    mock_service = MagicMock()
    mock_hdf = MagicMock()

    # Define test statistics that match our phenotype data
    # These values should align with the phenotype_service fixture
    test_stats = {
        ("GO:0030800", "250.2"): {
            "num_rp": 3.0,  # 3 cases in phenotype data
            "num_rn": 5.0,  # 5 controls in phenotype data
            "metric1_pvalue": 0.0026,
            "metric1_tp": 2,
            "metric1_tn": 685,
            "metric1_fp": 257045,
            "metric1_fn": 532,
            "metric1_threshold": 0.022,
            # MCC stats
            "roc_stats_mcc_pvalue": 0.0015,
            "roc_stats_mcc_or": 2.5,
            "roc_stats_mcc_threshold": 0.025,
            "roc_stats_mcc_tp": 2,
            "roc_stats_mcc_tn": 4,
            "roc_stats_mcc_fp": 1,
            "roc_stats_mcc_fn": 1,
            # YJS stats
            "roc_stats_yjs_pvalue": 0.045,
            "roc_stats_yjs_or": 3.0,
            "roc_stats_yjs_threshold": 0.02,
            "roc_stats_yjs_tp": 2,
            "roc_stats_yjs_tn": 3,
            "roc_stats_yjs_fp": 2,
            "roc_stats_yjs_fn": 1,
            # LRP stats
            "roc_stats_lrp_pvalue": 0.16,
            "roc_stats_lrp_or": 2.0,
            "roc_stats_lrp_threshold": 0.03,
            "roc_stats_lrp_tp": 1,
            "roc_stats_lrp_tn": 4,
            "roc_stats_lrp_fp": 1,
            "roc_stats_lrp_fn": 2,
            # Protective stats
            "roc_stats_lrn_protective_pvalue": 1.0,
            "roc_stats_lrn_protective_or": 0.5,
            "roc_stats_lrn_protective_threshold": 0.0303,
            "roc_stats_lrn_protective_tp": 0,
            "roc_stats_lrn_protective_tn": 247353,
            "roc_stats_lrn_protective_fp": 84,
            "roc_stats_lrn_protective_fn": 497,
        }
    }

    def mock_get_stats_by_term_phecode(term, phecode, statstype=None):
        """Return statistics for a specific term and phecode combination."""
        key = (term, phecode)
        if key not in test_stats:
            raise ValueError(f"No test stats for term={term}, phecode={phecode}")

        if statstype is None:
            return test_stats[key]

        if isinstance(statstype, str):
            return {statstype: test_stats[key][statstype]}

        return {
            stat: test_stats[key][stat] for stat in statstype if stat in test_stats[key]
        }

    # Set up the mock
    mock_hdf.get_stats_by_term_phecode = mock_get_stats_by_term_phecode

    # Add get_stats_by_phecode method for completeness
    def mock_get_stats_by_phecode(phecode, statstype=None):
        """Return statistics for a specific phecode."""
        # Create a simple DataFrame with stats for all terms for this phecode
        result = {}
        for key in test_stats:
            if key[1] == phecode:
                result[key[0]] = test_stats[key]

        if not result:
            raise ValueError(f"No test stats for phecode={phecode}")

        # Convert to DataFrame
        df = pd.DataFrame(result).T

        if statstype is None:
            return df
        elif isinstance(statstype, str):
            return df[statstype]
        else:
            return df[statstype]

    mock_hdf.get_stats_by_phecode = mock_get_stats_by_phecode

    mock_service._hdf = mock_hdf

    return mock_service


@pytest.fixture
def nomaly_scores_service():
    """Mock nomaly scores service for unit tests."""
    from unittest.mock import MagicMock
    import numpy as np

    # Create a mock service with the necessary methods
    mock_service = MagicMock()

    # Create a mock _hdf attribute
    mock_hdf = MagicMock()

    # Define the test case and control scores as in the original test
    case_scores = np.array([0.03, 0.025, 0.01])  # 2 out of 3 cases above threshold
    control_scores = np.array(
        [0.005, 0.015, 0.03, 0.01, 0.005]
    )  # 1 out of 5 controls above threshold

    # Map of eids to their scores - include all eids from our test cohort
    eid_to_score = {
        1001: 0.03,  # case
        1002: 0.025,  # case
        1003: 0.01,  # case
        1004: 0.005,  # control
        1005: 0.015,  # control
        1006: 0.03,  # control
        1007: 0.01,  # control
        1008: 0.005,  # control
        1009: 0.0,  # excluded - no score
        1010: 0.0,  # excluded - no score
    }

    def mock_get_scores_by_eids_unsorted(eids, terms=None):
        """Return scores for the given eids."""
        return np.array([eid_to_score.get(eid, 0.0) for eid in eids])

    # Assign method to mock
    mock_hdf.get_scores_by_eids_unsorted = mock_get_scores_by_eids_unsorted

    # Assign the mock HDF5 object to the service
    mock_service._hdf = mock_hdf

    return mock_service


@pytest.fixture
def unit_test_app_client(unit_test_app):
    """Client for unit tests."""
    return unit_test_app.test_client()


@pytest.fixture(scope="session")
def integration_app():
    """App configured for integration tests with real services."""
    # Use DEVELOPMENT instead of TESTING to allow ServiceRegistry to initialize
    _app = create_app("development")
    _app.config.update(
        {
            "SERVER_NAME": "localhost.localdomain",
            "PREFERRED_URL_SCHEME": "http",
            "APPLICATION_ROOT": "/",
        }
    )

    # # Ensure the app context is active when initializing services
    # with _app.app_context():
    #     # We need to manually initialize services because ServiceRegistry
    #     # doesn't initialize them when TESTING=True
    #     from data_services.registry import ServiceRegistry

    #     # Create service registry and force initialization
    #     services = ServiceRegistry()

    #     # Manually set up mock services for testing
    #     services.genotype = MagicMock()
    #     services.genotype._hdf = MagicMock()

    #     services.phenotype = MagicMock()
    #     services.phenotype._hdf = MagicMock()

    #     services.stats = MagicMock()
    #     services.stats._hdf = MagicMock()

    #     services.nomaly_score = MagicMock()
    #     services.nomaly_score._hdf = MagicMock()

    #     # Register with the app
    #     _app.extensions["nomaly_services"] = services

    return _app


@pytest.fixture
def integration_app_client(integration_app):
    """Client for integration tests."""
    return integration_app.test_client()


@pytest.fixture
def test_admin():
    """Create a test admin user with full permissions."""
    conn = get_db_connection()
    cursor = conn.cursor()

    # Hash the password
    hashed_password = generate_password_hash("test_password")

    # Create test admin user
    cursor.execute(
        """
        INSERT INTO users2 (username, password, email, is_active)
        VALUES ('test_admin', %s, 'test_admin@example.com', TRUE)
    """,
        (hashed_password,),
    )
    user_id = cursor.lastrowid

    # Grant admin permissions
    cursor.execute(
        """
        INSERT INTO user_permissions (user_id, allowed_paths)
        VALUES (%s, '*')
    """,
        (user_id,),
    )

    conn.commit()
    cursor.close()
    conn.close()

    yield {"id": user_id, "username": "test_admin", "password": "test_password"}

    # Cleanup
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute("DELETE FROM user_permissions WHERE user_id = %s", (user_id,))
    cursor.execute("DELETE FROM users2 WHERE id = %s", (user_id,))
    conn.commit()
    cursor.close()
    conn.close()


def cleanup_test_admin_after_test_timeout():
    """Cleanup the test admin user."""
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute(
        """
        DELETE
            user_permissions, users2
        FROM
            user_permissions 
        INNER JOIN
            users2
        ON
          user_permissions.user_id = users2.id
        WHERE
            username = "test_admin"
        """
    )
    conn.commit()
    cursor.close()
    conn.close()


@pytest.fixture
def auth_integration_app_client(integration_app_client, test_admin):
    """Authenticated client for integration tests."""
    response = integration_app_client.post(
        "/login",
        data={"username": test_admin["username"], "password": test_admin["password"]},
        follow_redirects=True,
    )
    assert response.status_code == 200
    assert b"welcome" in response.data.lower()
    return integration_app_client


@pytest.fixture
def mock_cache_dir_config(tmp_path):
    """Create a temporary directory for cache testing"""
    with patch("config.Config") as mock_config:
        mock_config.GWAS_PHENO_DIR = str(tmp_path)
        mock_config.PHEWAS_PHENO_DIR = str(tmp_path)
        mock_config.PHECODE_TERM_DIR = str(tmp_path)
        mock_config.VARIANT_SCORES_DIR = str(tmp_path)

        yield mock_config


@pytest.fixture
def unit_test_app_client_with_cache(unit_test_app_client, mock_cache_dir_config):
    """Client for unit tests with mock cache directories configured."""
    unit_test_app_client.application.config.update(
        {
            "GWAS_PHENO_DIR": mock_cache_dir_config.GWAS_PHENO_DIR,
            "PHEWAS_PHENO_DIR": mock_cache_dir_config.PHEWAS_PHENO_DIR,
            "PHECODE_TERM_DIR": mock_cache_dir_config.PHECODE_TERM_DIR,
            "VARIANT_SCORES_DIR": mock_cache_dir_config.VARIANT_SCORES_DIR,
        }
    )
    return unit_test_app_client


def main():
    # Disaster recovery...
    cleanup_test_admin_after_test_timeout()
    exit(0)


if __name__ == "__main__":
    main()
