import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
import pytest
from werkzeug.security import generate_password_hash

from app import create_app
from db import get_db_connection

""" Our lovely fixtures:

test_app (session scope)
└── unit_test_app
    ├── mock_genotype_hdf5_file (session scope)
    │   └── mock_genotype_hdf5_file_with_npy
    ├── mock_phenotype_hdf5_file
    └── unit_test_app_client

integration_app (session scope)
├── integration_app_client
└── auth_integration_app_client
    └── test_admin


TODOs for later cleanup:
1. The genotype fixture's matrix shape looks wrong (5x4 instead of 4x5 - variants should be columns)
2. Some fixture names are inconsistent (`mock_genotype_hdf5_file` vs `mock_phenotype_hdf5_file`)
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
            hdf.create_dataset("eid", data=eids)

            individual_sex = np.array(["M", "F", "M", "F"], dtype=np.bytes_)
            hdf.create_dataset("sex", data=individual_sex)

            ancestry = np.array(["EUR", "EUR", "EUR", "SAS"], dtype=np.bytes_)
            hdf.create_dataset("ancestry", data=ancestry)

            ref = np.array(
                [
                    b"A",
                    b"C",
                    b"A",
                    b"A",
                    b"C",
                ]
            )
            alt = np.array(
                [
                    b"T",
                    b"G",
                    b"G",
                    b"A",
                    b"T",
                ]
            )

            hdf.create_dataset("REF", data=ref)
            hdf.create_dataset("ALT", data=alt)

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

            rsIDs = np.array(
                [
                    b"rs123456789",
                    b"rs123456790",
                    b"rs123456791",
                    b"rs123456792",
                    b"Missing",
                ]
            )
            hdf.create_dataset("rsID", data=rsIDs)

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
            hdf.create_dataset("genotypes", data=genotypes)

    yield Path(f.name)
    Path(f.name).unlink()  # Clean up after tests


@pytest.fixture
def mock_genotype_hdf5_file_with_npy(mock_genotype_hdf5_file):
    """The current implementation of GenotypesHDF5 requires a .npy file to be
    present 'next to' the HDF5 file. This fixture creates that .npy file."""
    with h5py.File(mock_genotype_hdf5_file, "r") as f:
        matrix = f["genotypes"]
        assert isinstance(matrix, h5py.Dataset)
        np_matrix = matrix[:]
        np.save(mock_genotype_hdf5_file.with_suffix(".npy"), np_matrix)
        yield mock_genotype_hdf5_file
        os.unlink(mock_genotype_hdf5_file.with_suffix(".npy"))


@pytest.fixture
def mock_phenotype_hdf5_file():
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
    test_app,
    phenotype_service,
    genotype_service,
    nomaly_scores_service,
    nomaly_data_service,
    stats_registry,
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
        services.nomaly_score = nomaly_scores_service
        services.nomaly_data = nomaly_data_service

        services.stats_registry = stats_registry

        # Register the services with the app
        test_app.extensions["nomaly_services"] = services

        yield test_app


@pytest.fixture
def phenotype_service(mock_phenotype_hdf5_file):
    """Mock phenotype service for unit tests."""
    from data_services.phenotype import PhenotypeService

    # Create a real PhenotypeService instance with our mock file
    service = PhenotypeService(mock_phenotype_hdf5_file)

    return service


@pytest.fixture
def genotype_service(mock_genotype_hdf5_file_with_npy):
    """Mock genotype service for unit tests."""
    from data_services.genotype import GenotypeService

    # Create a real GenotypeService instance with our mock file
    service = GenotypeService(mock_genotype_hdf5_file_with_npy)

    return service


@pytest.fixture
def mock_stats_hdf5_file():
    """Create a temporary mock HDF5 file with test stats data."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Define our test data
            # The stats are stored in a 3D array: (term, phecode, stat_type)

            # Define the terms we support
            terms = np.array([b"GO:0030800", b"MP:0005179", b"GO:0006915"])

            # Define the phecodes we support
            phecodes = np.array([b"250.2", b"290.11", b"401.1", b"635.2", b"601.1"])

            # Define the stat types
            stats_types = np.array(
                [
                    b"num_rp",
                    b"num_rn",
                    b"metric1_pvalue",
                    b"metric1_tp",
                    b"metric1_tn",
                    b"metric1_fp",
                    b"metric1_fn",
                    b"metric1_threshold",
                    b"roc_stats_mcc_pvalue",
                    b"roc_stats_mcc_or",
                    b"roc_stats_mcc_threshold",
                    b"roc_stats_mcc_tp",
                    b"roc_stats_mcc_tn",
                    b"roc_stats_mcc_fp",
                    b"roc_stats_mcc_fn",
                    b"roc_stats_yjs_pvalue",
                    b"roc_stats_yjs_or",
                    b"roc_stats_yjs_threshold",
                    b"roc_stats_yjs_tp",
                    b"roc_stats_yjs_tn",
                    b"roc_stats_yjs_fp",
                    b"roc_stats_yjs_fn",
                    b"roc_stats_lrp_pvalue",
                    b"roc_stats_lrp_or",
                    b"roc_stats_lrp_threshold",
                    b"roc_stats_lrp_tp",
                    b"roc_stats_lrp_tn",
                    b"roc_stats_lrp_fp",
                    b"roc_stats_lrp_fn",
                    b"roc_stats_lrn_protective_pvalue",
                    b"roc_stats_lrn_protective_or",
                    b"roc_stats_lrn_protective_threshold",
                    b"roc_stats_lrn_protective_tp",
                    b"roc_stats_lrn_protective_tn",
                    b"roc_stats_lrn_protective_fp",
                    b"roc_stats_lrn_protective_fn",
                ]
            )

            # Create a 3D array to hold the data
            data = np.zeros((len(terms), len(phecodes), len(stats_types)))

            # Fill in test data for GO:0030800 + 250.2 combination
            # Get the indices for our combination
            term_idx = np.where(terms == b"GO:0030800")[0][0]
            phecode_idx = np.where(phecodes == b"250.2")[0][0]

            # Create a dictionary mapping stat type to value for our test combo
            test_stats = {
                "num_rp": 3.0,
                "num_rn": 5.0,
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

            # Populate data for our test combination
            for stat_name, value in test_stats.items():
                stat_idx = np.where(stats_types == stat_name.encode())[0][0]
                data[term_idx, phecode_idx, stat_idx] = value

            # Save datasets to file
            f.create_dataset("data", data=data)
            f.create_dataset("term", data=terms)
            f.create_dataset("phecode", data=phecodes)
            f.create_dataset("stats_type", data=stats_types)

        yield tmp.name


@pytest.fixture
def stats_registry(mock_stats_hdf5_file):
    """Mock stats registry using real StatsRegistry with mock file."""
    from pathlib import Path

    from data_services.stats import StatsRegistry

    # Create the selector dictionary in the format expected by StatsRegistry
    # {run_version: {ancestry: path}}
    selector = {
        "Run-v1": {"EUR": Path(mock_stats_hdf5_file)},
        "Run-v2": {"EUR": Path(mock_stats_hdf5_file)},
    }

    # Create a real StatsRegistry with our mock selector
    registry = StatsRegistry(selector)

    return registry


@pytest.fixture
def stats_service(stats_registry):
    """Get a stats service from the registry for tests that need it directly."""
    # Get the default service (Run-v1, EUR)
    return stats_registry.get("Run-v1", "EUR")


@pytest.fixture
def mock_nomaly_scores_hdf5_file():
    """Create a temporary mock HDF5 file with test nomaly scores data."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create required datasets

            # Individual IDs - should match those used in other fixtures
            test_eids = np.array(
                [1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010]
            )
            f.create_dataset("eid", data=test_eids)

            # Define some terms
            terms = np.array([b"GO:0030800", b"MP:0005179", b"GO:0006915"])
            f.create_dataset("term", data=terms)

            # Create score matrix (eids x terms)
            # Define scores that align with our test cases/controls
            scores = np.array(
                [
                    [0.03, 0.02, 0.01],  # 1001 (case)
                    [0.025, 0.015, 0.02],  # 1002 (case)
                    [0.01, 0.03, 0.015],  # 1003 (case)
                    [0.005, 0.01, 0.02],  # 1004 (control)
                    [0.015, 0.005, 0.01],  # 1005 (control)
                    [0.03, 0.015, 0.008],  # 1006 (control)
                    [0.01, 0.005, 0.002],  # 1007 (control)
                    [0.005, 0.01, 0.005],  # 1008 (control)
                    [0.0, 0.0, 0.0],  # 1009 (excluded)
                    [0.0, 0.0, 0.0],  # 1010 (excluded)
                ]
            )
            f.create_dataset("scores", data=scores)

        # Create the required .npy file
        np.save(f"{tmp.name}.npy", scores)

        yield Path(tmp.name)

        # Clean up the .npy file
        try:
            os.unlink(f"{tmp.name}.npy")
        except OSError:
            pass


@pytest.fixture
def nomaly_scores_service(mock_nomaly_scores_hdf5_file):
    """Mock nomaly scores service for unit tests."""
    from data_services.nomaly_score import NomalyScoreService

    # Create a real NomalyScoreService instance with our mock file
    service = NomalyScoreService(mock_nomaly_scores_hdf5_file)

    return service


@pytest.fixture
def nomaly_data_service():
    """Mock nomaly data service for unit tests."""
    from unittest.mock import MagicMock

    # Create a mock nomaly data service
    mock_service = MagicMock()
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
