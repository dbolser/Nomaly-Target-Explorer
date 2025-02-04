from unittest.mock import patch, MagicMock
import pytest
from werkzeug.security import generate_password_hash
import numpy as np
import h5py
import tempfile
from pathlib import Path

from app import create_app
from services import ServiceRegistry
from db import get_db_connection
import os


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
    _app = create_app()  # NOTE: Services get create as normal here!
    _app.config.update(
        {
            "TESTING": True,
            "SERVER_NAME": "localhost.localdomain",
            "PREFERRED_URL_SCHEME": "http",
            "APPLICATION_ROOT": "/",
        }
    )
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
                    b"1:186977780:G:A",
                    b"1:46813503:C:T",
                ]
            )
            hdf.create_dataset("bim", data=bim)

            # Use underscore format for nomaly_variant_ids
            nomaly_variant_id = np.array(
                [
                    b"1_100_A/T",
                    b"2_200_C/G",
                    b"1_186977737_A/G",
                    b"1_186977780_G/A",
                    b"1_46813503_C/T",
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
def unit_test_app(test_app, mock_genotype_hdf5_file_with_npy, mock_phenotype_file):
    """App configured for unit tests with mocked services."""

    test_app.config["LOGIN_DISABLED"] = True

    with test_app.app_context():
        from data_services.genotype import GenotypeService
        from data_services.phenotype import PhenotypeService

        # Create test services
        services = ServiceRegistry()

        # TODO: Mock these out using fixtures from specific tests.
        services.genotype = GenotypeService(mock_genotype_hdf5_file_with_npy)
        services.phenotype = PhenotypeService(mock_phenotype_file)
        services.stats = MagicMock()
        services.stats_v2 = MagicMock()

        test_app.extensions["nomaly_services"] = services

        yield test_app


@pytest.fixture
def phenotype_service(unit_test_app):
    return unit_test_app.extensions["nomaly_services"].phenotype._hdf


@pytest.fixture
def unit_test_app_client(unit_test_app):
    """Client for unit tests."""
    return unit_test_app.test_client()


@pytest.fixture(scope="session")
def integration_app():
    """App configured for integration tests with real services."""
    _app = create_app()
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
