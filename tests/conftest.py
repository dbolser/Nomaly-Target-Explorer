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
def unit_test_app(test_app, mock_genotype_hdf5_file_with_npy, mock_phenotype_file, stats_service, stats_service_v2):
    """App configured for unit tests with mocked services."""

    # Override config paths with actual mock files
    test_app.config.update({
        "LOGIN_DISABLED": True,
        "GENOTYPES_H5": mock_genotype_hdf5_file_with_npy,
        "PHENOTYPES_H5": mock_phenotype_file,
    })

    with test_app.app_context():
        from data_services.genotype import GenotypeService
        from data_services.phenotype import PhenotypeService

        from data_services.stats import StatsService
        from data_services.nomaly_score import NomalyScoreService

        # Create test services
        services = ServiceRegistry()

        services.genotype = GenotypeService(mock_genotype_hdf5_file_with_npy)
        services.phenotype = PhenotypeService(mock_phenotype_file)

        # Use the mock stats services
        services.stats = stats_service
        services.stats_v2 = stats_service_v2
        services.nomaly_score = MagicMock()

        test_app.extensions["nomaly_services"] = services

        yield test_app


@pytest.fixture
def phenotype_service(unit_test_app):
    return unit_test_app.extensions["nomaly_services"].phenotype._hdf


@pytest.fixture
def genotype_service(unit_test_app):
    """Return the underlying GenotypesHDF5 object that has the individual property."""
    from unittest.mock import MagicMock
    import numpy as np
    
    # Create a mock service with the necessary methods
    mock_service = MagicMock()
    
    # Mock the individual property
    mock_service.individual = np.array([1001, 1002, 1003, 1004])
    
    # Mock the query_variantID_genotypes method
    def mock_query_variantID_genotypes(variant):
        # Return predefined results based on the variant
        if variant == "1_186977737_A/G":
            return np.array([1001, 1002, 1003, 1004]), np.array([2, 2, 2, 0])
        elif variant == "1_186977780_G/A":
            return np.array([1001, 1002, 1003, 1004]), np.array([-1, -1, 2, 1])
        elif variant == "1_46813503_C/T":
            return np.array([1001, 1002, 1003, 1004]), np.array([0, 1, 0, 2])
        else:
            return None
    
    mock_service.query_variantID_genotypes = mock_query_variantID_genotypes
    
    # Add any other methods that might be needed by tests
    mock_service.get_genotypes = MagicMock(return_value=np.zeros((3, 3)))
    
    return mock_service


@pytest.fixture
def stats_service():
    """Mock stats service for unit tests."""
    from unittest.mock import MagicMock
    import numpy as np
    import pandas as pd
    
    # Create a mock service with the necessary methods
    mock_service = MagicMock()
    
    # Create a mock _hdf attribute
    mock_hdf = MagicMock()
    
    # Set up properties on the mock HDF5 object
    mock_hdf.phecodes = np.array(["250.2", "290.11", "401.1", "635.2", "601.1"])
    mock_hdf.terms = np.array(["GO:0030800", "GO:0016403", "GO:0036265", "GO:0048712"])
    mock_hdf.statistics = np.array(["num_rp", "num_rn", "mwu_pvalue", "tti_pvalue", "metric1_pvalue"])
    
    # Sample 3D data array (terms x phecodes x statistics)
    mock_hdf.data = np.random.rand(4, 5, 5)  # 4 terms, 5 phecodes, 5 statistics
    
    # Mock get_stats_by_phecode method
    def mock_get_stats_by_phecode(phecode, statstype=None):
        if phecode not in mock_hdf.phecodes:
            raise ValueError(f"Disease {phecode} not found in dataset")
            
        # Find the index of the phecode
        phecode_idx = np.where(mock_hdf.phecodes == phecode)[0][0]
        
        if statstype is None:
            # Return all statistics for this phecode
            return pd.DataFrame(
                mock_hdf.data[:, phecode_idx, :],
                index=mock_hdf.terms.astype(str),
                columns=mock_hdf.statistics
            )
        elif isinstance(statstype, str):
            # Return specific statistic
            stat_idx = np.where(mock_hdf.statistics == statstype)[0][0]
            return mock_hdf.data[:, phecode_idx, stat_idx]
        else:
            # Return multiple statistics
            stat_indices = [np.where(mock_hdf.statistics == s)[0][0] for s in statstype]
            return mock_hdf.data[:, phecode_idx, stat_indices]
    
    # Mock get_stats_by_term_phecode method
    def mock_get_stats_by_term_phecode(term, phecode, statstype=None):
        if term not in mock_hdf.terms:
            raise ValueError(f"Term {term} not found in dataset")
        if phecode not in mock_hdf.phecodes:
            raise ValueError(f"Disease {phecode} not found in dataset")
            
        # Find indices
        term_idx = np.where(mock_hdf.terms == term)[0][0]
        phecode_idx = np.where(mock_hdf.phecodes == phecode)[0][0]
        
        if statstype is None:
            statstype = mock_hdf.statistics
        elif isinstance(statstype, str):
            statstype = np.array([statstype])
            
        # Create result dictionary
        result = {}
        for stat in statstype:
            stat_idx = np.where(mock_hdf.statistics == stat)[0][0]
            result[stat] = mock_hdf.data[term_idx, phecode_idx, stat_idx]
            
        return result
    
    # Assign methods to mock
    mock_hdf.get_stats_by_phecode = mock_get_stats_by_phecode
    mock_hdf.get_stats_by_term_phecode = mock_get_stats_by_term_phecode
    
    # Add get_stats_wip method if needed
    def mock_get_stats_wip(phecode=None, term=None, statstype=None):
        if not phecode and not term:
            raise ValueError("At least one of phecode or term must be provided!")
            
        # For simplicity, just return a DataFrame with sample data
        return pd.DataFrame(
            np.random.rand(5, 5),
            columns=mock_hdf.statistics
        )
    
    mock_hdf.get_stats_wip = mock_get_stats_wip
    
    # Assign the mock HDF5 object to the service
    mock_service._hdf = mock_hdf
    
    return mock_service


@pytest.fixture
def stats_service_v2():
    """Mock stats service V2 for unit tests."""
    # Reuse the same implementation as stats_service for now
    # In a real scenario, you might want to customize this for V2-specific behavior
    from unittest.mock import MagicMock
    import numpy as np
    import pandas as pd
    
    # Create a mock service with the necessary methods
    mock_service = MagicMock()
    
    # Create a mock _hdf attribute
    mock_hdf = MagicMock()
    
    # Set up properties on the mock HDF5 object
    mock_hdf.phecodes = np.array(["250.2", "290.11", "401.1", "635.2", "601.1"])
    mock_hdf.terms = np.array(["GO:0030800", "GO:0016403", "GO:0036265", "GO:0048712"])
    mock_hdf.statistics = np.array(["num_rp", "num_rn", "mwu_pvalue", "tti_pvalue", "metric1_pvalue"])
    
    # Sample 3D data array (terms x phecodes x statistics)
    mock_hdf.data = np.random.rand(4, 5, 5)  # 4 terms, 5 phecodes, 5 statistics
    
    # Mock get_stats_by_phecode method
    def mock_get_stats_by_phecode(phecode, statstype=None):
        if phecode not in mock_hdf.phecodes:
            raise ValueError(f"Disease {phecode} not found in dataset")
            
        # Find the index of the phecode
        phecode_idx = np.where(mock_hdf.phecodes == phecode)[0][0]
        
        if statstype is None:
            # Return all statistics for this phecode
            return pd.DataFrame(
                mock_hdf.data[:, phecode_idx, :],
                index=mock_hdf.terms.astype(str),
                columns=mock_hdf.statistics
            )
        elif isinstance(statstype, str):
            # Return specific statistic
            stat_idx = np.where(mock_hdf.statistics == statstype)[0][0]
            return mock_hdf.data[:, phecode_idx, stat_idx]
        else:
            # Return multiple statistics
            stat_indices = [np.where(mock_hdf.statistics == s)[0][0] for s in statstype]
            return mock_hdf.data[:, phecode_idx, stat_indices]
    
    # Mock get_stats_by_term_phecode method
    def mock_get_stats_by_term_phecode(term, phecode, statstype=None):
        if term not in mock_hdf.terms:
            raise ValueError(f"Term {term} not found in dataset")
        if phecode not in mock_hdf.phecodes:
            raise ValueError(f"Disease {phecode} not found in dataset")
            
        # Find indices
        term_idx = np.where(mock_hdf.terms == term)[0][0]
        phecode_idx = np.where(mock_hdf.phecodes == phecode)[0][0]
        
        if statstype is None:
            statstype = mock_hdf.statistics
        elif isinstance(statstype, str):
            statstype = np.array([statstype])
            
        # Create result dictionary
        result = {}
        for stat in statstype:
            stat_idx = np.where(mock_hdf.statistics == stat)[0][0]
            result[stat] = mock_hdf.data[term_idx, phecode_idx, stat_idx]
            
        return result
    
    # Assign methods to mock
    mock_hdf.get_stats_by_phecode = mock_get_stats_by_phecode
    mock_hdf.get_stats_by_term_phecode = mock_get_stats_by_term_phecode
    
    # Add get_stats_wip method if needed
    def mock_get_stats_wip(phecode=None, term=None, statstype=None):
        if not phecode and not term:
            raise ValueError("At least one of phecode or term must be provided!")
            
        # For simplicity, just return a DataFrame with sample data
        return pd.DataFrame(
            np.random.rand(5, 5),
            columns=mock_hdf.statistics
        )
    
    mock_hdf.get_stats_wip = mock_get_stats_wip
    
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
