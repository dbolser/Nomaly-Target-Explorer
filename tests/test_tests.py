from flask import current_app
import os


def test_test_app(test_app):
    """Test that test_app fixture provides a working Flask app."""
    assert test_app.config["TESTING"]


def test_unit_test_app(unit_test_app):
    """Test that unit_test_app provides an app with mocked services."""
    with unit_test_app.app_context():
        # Check config
        assert current_app.config["TESTING"]
        assert current_app.config[
            "LOGIN_DISABLED"
        ]  # This should be True for unit tests

        # Check services
        assert "nomaly_services" in current_app.extensions
        services = current_app.extensions["nomaly_services"]

        # Check that all services are properly mocked
        assert services.genotype._hdf is not None
        assert services.phenotype._hdf is not None
        assert services.stats._hdf is not None
        assert services.nomaly_score._hdf is not None

        # Check that our mock services have the expected methods
        assert hasattr(services.phenotype._hdf, "get_cases_for_phecode")
        assert hasattr(services.genotype._hdf, "query_variantID_genotypes")
        assert hasattr(services.genotype._hdf, "get_genotypes")
        assert hasattr(services.stats._hdf, "get_stats_by_term_phecode")
        assert hasattr(services.nomaly_score._hdf, "get_scores_by_eids_unsorted")


def test_unit_test_app_client(unit_test_app_client):
    """Test that unit_test_app_client provides a working test client."""
    # With LOGIN_DISABLED=True, we should be able to access protected routes
    response = unit_test_app_client.get("/search")
    assert response.status_code == 200  # Should work without auth due to LOGIN_DISABLED


def test_mock_genotype_hdf5_file_with_npy(mock_genotype_hdf5_file_with_npy):
    """Test that the mock_genotype_hdf5_file_with_npy fixture creates a valid HDF5 file."""
    assert os.path.exists(mock_genotype_hdf5_file_with_npy)
    assert mock_genotype_hdf5_file_with_npy.endswith(".h5")
    assert os.path.exists(f"{mock_genotype_hdf5_file_with_npy}.npy")


def test_integration_app(integration_app):
    """Test that integration_app provides an app with real services."""
    with integration_app.app_context():
        # Check that the app is configured for development (not testing)
        assert integration_app.config["TESTING"] is False
        assert integration_app.config["DEBUG"] is True  # Development mode

        # Check that services are manually initialized
        assert "nomaly_services" in current_app.extensions
        services = current_app.extensions["nomaly_services"]

        # Check that all required services exist
        assert hasattr(services, "genotype")
        assert hasattr(services, "phenotype")
        assert hasattr(services, "stats")
        assert hasattr(services, "nomaly_score")

        # Check that each service has the _hdf attribute
        assert hasattr(services.genotype, "_hdf")
        assert hasattr(services.phenotype, "_hdf")
        # assert hasattr(services.stats, "_hdf")
        assert hasattr(services.nomaly_score, "_hdf")



def test_integration_app_client(integration_app_client):
    """Test that integration_app_client provides a working test client."""
    response = integration_app_client.get("/search")
    assert response.status_code == 302  # Redirect to login


def test_test_admin(test_admin):
    """Test that test_admin fixture creates a valid admin user."""
    assert "id" in test_admin
    assert "username" in test_admin
    assert "password" in test_admin
    assert test_admin["username"] == "test_admin"


def test_auth_integration_app_client(auth_integration_app_client):
    """Test that auth_integration_app_client provides an authenticated client."""
    response = auth_integration_app_client.get("/search")
    assert response.status_code == 200  # Should be authenticated


def test_mock_config(mock_cache_dir_config):
    """Test that mock_config provides a temporary directory."""
    assert os.path.exists(mock_cache_dir_config.VARIANT_SCORES_DIR)
    assert os.path.isdir(mock_cache_dir_config.VARIANT_SCORES_DIR)
