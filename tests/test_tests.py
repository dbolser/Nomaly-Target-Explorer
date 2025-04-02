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
        assert services.nomaly_score._hdf is not None

        # Check that our mock services have the expected methods
        assert hasattr(services.phenotype._hdf, "get_cases_for_phecode")
        assert hasattr(services.genotype._hdf, "query_variantID_genotypes")
        assert hasattr(services.genotype._hdf, "get_genotypes")
        assert hasattr(services.nomaly_score._hdf, "get_scores_by_eids_unsorted")


def test_unit_test_app_client(unit_test_app_client):
    """Test that unit_test_app_client provides a working test client."""
    # With LOGIN_DISABLED=True, we should be able to access protected routes
    response = unit_test_app_client.get("/search")
    assert response.status_code == 200  # Should work without auth due to LOGIN_DISABLED


def test_mock_genotype_hdf5_file_with_npy(mock_genotype_hdf5_file_with_npy):
    """Test that the mock_genotype_hdf5_file_with_npy fixture creates a valid HDF5 file."""
    assert mock_genotype_hdf5_file_with_npy.exists()
    assert mock_genotype_hdf5_file_with_npy.with_suffix(".npy").exists()


def test_stats_registry(stats_registry):
    """Test that the stats_registry fixture provides a working registry."""
    # Check that the registry is initialized
    assert stats_registry.initialized is True

    # Check that we can get services for different versions
    service_v1 = stats_registry.get("Run-v1", "EUR")
    assert service_v1 is not None
    assert hasattr(service_v1, "_hdf")

    service_v2 = stats_registry.get("Run-v2", "EUR")
    assert service_v2 is not None
    assert hasattr(service_v2, "_hdf")


def test_stats_service(stats_service):
    """Test that the stats_service fixture provides a working service."""
    # Verify basic properties
    assert stats_service is not None
    assert hasattr(stats_service, "_hdf")

    # Test accessing the test data through the service
    df = stats_service.get_phecode_stats("250.2", term="GO:0030800")
    assert df is not None
    assert not df.empty
    assert "num_rp" in df.columns
    assert df["num_rp"].iloc[0] == 3.0

    # Test the term_stats function
    df2 = stats_service.get_term_stats("GO:0030800", phecode="250.2")
    assert df2 is not None
    assert not df2.empty
    assert "num_rp" in df2.columns
    assert df2["num_rp"].iloc[0] == 3.0


def test_phenotype_service(phenotype_service):
    """Test that the phenotype_service fixture provides a working service."""
    # Verify basic properties
    assert phenotype_service is not None
    assert hasattr(phenotype_service, "_hdf")

    # Test the real service implementations using our mock data
    df = phenotype_service.get_cases_for_phecode("250.2")
    assert df is not None
    assert not df.empty
    assert "eid" in df.columns
    assert "sex" in df.columns
    assert "phenotype" in df.columns


def test_genotype_service(genotype_service):
    """Test that the genotype_service fixture provides a working service."""
    # Verify basic properties
    assert genotype_service is not None
    assert hasattr(genotype_service, "_hdf")

    # Check we can access the individuals
    assert hasattr(genotype_service, "individual")
    assert len(genotype_service.individual) > 0


def test_nomaly_scores_service(nomaly_scores_service):
    """Test that the nomaly_scores_service fixture provides a working service."""
    # Verify basic properties
    assert nomaly_scores_service is not None
    assert hasattr(nomaly_scores_service, "_hdf")

    # Test the real service implementation
    import numpy as np

    test_eids = np.array([1001, 1002, 1003])
    scores = nomaly_scores_service.get_scores_by_eids_unsorted(test_eids)
    assert scores is not None
    assert len(scores) == len(test_eids)
    assert np.array_equal(scores[0], np.array([0.030, 0.020, 0.010, 0.001]))
    assert np.array_equal(scores[1], np.array([0.025, 0.015, 0.020, 0.002]))

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
        assert hasattr(services, "nomaly_score")
        assert hasattr(services, "nomaly_data")

        assert hasattr(services, "stats_registry")

        # Check that each service has the 'initalised' attribute
        assert services.genotype._check_initialized
        assert services.phenotype._check_initialized
        assert services.nomaly_score._check_initialized
        assert services.nomaly_data._check_initialized
        assert services.stats_registry._check_initialized


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
    # Try different routes and print the results
    test_routes = [
        "/",  # Home
        "/search",  # Search page
        "/phecode/250.2",  # Phecode from test data
        "/phecode/256",  # Phecode from test_logout
        "/phenotype",  # Phenotype page
        "/profile",  # Profile page
        "/variant/1_100_A/T",  # Variant page
    ]

    print("\nAccessibility of routes in actual test:")
    for route in test_routes:
        response = auth_integration_app_client.get(route, follow_redirects=False)
        status = response.status_code
        if status == 302:
            redirect_to = response.headers.get("Location", "unknown")
            print(f"  {route}: {status} -> redirects to {redirect_to}")
        else:
            print(f"  {route}: {status}")

    # Verify we can at least access the home page
    home_response = auth_integration_app_client.get("/")
    assert home_response.status_code == 200, "Could not access home page"

    # Check session state
    with auth_integration_app_client.session_transaction() as session:
        print(f"Session in test: {session}")
        # Check there's something in the session
        assert len(session) > 0, "Empty session - not logged in"


def test_integration_services_initialized(auth_integration_app_client):
    """Test that the real data services are properly initialized."""
    with auth_integration_app_client.application.app_context():
        services = auth_integration_app_client.application.extensions.get(
            "nomaly_services"
        )
        assert services is not None, "No services registered in the app"

        # Check if real data services are initialized
        assert services.phenotype.initialized, "Phenotype service not initialized"
        assert services.genotype.initialized, "Genotype service not initialized"
        assert services.nomaly_score.initialized, "Nomaly score service not initialized"
        assert services.stats_registry.initialized, "Stats registry not initialized"

        # Verify we can access real data
        phecode_df = services.phenotype.get_cases_for_phecode("250.2")
        assert phecode_df is not None, "Could not get phenotype data"
        assert len(phecode_df) > 0, "Phenotype data is empty"


def test_mock_config(mock_cache_dir_config):
    """Test that mock_config provides a temporary directory."""
    assert os.path.exists(mock_cache_dir_config.VARIANT_SCORES_DIR)
    assert os.path.isdir(mock_cache_dir_config.VARIANT_SCORES_DIR)


def test_debug_accessible_routes(auth_integration_app_client):
    """Debug test to check which routes are accessible when logged in."""
    # Try different routes and print the results
    test_routes = [
        "/",  # Home
        "/search",  # Search page
        "/phecode/250.2",  # Phecode from test data
        "/phecode/256",  # Phecode from test_logout
        "/phenotype",  # Phenotype page
        "/profile",  # Profile page
        "/variant/1_100_A/T",  # Variant page
    ]

    print("\nChecking accessibility of routes:")
    for route in test_routes:
        response = auth_integration_app_client.get(route, follow_redirects=False)
        status = response.status_code
        if status == 302:
            redirect_to = response.headers.get("Location", "unknown")
            print(f"  {route}: {status} -> redirects to {redirect_to}")
        else:
            print(f"  {route}: {status}")

    # Also check if the login succeeded at all
    with auth_integration_app_client.session_transaction() as session:
        print(f"Session contents: {session}")
