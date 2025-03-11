import json

# TODO: WHY DOES THE unit_test_app_client NEED .../genotypes_with_metadata.h5 TO EXISTS?
# TODO: WHY DOES THE unit_test_app_client NEED .../phecode_cases_excludes_with_metadata.h5 TO EXISTS?
# TODO: WHY DOES THE unit_test_app_client NEED .../stats-All-2025-02-10.h5 TO EXISTS?
# TODO: WHY DOES THE unit_test_app_client NEED .../float16_scores.h5 TO EXISTS?
# NOTE: This is 'filed' as an INTEGRATION test?

def test_variant_api(unit_test_app_client):
    """Test the variant API endpoints."""
    test_variant = "11:69083946:T:C"

    # Test initiating PheWAS analysis
    response = unit_test_app_client.post(f"/run-phewas/{test_variant}")
    assert response.status_code == 202
    data = json.loads(response.data)
    assert data["status"] == "Task started"

    # Test getting PheWAS results
    response = unit_test_app_client.get(f"/phewas-result/{test_variant}")
    assert response.status_code == 200
    data = json.loads(response.data)
    assert "result" in data
    assert "associations" in data


def test_nomaly_api(auth_integration_app_client):
    """Test the Nomaly-related endpoints."""
    test_phecodes = ["250.2", "401.1", "296.2", "038", "008"]

    for test_phecode in test_phecodes:
        # Test getting Nomaly results
        response = auth_integration_app_client.get(f"/phecode/{test_phecode}")
        if response.status_code != 200:
            print(f"Unexpected response: {response.status_code}")
            print(
                f"Response data: {response.data}"
            )  # This will help debug the error message
        assert response.status_code == 200
        # Add specific assertions based on expected response structure

    test_phecode = "banana"
    response = auth_integration_app_client.get(f"/phecode/{test_phecode}")
    assert response.status_code == 404


# Plotting isn't called like this!
# def test_plot_generation(client):
#     """Test plot generation endpoints."""
#     test_phecode = '250.2'

#     response = client.get(f'/plot/{test_phecode}')
#     assert response.status_code == 200
#     # Add assertions for plot data structure in
