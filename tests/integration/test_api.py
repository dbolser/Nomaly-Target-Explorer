import json

# the client is defined in tests/conftest.py


def test_variant_api(client):
    """Test the variant API endpoints."""
    test_variant = "11:69083946:T:C"

    # Test initiating PheWAS analysis
    response = client.post(f"/run-phewas/{test_variant}")
    assert response.status_code == 202
    data = json.loads(response.data)
    assert data["status"] == "Task started"

    # Test getting PheWAS results
    response = client.get(f"/phewas-result/{test_variant}")
    assert response.status_code == 200
    data = json.loads(response.data)
    assert "result" in data
    assert "associations" in data


def test_nomaly_api(client):
    """Test the Nomaly-related endpoints."""
    test_phecode = "250.2"

    # Test getting Nomaly results
    response = client.get(f"/phecode/{test_phecode}")
    if response.status_code != 200:
        print(f"Unexpected response: {response.status_code}")
        print(
            f"Response data: {response.data}"
        )  # This will help debug the error message
    assert response.status_code == 200
    # Add specific assertions based on expected response structure


# Plotting isn't called like this!
# def test_plot_generation(client):
#     """Test plot generation endpoints."""
#     test_phecode = '250.2'

#     response = client.get(f'/plot/{test_phecode}')
#     assert response.status_code == 200
#     # Add assertions for plot data structure in
