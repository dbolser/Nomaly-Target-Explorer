

def test_phecode_term_route_works(auth_integration_app_client):
    """Test that the phecode term route works for a known good example."""
    phecode = "705"
    term = "CC:MESH:C020806"

    response = auth_integration_app_client.get(f"/phecode/{phecode}/term/{term}")

    # This page should work - if it returns 500, the test will fail
    assert response.status_code == 200

    # We should see some expected content in the response
    assert bytes(phecode, "utf-8") in response.data
    assert bytes(term, "utf-8") in response.data
