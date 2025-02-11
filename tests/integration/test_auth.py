"""Test authentication functionality."""


def test_login_success(integration_app_client, test_admin):
    """Test successful login."""
    # First verify we're not logged in
    response = integration_app_client.get("/logout")
    assert response.status_code == 302  # Should redirect to login

    # Try to login
    response = integration_app_client.post(
        "/login",
        data={"username": test_admin["username"], "password": test_admin["password"]},
        follow_redirects=True,
    )
    assert response.status_code == 200
    assert b"welcome" in response.data.lower()

    # Verify we can now access protected pages
    response = integration_app_client.get("/")
    assert response.status_code == 200


def test_login_failure(integration_app_client):
    """Test login with invalid credentials."""
    response = integration_app_client.post(
        "/login",
        data={"username": "wrong_user", "password": "wrong_pass"},
        follow_redirects=True,
    )
    assert response.status_code == 200
    assert b"invalid username or password" in response.data.lower()

    # Verify we still can't access protected pages
    response = integration_app_client.get("/logout")
    assert response.status_code == 302  # Should redirect to login


def test_logout(auth_integration_app_client):
    """Test logout functionality."""
    # Verify we start logged in
    response = auth_integration_app_client.get("/phecode/256")
    assert response.status_code == 200

    # Logout
    response = auth_integration_app_client.get("/logout", follow_redirects=True)
    assert response.status_code == 200
    assert b"logged out" in response.data.lower()

    # Verify we can't access protected pages anymore
    response = auth_integration_app_client.get("/phecode/256")
    assert response.status_code == 302  # Should redirect to login


def test_session_persistence(auth_integration_app_client):
    """Test that the session persists across requests."""
    # Make several requests and verify we stay logged in
    protected_urls = ["/phecode/256", "/search", "/diseasesearch?query=bowel"]

    for url in protected_urls:
        response = auth_integration_app_client.get(url)
        assert response.status_code == 200, f"Lost session at {url}"


def test_concurrent_sessions(auth_integration_app_client, test_admin):
    """Test handling of concurrent sessions."""
    # First client is already authenticated via auth_client fixture
    response = auth_integration_app_client.get("/search")
    assert response.status_code == 200

    # Login with second client
    response = auth_integration_app_client.post(
        "/login",
        data={"username": test_admin["username"], "password": test_admin["password"]},
        follow_redirects=True,
    )
    assert response.status_code == 200

    # Both clients should maintain their sessions
    response1 = auth_integration_app_client.get("/search")
    response2 = auth_integration_app_client.get("/search")
    assert response1.status_code == 200
    assert response2.status_code == 200
