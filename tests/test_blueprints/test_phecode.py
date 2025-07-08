"""Tests for the phecode blueprint and related functionality."""

import pytest
from flask import url_for

from tests.data_utils import production_data_available

pytestmark = pytest.mark.skipif(
    not production_data_available(),
    reason="requires production data services",
)


def test_home_route_works(auth_integration_app_client):
    """Test that the home page is accessible when logged in."""
    # We can see from our debug output that the home page is accessible
    response = auth_integration_app_client.get(url_for("index"))

    # This page should work
    assert response.status_code == 200, "Home page should be accessible"

    # We should see some expected content on the home page
    assert b"Nomaly" in response.data, "Home page should contain 'Nomaly'"


def test_search_route_works(auth_integration_app_client):
    """Test that the search page is accessible when logged in."""
    # We can see from our debug output that the search page is accessible
    response = auth_integration_app_client.get(url_for("search.show"))

    # This page should work
    assert response.status_code == 200, "Search page should be accessible"

    # We should see some expected content on the search page
    assert b"search" in response.data.lower(), "Search page should contain 'search'"


def test_admin_phecode_access(auth_integration_app_client):
    """Test that an admin can access phecode pages."""
    # Try a phecode that should exist in the test data
    phecode = "250.2"  # Type 2 diabetes

    # The admin should be able to access the phecode page directly
    response = auth_integration_app_client.get(
        url_for("phecode.show_phecode", phecode=phecode),
        follow_redirects=False,
    )

    # Print the result for debugging
    print(f"\nAdmin access to phecode {phecode}: {response.status_code}")
    if response.status_code == 302:
        redirect_to = response.headers.get("Location", "unknown")
        print(f"  Redirects to: {redirect_to}")

    # If the admin has proper permissions, this should work without redirect
    assert response.status_code == 200, "Admin should be able to access phecode pages"

    # We should see the phecode on the page
    assert bytes(phecode, "utf-8") in response.data, (
        f"Should show phecode {phecode} content"
    )

    # Now try a term page as well
    term = "GO:0030800"
    term_response = auth_integration_app_client.get(
        url_for("phecode_term.show_phecode_term", phecode=phecode, term=term),
        follow_redirects=False,
    )

    # Print the result for debugging
    print(f"Admin access to phecode/term {phecode}/{term}: {term_response.status_code}")

    # This should also be accessible
    assert term_response.status_code == 200, (
        "Admin should be able to access phecode term pages"
    )

    # We should see both the phecode and term on the page
    assert bytes(phecode, "utf-8") in term_response.data
    assert bytes(term, "utf-8") in term_response.data


@pytest.mark.skip(reason="We need to set up a limited user to test this.")
def test_phecode_route_permission(auth_integration_app_client):
    """Test that phecode routes redirect for users without permissions."""
    # Using the phecode from our tests
    phecode = "250.2"

    # The phecode page should redirect
    response = auth_integration_app_client.get(
        url_for("phecode.show_phecode", phecode=phecode),
        follow_redirects=False,
    )

    # Check that we get a redirect
    assert response.status_code == 302, "Expected redirect for restricted phecode page"

    # Check that the redirect location is the home page
    assert response.headers.get("Location") == "/", "Should redirect to home page"

    # Now follow the redirect and check for the error message
    response = auth_integration_app_client.get(
        url_for("phecode.show_phecode", phecode=phecode),
        follow_redirects=True,
    )
    assert response.status_code == 200

    # We should see an access denied message
    assert (
        b"access denied" in response.data.lower()
        or b"insufficient permissions" in response.data.lower()
        or b"not authorized" in response.data.lower()
    ), "Should show access denied message"
