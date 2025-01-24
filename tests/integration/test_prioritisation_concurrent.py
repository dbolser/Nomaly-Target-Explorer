"""Test concurrent request handling in the prioritisation blueprint."""

import pytest
from flask import url_for
import json
from concurrent.futures import ThreadPoolExecutor
from db import get_db_connection


@pytest.fixture(autouse=True)
def cleanup_test_users():
    """Ensure test users are cleaned up even if tests fail or hang."""
    yield
    # Force cleanup of any remaining test users after each test
    conn = get_db_connection()
    with conn.cursor() as cursor:
        cursor.execute(
            "DELETE FROM user_permissions WHERE user_id IN (SELECT id FROM users2 WHERE username = 'test_admin')"
        )
        cursor.execute("DELETE FROM users2 WHERE username = 'test_admin'")
        conn.commit()
    conn.close()


@pytest.fixture
def configure_app(app):
    """Configure the app for URL generation in tests."""
    app.config.update(
        {
            "SERVER_NAME": "localhost.localdomain",
            "PREFERRED_URL_SCHEME": "http",
        }
    )
    return app


@pytest.mark.usefixtures("app")
def test_stream_isolation(auth_client, configure_app):
    """Test that streams are properly isolated and don't leak between requests."""
    app = configure_app
    client = auth_client

    # Verify we're logged in first
    response = client.get("/")
    assert response.status_code == 200, "Not logged in at start of test"

    # Test pairs with different characteristics
    test_pairs = [
        ("571.5", "UP:UPA00240"),  # potentially slower
        ("250.2", "UP:UPA00246"),  # potentially faster
    ]

    def collect_messages(response):
        """Collect messages from a streaming response."""
        messages = []
        for line in response.response:  # Flask test client provides an iterator
            if line.startswith(b"data: "):
                message = json.loads(line[6:].decode())
                messages.append(message)
                # Stop if we get a done message
                if message.get("type") == "done":
                    break
        return messages

    def run_request(disease_code, term):
        """Run a single request and collect its messages."""
        with app.app_context():
            response = client.get(
                url_for(
                    "prioritisation.stream_progress",
                    disease_code=disease_code,
                    term=term,
                ),
                follow_redirects=True,
            )
            assert response.status_code == 200
            return {
                "disease_code": disease_code,
                "term": term,
                "messages": collect_messages(response),
            }

    # Run requests concurrently
    with ThreadPoolExecutor(max_workers=len(test_pairs)) as executor:
        futures = [
            executor.submit(run_request, code, term) for code, term in test_pairs
        ]
        results = [f.result() for f in futures]

    # Verify each stream's integrity
    for result in results:
        messages = result["messages"]

        # Check message sequence
        types = [m["type"] for m in messages]
        assert "progress" in types, "Should have progress messages"
        assert "results" in types, "Should have results"
        assert "done" in types, "Should end with done"

        # Check data isolation
        progress_msgs = [m for m in messages if m["type"] == "progress"]
        for msg in progress_msgs:
            # Verify this stream only contains its own data
            assert result["disease_code"] in str(msg) or result["term"] in str(msg), (
                "Progress message should reference its own request"
            )


@pytest.mark.usefixtures("app")
def test_error_handling(auth_client, configure_app):
    """Test error handling for invalid requests."""
    #app = configure_app
    client = auth_client

    # Invalid request
    response = client.get(
        url_for(
            "prioritisation.stream_progress", disease_code="invalid", term="invalid"
        ),
        follow_redirects=True,
    )

    messages = []
    for line in response.response:
        if line.startswith(b"data: "):
            messages.append(json.loads(line[6:].decode()))

    # Should have error and done messages
    types = [m["type"] for m in messages]
    assert "error" in types, "Should have error message"
    assert "done" in types, "Should end with done message"
