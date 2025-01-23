"""Test concurrent request handling in the prioritisation blueprint."""

import pytest
from flask import url_for, current_app
import threading
import queue
import json
import atexit
from concurrent.futures import ThreadPoolExecutor, as_completed
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

    # We'll use two requests: one slow (more variants) and one fast (fewer variants)
    test_pairs = [
        ("571.5", "UP:UPA00240"),  # potentially slower
        ("250.2", "UP:UPA00246"),  # potentially faster
    ]

    results = queue.Queue()

    def process_request(disease_code, term):
        """Process a single request and collect its output stream."""
        with app.app_context():
            response = client.get(
                url_for(
                    "prioritisation.stream_progress",
                    disease_code=disease_code,
                    term=term,
                )
            )

            assert response.status_code == 200

            # Collect all SSE messages with timestamps
            messages = []
            for line in response.data.decode().split("\n\n"):
                if line.startswith("data: "):
                    message = line[6:]  # Remove 'data: ' prefix
                    messages.append(message)

            results.put(
                {
                    "disease_code": disease_code,
                    "term": term,
                    "messages": messages,
                }
            )

    # Start all requests concurrently
    with ThreadPoolExecutor(max_workers=len(test_pairs)) as executor:
        futures = [
            executor.submit(process_request, disease_code, term)
            for disease_code, term in test_pairs
        ]
        for future in as_completed(futures):
            future.result()

    # Verify results
    processed_results = []
    while not results.empty():
        processed_results.append(results.get())

    assert len(processed_results) == len(test_pairs)

    # For each stream, verify it only contains messages relevant to its disease/term
    for result in processed_results:
        disease_code = result["disease_code"]
        term = result["term"]
        messages = result["messages"]

        # Check that each non-DONE message mentions the correct disease code or term
        for msg in messages:
            if msg != "DONE":
                other_pairs = [p for p in test_pairs if p != (disease_code, term)]
                for other_disease, other_term in other_pairs:
                    assert other_disease not in msg, (
                        f"Found disease {other_disease} in stream for {disease_code}"
                    )
                    assert other_term not in msg, (
                        f"Found term {other_term} in stream for {term}"
                    )


@pytest.mark.usefixtures("app")
def test_concurrent_requests(auth_client, configure_app):
    """Test that multiple concurrent requests don't interfere with each other."""
    app = configure_app  # Use the configured app
    client = auth_client  # Use the authenticated client

    # Different disease/term pairs to test
    test_pairs = [
        ("571.5", "UP:UPA00240"),
        ("332", "GO:0030800"),
        ("250.2", "UP:UPA00246"),
    ]

    results = queue.Queue()

    def process_request(disease_code, term):
        """Process a single request and collect its output stream."""
        with app.app_context():
            # Start the processing
            response = client.get(
                url_for(
                    "prioritisation.stream_progress",
                    disease_code=disease_code,
                    term=term,
                )
            )

            assert response.status_code == 200

            # Collect all SSE messages
            messages = []
            for line in response.data.decode().split("\n\n"):
                if line.startswith("data: "):
                    message = line[6:]  # Remove 'data: ' prefix
                    messages.append(message)
                    if message == "DONE":
                        break

            # Get the results
            results_response = client.get(
                url_for(
                    "prioritisation.get_variant_scores_data",
                    disease_code=disease_code,
                    term=term,
                )
            )

            assert results_response.status_code == 200
            results_data = json.loads(results_response.data)

            # Store results for verification
            results.put(
                {
                    "disease_code": disease_code,
                    "term": term,
                    "messages": messages,
                    "results": results_data,
                }
            )

    # Start all requests concurrently
    with ThreadPoolExecutor(max_workers=len(test_pairs)) as executor:
        futures = [
            executor.submit(process_request, disease_code, term)
            for disease_code, term in test_pairs
        ]
        # Wait for all to complete
        for future in as_completed(futures):
            # This will raise any exceptions that occurred
            future.result()

    # Verify results
    processed_results = []
    while not results.empty():
        processed_results.append(results.get())

    assert len(processed_results) == len(test_pairs)

    # Verify each result has its own unique stream and correct data
    for result in processed_results:
        # Check we got some output messages
        assert len(result["messages"]) > 0
        # Check the last message is DONE
        assert result["messages"][-1] == "DONE"
        # Verify we got results data
        assert "top_variants" in result["results"]
        assert "top_gene_set" in result["results"]

        # Verify results match the request
        disease_term_pair = (result["disease_code"], result["term"])
        assert disease_term_pair in test_pairs


@pytest.mark.usefixtures("app")
def test_error_handling_concurrent(auth_client, configure_app):
    """Test error handling during concurrent processing."""
    app = configure_app
    client = auth_client

    # Verify we're actually logged in
    response = client.get("/")
    assert response.status_code == 200, "Not logged in at test start"

    # Mix of valid and invalid requests
    test_pairs = [
        ("571.5", "UP:UPA00240"),  # valid
        ("999.9", "INVALID"),  # invalid
        ("332", "GO:0030800"),  # valid
    ]

    results = queue.Queue()

    def process_request(disease_code, term):
        """Process a single request and collect its output stream."""
        with app.app_context():
            # Use the same session for all requests
            with client.session_transaction() as sess:
                print(f"Session before request: {sess}")

            response = client.get(
                url_for(
                    "prioritisation.stream_progress",
                    disease_code=disease_code,
                    term=term,
                ),
                follow_redirects=True,  # Follow any redirects to help debug
            )

            if response.status_code != 200:
                print(
                    f"Unexpected status code {response.status_code} for {disease_code}, {term}"
                )
                print(f"Response data: {response.data.decode()}")

            # Collect all messages from the stream
            messages = []
            try:
                for line in response.data.decode().split("\n\n"):
                    if line.startswith("data: "):
                        message = line[6:]  # Remove 'data: ' prefix
                        messages.append(message)
                        # Don't break on DONE - collect all messages
            except Exception as e:
                print(f"Error processing stream for {disease_code}, {term}: {e}")
                messages.append(f"ERROR: Stream processing failed: {str(e)}")

            # Get the results regardless of stream status
            results_response = client.get(
                url_for(
                    "prioritisation.get_variant_scores_data",
                    disease_code=disease_code,
                    term=term,
                ),
                follow_redirects=True,
            )

            # Check session after requests
            with client.session_transaction() as sess:
                print(f"Session after requests: {sess}")

            results.put(
                {
                    "disease_code": disease_code,
                    "term": term,
                    "messages": messages,
                    "results_status": results_response.status_code,
                    "stream_status": response.status_code,
                    "redirect_chain": getattr(response, "history", []),
                }
            )

    # Process all requests concurrently
    with ThreadPoolExecutor(max_workers=len(test_pairs)) as executor:
        futures = [
            executor.submit(process_request, disease_code, term)
            for disease_code, term in test_pairs
        ]
        for future in as_completed(futures):
            future.result()

    # Verify results
    processed_results = []
    while not results.empty():
        processed_results.append(results.get())

    assert len(processed_results) == len(test_pairs)

    # Check that valid requests succeeded and invalid ones failed appropriately
    for result in processed_results:
        print(f"Checking result for {result['disease_code']}, {result['term']}")
        print(f"Messages: {result['messages']}")
        print(
            f"Status codes - Stream: {result['stream_status']}, Results: {result['results_status']}"
        )
        if result.get("redirect_chain"):
            print(f"Redirect chain: {[r.location for r in result['redirect_chain']]}")

        if result["disease_code"] == "999.9" or result["term"] == "INVALID":
            # For invalid requests:
            # - Either we should see an error message in the stream
            # - Or the results endpoint should return 404
            # - Or both
            has_error = any("ERROR" in msg for msg in result["messages"])
            is_404 = result["results_status"] == 404
            assert has_error or is_404, (
                f"Invalid request didn't fail properly: {result}"
            )
        else:
            # For valid requests:
            # - Stream should complete with DONE
            # - Results should be available
            assert result["stream_status"] == 200, (
                f"Stream failed for valid request: {result}"
            )
            assert result["results_status"] == 200, (
                f"Results not available for valid request: {result}"
            )
            if not any(msg == "DONE" for msg in result["messages"]):
                print(
                    f"Warning: No DONE message for valid request {result['disease_code']}, {result['term']}"
                )
                print(f"Messages received: {result['messages']}")
                assert False, "Valid request didn't complete properly"
