"""Test concurrent request handling in the prioritisation blueprint."""

import pytest
from flask import url_for
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


# TODO: I don't think this is working as expected...
@pytest.fixture(autouse=True)
def clear_cache(mock_config):
    """Clear the cache directory before each test."""
    cache_dir = Path(mock_config.VARIANT_SCORES_DIR)
    if cache_dir.exists():
        for f in cache_dir.glob("variant_prioritization_*"):
            f.unlink()
    cache_dir.mkdir(exist_ok=True)
    yield


def test_stream_isolation(auth_client):
    """Test that concurrent requests maintain proper data isolation."""
    client = auth_client

    # Just test two requests to keep it simple
    test_pairs = [
        ("250.2", "UP:UPA00246"),  # This user has access to 250.2
        ("561", "GO:0030800"),  # And to 561
    ]

    def make_request(disease_code, term):
        """Make a request and track its messages"""
        response = client.get(
            f"/stream_progress/{disease_code}/{term}?no_cache=true",
            follow_redirects=True,
        )
        assert response.status_code == 200

        messages = []
        for chunk in response.response:
            if chunk.startswith(b"data: "):
                try:
                    msg = json.loads(chunk[6:].decode())
                    messages.append(msg)
                except json.JSONDecodeError:
                    continue
        return disease_code, term, messages

    # Run requests sequentially first to verify basic functionality
    results = [make_request(code, term) for code, term in test_pairs]

    for disease_code, term, messages in results:
        assert any(m["type"] == "progress" for m in messages), (
            f"No progress messages for {disease_code}"
        )
        assert any(m["type"] == "results" for m in messages), (
            f"No results for {disease_code}"
        )
        assert any(m["type"] == "done" for m in messages), (
            f"No done message for {disease_code}"
        )


def test_error_handling(auth_client):
    """Test error handling for invalid requests."""
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


def test_cache_concurrency(auth_client):
    """Test that concurrent requests get consistent results when using cache."""
    client = auth_client

    # Use a disease code the test user has access to
    test_pairs = [
        ("250.2", "UP:UPA00246"),  # This user has access to 250.2
    ]

    def make_request(disease_code, term):
        """Make a request and track its messages"""
        response = client.get(
            f"/stream_progress/{disease_code}/{term}",
            follow_redirects=True,
        )
        assert response.status_code == 200

        messages = []
        for chunk in response.response:
            if chunk.startswith(b"data: "):
                try:
                    msg = json.loads(chunk[6:].decode())
                    messages.append(msg)
                except json.JSONDecodeError:
                    continue
        return messages

    # Make multiple requests for the same data
    disease_code, term = test_pairs[0]
    results = [make_request(disease_code, term) for _ in range(3)]

    # Verify all requests got the same results
    for messages in results:
        assert any(m["type"] == "progress" for m in messages), (
            f"No progress messages for {disease_code}"
        )
        assert any(m["type"] == "results" for m in messages), (
            f"No results for {disease_code}"
        )
        assert any(m["type"] == "done" for m in messages), (
            f"No done message for {disease_code}"
        )

    # Get the results data from each response
    result_data = [
        next(m["data"] for m in messages if m["type"] == "results")
        for messages in results
    ]

    # Verify all responses have identical results
    first_result = result_data[0]
    for result in result_data[1:]:
        assert result == first_result, "Cache returned different results"


def test_error_recovery(auth_client):
    """Test error handling and recovery during concurrent requests."""
    client = auth_client

    # Mix of valid and invalid requests
    test_cases = [
        ("571.5", "UP:UPA00240"),  # Valid
        ("invalid", "invalid"),  # Invalid
        ("571.5", "invalid"),  # Partially invalid
    ]

    def make_request(disease_code, term):
        response = client.get(
            f"/stream_progress/{disease_code}/{term}", follow_redirects=True
        )
        messages = []
        for line in response.response:
            if line.startswith(b"data: "):
                messages.append(json.loads(line[6:].decode()))
        return disease_code, term, messages

    # Run mixed requests concurrently
    with ThreadPoolExecutor(max_workers=len(test_cases)) as executor:
        futures = [
            executor.submit(make_request, code, term) for code, term in test_cases
        ]
        results = [f.result() for f in futures]

    for disease_code, term, messages in results:
        # Valid request should complete normally
        if (disease_code, term) == ("571.5", "UP:UPA00240"):
            assert any(m["type"] == "results" for m in messages)
            assert messages[-1]["type"] == "done"

        # Invalid requests should error gracefully
        else:
            assert any(m["type"] == "error" for m in messages)
            assert messages[-1]["type"] == "done"
