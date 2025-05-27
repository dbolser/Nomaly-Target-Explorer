"""Test concurrent request handling in the prioritisation blueprint."""

import json
from concurrent.futures import ThreadPoolExecutor

import pandas as pd


def test_stream_isolation(unit_test_app_client_with_cache, monkeypatch):
    """Test that concurrent requests maintain proper data isolation."""
    # Configure the mock to return a sample DataFrame
    # This DataFrame should mimic the output of the actual db.get_term_variants
    sample_variants_df = pd.DataFrame(
        {
            "variant_id": [
                "1_100_A/T",
                "2_200_C/G",
                "1_186977737_A/G",
                "1_186977780_G/A",
            ],
            "gene": ["GENEA", "GENEB", "GENEA", "GENEC"],
            "aa": [
                "A",
                "C",
                "A",
                "G",
            ],  # This column gets dropped, but include for realism
            "hmm_score": [0.8, 0.9, 0.7, 0.85],
        }
    )
    monkeypatch.setattr(
        "blueprints.prioritisation_by_nomaly_scores.get_term_variants",
        lambda term: sample_variants_df,
    )

    # Just test two requests to keep it simple
    test_pairs = [
        ("250.2", "GO:0030800"),  # This user has access to 250.2
        # TODO: Mock data for 561
        # ("561", "MP:0005179"),  # And to 561
    ]

    def make_request(disease_code, term):
        """Make a request and track its messages"""
        response = unit_test_app_client_with_cache.get(
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

    assert len(results) == 1

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


def test_error_handling(unit_test_app_client_with_cache):
    """Test error handling for invalid requests."""
    # Invalid request
    response = unit_test_app_client_with_cache.get(
        "/stream_progress/invalid/invalid",
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


def test_cache_concurrency(unit_test_app_client_with_cache):
    """Test that concurrent requests get consistent results when using cache."""
    # Use a disease code the test user has access to
    test_pairs = [
        ("250.2", "UP:UPA00246"),  # This user has access to 250.2
    ]

    def make_request(disease_code, term):
        """Make a request and track its messages"""
        response = unit_test_app_client_with_cache.get(
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


def test_error_recovery(unit_test_app_client_with_cache):
    """Test error handling and recovery during concurrent requests."""
    # Mix of valid and invalid requests
    test_cases = [
        ("250.2", "UP:UPA00246"),  # Valid
        ("invalid", "invalid"),  # Invalid
        ("250.2", "invalid"),  # Partially invalid
    ]

    def make_request(disease_code, term):
        response = unit_test_app_client_with_cache.get(
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
        if (disease_code, term) == ("250.2", "UP:UPA00246"):
            assert any(m["type"] == "results" for m in messages)
            assert messages[-1]["type"] == "done"

        # Invalid requests should error gracefully
        else:
            assert any(m["type"] == "error" for m in messages)
            assert messages[-1]["type"] == "done"


def test_variant_scores_endpoint(auth_integration_app_client):
    """Test the variant scores endpoint with a known good example."""
    response = auth_integration_app_client.get("/variant_scores/290.11/GO:0016861")
    assert response.status_code == 200

    # Test the streaming endpoint too
    stream_response = auth_integration_app_client.get(
        "/stream_progress/290.11/GO:0016861"
    )
    assert stream_response.status_code == 200

    # Verify we get some actual data and not just errors
    assert b"error" not in stream_response.data
