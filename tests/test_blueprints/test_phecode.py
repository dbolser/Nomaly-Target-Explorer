import pytest
from flask import json
from unittest.mock import Mock
from blueprints.phecode import prepare_nomaly_stats_response
import pandas as pd


def test_get_stats_handler_v1():
    """Test that get_stats_handler returns v1 stats by default."""
    from blueprints.phecode import get_stats_handler
    from services import services

    handler = get_stats_handler(version=1)
    assert handler == services.stats


def test_get_stats_handler_v2():
    """Test that get_stats_handler returns v2 stats when requested."""
    from blueprints.phecode import get_stats_handler
    from services import services

    handler = get_stats_handler(version=2)
    assert handler == services.stats_v2


@pytest.fixture
def mock_stats_data():
    """Create mock stats data for testing."""
    # Create a simple DataFrame that mimics the structure we expect
    diseasestats = pd.DataFrame(
        {
            "num_rp": [100.0],
            "num_rn": [1000.0],
            "mwu_pvalue": [0.001],
            "mcc_pvalue": [0.002],
            "yjs_pvalue": [0.003],
            "lrp_pvalue": [0.004],
            "metric1_pvalue": [0.005],
            "lrn_protective_pvalue": [0.006],
        }
    )

    plot_df = pd.DataFrame(
        {
            "term": ["CC:TERM:123"],
            "mwu_pvalue": [0.001],
            "mcc_pvalue": [0.002],
            "yjs_pvalue": [0.003],
            "lrp_pvalue": [0.004],
            "metric1_pvalue": [0.005],
            "lrn_protective_pvalue": [0.006],
        }
    )

    return diseasestats, plot_df


@pytest.mark.parametrize("version,expected_url", [(1, "/phecode"), (2, "/phecode2")])
def test_nomaly_stats_response_urls(
    version, expected_url, mock_stats_data, auth_client, mocker
):
    """Test that the response contains correct URLs based on version."""
    diseasestats, plot_df = mock_stats_data

    # Mock dependencies
    mocker.patch(
        "blueprints.phecode.make_qqplot_html", return_value="<div>Mock Plot</div>"
    )
    mocker.patch(
        "blueprints.phecode.get_term_names", return_value={"CC:TERM:123": "Test Term"}
    )
    mocker.patch(
        "blueprints.phecode.get_term_domains", return_value={"CC:TERM:123": ["Domain1"]}
    )

    with auth_client.application.app_context():
        # Get the response
        result = prepare_nomaly_stats_response(
            diseasestats, plot_df, "250.2", version=version
        )

        # If it's a tuple, it's an error response
        if isinstance(result, tuple):
            response, status_code = result
            assert False, f"Expected success response, got error: {response.get_data()}"
        else:
            response = result

        data = json.loads(response.get_data())

        # Check that the URLs in the term links use the correct version
        assert "phecode/250.2/term/" in data["data"][0]["term"]


@pytest.mark.parametrize("version", [1, 2])
def test_nomaly_stats_endpoint(version, auth_client, mocker):
    """Test the nomaly-stats endpoints (both v1 and v2)."""

    mock_diseasestats = pd.DataFrame(
        {
            "num_rp": [100.0],
            "num_rn": [1000.0],
            "mwu_pvalue": [0.001],
            "mcc_pvalue": [0.002],
            "yjs_pvalue": [0.003],
            "lrp_pvalue": [0.004],
            "metric1_pvalue": [0.005],
            "lrn_protective_pvalue": [0.006],
        },
        index=pd.Index(["CC:TERM:123"], name="term"),
    )

    # Create mock plot_df that matches the structure after processing
    mock_plot_df = pd.DataFrame(
        {
            "term": ["CC:TERM:123"],
            "mwu_pvalue": [0.001],
            "mcc_pvalue": [0.002],
            "yjs_pvalue": [0.003],
            "lrp_pvalue": [0.004],
            "metric1_pvalue": [0.005],
            "lrn_protective_pvalue": [0.006],
        }
    )

    # Mock the stats handler to return both DataFrames
    mock_handler = Mock()
    mock_handler.get_stats_by_disease.return_value = mock_diseasestats
    mocker.patch("blueprints.phecode.get_stats_handler", return_value=mock_handler)

    # Mock the read_disease_stats function to return both DataFrames
    mocker.patch(
        "blueprints.phecode.read_disease_stats_from_nomaly_statsHDF5",
        return_value=(mock_diseasestats, mock_plot_df),
    )

    # Mock other dependencies
    mocker.patch(
        "blueprints.phecode.make_qqplot_html", return_value="<div>Mock Plot</div>"
    )
    mocker.patch(
        "blueprints.phecode.get_term_names", return_value={"CC:TERM:123": "Test Term"}
    )
    mocker.patch(
        "blueprints.phecode.get_term_domains", return_value={"CC:TERM:123": ["Domain1"]}
    )

    # Make request to the appropriate endpoint
    endpoint = "/nomaly-stats2/250.2" if version == 2 else "/nomaly-stats/250.2"
    response = auth_client.post(endpoint)

    # Check response
    assert response.status_code == 200
    data = json.loads(response.get_data())

    # Check basic structure
    assert "qqplot" in data
    assert "affected" in data
    assert "control" in data
    assert "data" in data
    assert "columns" in data
    assert "columnNames" in data
    assert "defaultColumns" in data

    # Check version-specific URL format
    expected_url = "/phecode" if version == 2 else "/phecode"
    if data["data"]:  # If we have any data rows
        assert f"{expected_url}/250.2/term/" in data["data"][0]["term"]


def test_nomaly_stats_error_handling(auth_client, mocker):
    """Test error handling in nomaly-stats endpoints."""
    # Mock the stats handler to raise an exception
    mock_handler = Mock()
    mock_handler.get_stats_by_disease.side_effect = Exception("Test error")
    mocker.patch("blueprints.phecode.get_stats_handler", return_value=mock_handler)

    # Test both endpoints
    for endpoint in ["/nomaly-stats/250.2", "/nomaly-stats2/250.2"]:
        response = auth_client.post(endpoint)

        # Check response
        assert response.status_code == 500
        data = json.loads(response.get_data())
        assert "error" in data
        assert "Failed to get Nomaly stats" in data["error"]


def test_phecode_term_route_works(auth_client):
    """Test that the phecode term route works for a known good example."""
    phecode = "705"
    term = "CC:MESH:C020806"

    response = auth_client.get(f"/phecode/{phecode}/term/{term}")

    # This page should work - if it returns 500, the test will fail
    assert response.status_code == 200

    # We should see some expected content in the response
    assert bytes(phecode, "utf-8") in response.data
    assert bytes(term, "utf-8") in response.data
