import json


def test_phecode_plot_data(client):
    # Assuming your API endpoint is something like /phecode/{id}/data
    response = client.get("/phecode/123.4/data")
    assert response.status_code == 200

    data = response.json
    assert "plotData" in data, "Response is missing plot data"

    plot_data = json.loads(data["plotData"])

    # Same assertions as before
    assert "data" in plot_data
    assert len(plot_data["data"]) > 0

    # Check for points
    point_count = 0
    for trace in plot_data["data"]:
        if "x" in trace and "y" in trace:
            point_count += len(trace["x"])

    assert point_count > 0, f"Plot should have data points, found {point_count}"
