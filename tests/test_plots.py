import json
from blueprints.phecode import get_qqplot_data
import pandas as pd


def test_qqplot_data_has_points():
    # Create test dataframe
    test_df = pd.DataFrame(
        {
            "term": ["A", "B", "C"],
            "mwu_pvalue": [0.01, 0.05, 0.001],
            "mcc_pvalue": [0.01, 0.05, 0.001],
            "description": ["Description A", "Description B", "Description C"],
            # Add other columns as needed
        }
    )
    # Generate the plot data
    plot_data_json = get_qqplot_data(test_df)
    plot_data = json.loads(plot_data_json)

    # Basic structure checks
    assert "data" in plot_data, "Plot JSON is missing 'data' key"
    assert len(plot_data["data"]) > 0, "Plot has no data traces"

    # Check for actual data points
    has_points = False
    for trace in plot_data["data"]:
        # Check if this is a data trace with points (not just a line)
        if "x" in trace and "y" in trace and len(trace["x"]) > 0:
            has_points = True
            break

    assert has_points, "Plot has no data points"

    # Check for legend
    assert "layout" in plot_data, "Plot JSON is missing 'layout' key"
    assert (
        "showlegend" not in plot_data["layout"] or plot_data["layout"]["showlegend"]
    ), "Legend is hidden"
