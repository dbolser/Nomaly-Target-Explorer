from blueprints.nomaly import make_qqplot
import pandas as pd


def test_qqplot_figure_has_data():
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

    # Get the figure
    fig = make_qqplot(test_df)

    # Check the figure has data
    assert fig is not None, "Figure should not be None"
    assert fig.data is not None, "Figure data should not be None"

    assert len(fig.data) > 0, "Figure should have at least one trace"

    # Check for data points
    has_points = False
    for trace in fig.data:
        if hasattr(trace, "x") and hasattr(trace, "y") and len(trace.x) > 0:
            has_points = True
            break

    assert has_points, "Figure has no data points"

    # Check for legend
    assert hasattr(fig.layout, "showlegend")
    # assert fig.layout.showlegend is True, "Legend should be visible"
