import plotly.express as px
import plotly.io as pio
from flask import Blueprint, render_template, current_app

from blueprints.phecode import (
    read_disease_stats_from_nomaly_statsHDF5,
    get_stats_handler,
)
from blueprints.nomaly import make_qqplot

test_plot_bp = Blueprint("test_plot", __name__, template_folder="../templates")


@test_plot_bp.route("/test_plot/<string:phecode>", methods=["GET"])
def test_plot(phecode):
    """A minimal test page to display just the Plotly plot."""
    stats_handler = get_stats_handler(version=1)
    disease_stats, plot_df = read_disease_stats_from_nomaly_statsHDF5(
        stats_handler, phecode
    )

    plot_df["description"] = "test description"

    # Create the figure directly
    fig = make_qqplot(plot_df)

    # Generate the plot HTML
    plot_html = pio.to_html(
        fig,
        full_html=False,
        include_plotlyjs=True,  # Include plotly.js for standalone test
        config={"responsive": True},
    )

    return render_template("test_plot.html", plot_html=plot_html, phecode=phecode)
