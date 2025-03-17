import numpy as np
import pandas as pd
import plotly.express as px

from config import Config


# Yay globals!

RESOURCE_DATA_DIR = Config.RESOURCE_DATA_DIR
pharos_path = RESOURCE_DATA_DIR / "pharos_api_query.out"
pp_path = RESOURCE_DATA_DIR / "pp_by_gene.tsv"

pharos = pd.read_csv(pharos_path, sep="\t", encoding="ISO-8859-1").rename(
    columns={"#symbol": "gene"}
)

pp = pd.read_csv(pp_path, sep="\t", index_col=0).rename(
    columns={"Description": "drug_program_indication"}
)
pp["gene"] = pp.index


# ------------------------------------------------------------------------------#
# plotting functions
# ------------------------------------------------------------------------------#


def qqstats(dfstats):
    tags = dfstats.columns[dfstats.columns.str.endswith("_pvalue")]

    melt_stats = pd.melt(
        dfstats, id_vars=["term"], value_vars=tags, var_name="tag", value_name="P_obs"
    )
    melt_stats.dropna(inplace=True)
    # add metric column as tag
    melt_stats["test"] = melt_stats["tag"].apply(lambda x: x.split("_pvalue")[0])
    # # drop the tag column
    # melt_stats.drop(columns='tag', inplace=True)
    # add -log10(P_obs) column
    try:
        melt_stats["-log10(observed)"] = -np.log10(melt_stats["P_obs"])
    except Exception as e:
        print(f'Exception "{e}" encourterd for {melt_stats["term"][0]}')
    # sort the table by P_obs
    melt_stats.sort_values("P_obs", inplace=True)
    # add -log10(expected) column
    for tag in tags:
        len_tag = len(melt_stats[melt_stats["tag"] == tag])
        melt_stats.loc[melt_stats["tag"] == tag, "-log10(expected)"] = -np.log10(
            np.linspace(0 + 1 / len_tag, 1 - 1 / len_tag, len_tag)
        )
    # Merge description column if available
    if "description" in dfstats.columns:
        mapping = dfstats.set_index("term")["description"].to_dict()
        melt_stats["description"] = melt_stats["term"].map(mapping)
    return melt_stats


def make_qqplot(plot_df):
    try:
        melt_stats = qqstats(plot_df)
    except Exception as e:
        print(f'Exception "{e}" encountered')
        # Create a simple empty figure instead of returning None
        fig = px.scatter(x=[0], y=[0])
        fig.update_layout(
            title="Error generating QQ plot",
            annotations=[
                dict(
                    text=f"Error: {str(e)}",
                    showarrow=False,
                    xref="paper",
                    yref="paper",
                    x=0.5,
                    y=0.5,
                )
            ],
        )
        return fig

    # Define fixed colors for each test statistic
    color_map = {
        "mwu": "#636EFA",  # blue
        "mcc": "#EF553B",  # red
        "yjs": "#00CC96",  # green
        "lrp": "#AB63FA",  # purple
        "metric1": "#FFA15A",  # orange
        "lrn_protective": "#19D3F3",  # light blue
    }
    # Define legend order based on keys order in color_map
    legend_order = ["mwu", "mcc", "yjs", "lrp", "metric1", "lrn_protective"]

    # Add a scatter plot with plotly using consistent colors and legend order
    xlabel = "-log10(expected)"
    ylabel = "-log10(observed)"
    fig = px.scatter(
        melt_stats,
        x=xlabel,
        y=ylabel,
        color="test",
        color_discrete_map=color_map,
        category_orders={"test": legend_order},
        hover_name="term",
        hover_data={
            "description": True if "description" in melt_stats.columns else False
        },
        # title=f'{disease_select} QQ plot'
    )
    # add the diagonal line
    lims = [0, melt_stats["-log10(expected)"].max()]
    fig.add_scatter(
        x=lims,
        y=lims,
        mode="lines",
        name="Expected",
        line=dict(color="gray", dash="dash"),
    )

    # figure size and make it responsive
    fig.update_layout(
        width=600,
        height=400,
        autosize=True,
        margin=dict(l=50, r=50, t=30, b=50),
    )

    return fig
