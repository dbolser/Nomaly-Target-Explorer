import json
import logging

import numpy as np
import pandas as pd
from flask import (
    Blueprint,
    current_app,
    jsonify,
    render_template,
    request,
    session,
    url_for,
    redirect,
)

from blueprints.gwas import format_gwas_results, run_gwas
from blueprints.nomaly import make_qqplot
from data_services import (
    PhenotypeService,
    ServiceRegistry,
    StatsRegistry,
    StatsService,
)
from db import (
    get_all_phecodes,
    get_phecode_info,
    get_term_domains,
    get_term_genes,
    get_term_names,
)

logger = logging.getLogger(__name__)
phecode_bp = Blueprint("phecode", __name__, template_folder="../templates")


@phecode_bp.route("/random_phecode", methods=["GET"])
def get_random_phecode():
    """Get a random phecode."""
    phecodes = get_all_phecodes()
    return phecodes.sample(1).iloc[0]["phecode"]


@phecode_bp.route("/phecode/<string:phecode>", methods=["GET"])
def show_phecode(phecode):
    """Show the phecode page."""

    # Get Run version and ancestry from session
    run_version = session.get("run_version", "Run-v1")
    ancestry = session.get("ancestry", "EUR")

    # Store values in session to ensure they persist
    session["run_version"] = run_version
    session["ancestry"] = ancestry

    try:
        # Get services while we're in 'app context'
        services: ServiceRegistry = current_app.extensions["nomaly_services"]

        # NOTE: We inject the appropriate services into 'backend' functions
        # (dependency injection)
        phecode_data = get_phecode_data(phecode, services.phenotype, ancestry)

        # Add runbatch and ancestry to the data for display AFTER getting phecode data
        phecode_data["runbatch"] = run_version
        phecode_data["ancestry"] = ancestry

        # Cobble together a half functional systsem..
        phecode_data["show_gwas"] = request.args.get("gwas") == "1"

        # And finally shove the whole thing into the template
        return render_template("phecode.html", data=phecode_data)
    except Exception as e:
        logger.error(f"Failed to get Phecode data for {phecode}: {e}")
        return jsonify({"error": "Failed to get Phecode data"}), 500


def get_phecode_data(
    phecode,
    phenotype_service: PhenotypeService,
    ancestry: str,
) -> dict:
    """Get the data for a phecode."""
    data = get_phecode_info(phecode)

    # Get case counts for the phecode
    case_counts = phenotype_service.get_case_counts_for_phecode(phecode, ancestry)

    data["population"] = ancestry or "EUR"
    data["affected"] = case_counts["affected"]
    data["excluded"] = case_counts["excluded"]
    data["control"] = case_counts["control"]

    return data


# Called by phecode.html
@phecode_bp.route("/nomaly-stats/<string:phecode>", methods=["POST"])
def get_nomaly_stats(phecode, run_version=None, ancestry=None):
    """Get nomaly stats"""

    if run_version is None or ancestry is None:
        # Get Run version and ancestry from the session
        run_version = session.get("run_version", "Run-v1")
        ancestry = session.get("ancestry", "EUR")

    try:
        services: ServiceRegistry = current_app.extensions["nomaly_services"]

        stats_registry: StatsRegistry = services.stats_registry
        stats_handler = stats_registry.get(run_version, ancestry)

        phecode_stats = stats_handler.get_phecode_stats(phecode)

        logger.info(f"Got {len(phecode_stats)} stats for {phecode}")
    except Exception as e:
        logger.error(f"Failed to get Nomaly stats for {phecode}: {e}")
        return jsonify({"error": "Failed to get Nomaly stats"}), 500

    return prepare_nomaly_stats_response(phecode, phecode_stats)


def prepare_nomaly_stats_response(phecode, phecode_stats):
    """Prepare the JSON response for nomaly stats."""

    # Rename columns (for some reason)
    renames = {
        x: x.replace("roc_stats_", "")
        for x in phecode_stats.columns
        if x.startswith("roc_stats_")
    }
    phecode_stats.rename(columns=renames, inplace=True)

    # Pull out pvalue_columns into a new dataframe
    pvalue_columns = [x for x in phecode_stats.columns if x.endswith("_pvalue")]
    plot_df = phecode_stats.loc[:, pvalue_columns]

    # NOTE: There are two pvalues we don't like...
    plot_df.drop(columns=["lrp_protective_pvalue", "tti_pvalue"], inplace=True)

    # Set the metric1_pvalue to na if it is 1
    # plot_df[plot_df["metric1_pvalue"].isna(), "metric1_pvalue"] = 1

    # First get the term 'names' for each term (term_id)

    plot_df["term"] = plot_df.index

    # NOTE: term = term_id = short code, term_name = term description)
    term_list = list(plot_df["term"])
    term_name_dict = get_term_names(term_list)

    # Add description column before plotting
    plot_df = plot_df.assign(
        description=plot_df["term"].map(lambda x: term_name_dict.get(x, "-"))
    )

    # Instead of HTML, send the plot data
    plot_data = get_qqplot_data(plot_df)

    # Format data for table display... why?
    plot_df = show_datatable_nomaly_stats(plot_df, phecode)

    # Add a link to the term page for each term
    # TODO: This should probably be done inside the template!
    plot_df["term"] = plot_df["term"].map(
        lambda x: f'<a href="{url_for("phecode_term.show_phecode_term", phecode=phecode, term=x)}" target="_blank">{x}</a>'
    )

    pval_nondirect = ["mwu_pvalue", "mcc_pvalue", "yjs_pvalue", "lrp_pvalue"]
    pval_pos = ["metric1_pvalue"]
    pval_neg = ["lrn_protective_pvalue"]
    columns_pval = pval_nondirect + pval_pos + pval_neg

    column_display_names = get_column_display_names()
    base_columns = ["minrank", "term", "name", "domain"]

    response = {
        "plotData": plot_data,  # Send JSON data instead of HTML
        "affected": phecode_stats["num_rp"].values[0],
        "control": phecode_stats["num_rn"].values[0],
        "data": plot_df.replace("nan", "1.00e+00").to_dict(orient="records"),
        "columns": base_columns + columns_pval,
        "columnNames": [
            column_display_names[col]["display"] for col in base_columns + columns_pval
        ],
        "columnTooltips": [
            column_display_names[col]["tooltip"] for col in base_columns + columns_pval
        ],
        "defaultColumns": base_columns + ["mwu_pvalue", "metric1_pvalue", "mcc_pvalue"],
        "numColumns": columns_pval,
    }

    return jsonify(response)


def get_qqplot_data(plot_df_pval):
    """Return the data needed to create the plot in JavaScript"""
    fig = make_qqplot(plot_df_pval)

    # Check if figure was created successfully
    if fig is None:
        # Return an empty plot to avoid errors
        return json.dumps({"data": [], "layout": {"title": "No data available"}})

    # Convert numpy arrays to Python native types
    def convert_numpy(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, dict):
            return {k: convert_numpy(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy(i) for i in obj]
        else:
            return obj

    # Convert the figure to dict and then clean up numpy types
    fig_dict = fig.to_dict()
    fig_dict_clean = convert_numpy(fig_dict)

    # Now we can safely serialize to JSON
    return json.dumps(fig_dict_clean)


def show_datatable_nomaly_stats(plot_df, phecode, addgene=False):
    """Format and prepare nomaly stats for datatable display."""
    pval_nondirect = ["mwu_pvalue", "mcc_pvalue", "yjs_pvalue", "lrp_pvalue"]
    pval_pos = ["metric1_pvalue"]
    pval_neg = ["lrn_protective_pvalue"]

    # add minimum rank from any of the pvalues
    columns_rank = []
    for col in pval_nondirect + pval_pos:
        plot_df[col.replace("_pvalue", "_minrank")] = plot_df[col].rank(method="min")
        columns_rank.append(col.replace("_pvalue", "_minrank"))

    plot_df["minrank"] = plot_df[columns_rank].min(axis=1)
    plot_df.drop(columns=columns_rank, inplace=True)
    plot_df.sort_values("minrank", inplace=True)

    # limit to top 1000
    plot_df_filtered = plot_df.iloc[:1000]

    # But add the top 50 from lrn_protective_pvalue
    plot_df_filtered = pd.concat(
        [
            plot_df_filtered,
            plot_df.sort_values("lrn_protective_pvalue", ascending=True)[:50],
        ]
    )

    # Remove duplicate terms
    plot_df_filtered = plot_df_filtered.drop_duplicates(subset="term")

    # get term names and domains
    term_name_dict = get_term_names(plot_df_filtered["term"].tolist())
    term_domain_dict = get_term_domains(plot_df_filtered["term"].tolist())

    plot_df_filtered = plot_df_filtered.assign(
        name=plot_df_filtered["term"].map(lambda x: term_name_dict.get(x, "-"))
    )

    # Create the domain column...

    def map_term_to_domain(term):
        domains = term_domain_dict[term]
        if len(domains) < 10:
            return ", ".join(domains)
        else:
            return f"{len(domains)} domains"

    plot_df_filtered["domain"] = plot_df_filtered["term"].map(map_term_to_domain)

    if addgene:
        plot_df_filtered = add_gene_info_to_DataTable(plot_df_filtered, phecode)

    # Round all P-values to scientific notation
    pval_columns = pval_nondirect + pval_pos + pval_neg
    for col in pval_columns:
        plot_df_filtered[col] = plot_df_filtered[col].apply(
            lambda x: f"{float(x):0.2e}" if x != "None" else x
        )

    return plot_df_filtered.fillna("None")


def add_gene_info_to_DataTable(plot_df, phecode):
    term_gene_df = get_term_genes(plot_df["term"].tolist())

    # Get GWAS results for gene filtering
    gwas_data = run_gwas(phecode)
    sig_variants = format_gwas_results(gwas_data)
    genefilter = set(v["Gene"] for v in sig_variants)

    # Filter and format gene information
    term_gene_df_sig = term_gene_df[term_gene_df["gene"].isin(genefilter)].rename(
        columns={"gene": "sig gene"}
    )

    # Group genes by term
    term_gene_df = (
        term_gene_df.groupby("term")["gene"]
        .apply(lambda x: ", ".join(x) if len(x) < 5 else f"{len(x)} genes")
        .reset_index()
    )

    term_gene_df_sig = (
        term_gene_df_sig.groupby("term")["sig gene"]
        .apply(lambda x: ", ".join(x) if len(x) < 50 else f"{len(x)} genes")
        .reset_index()
    )

    # Merge gene information with main dataframe
    plot_df = plot_df.merge(term_gene_df, on="term", how="left")
    plot_df = plot_df.merge(term_gene_df_sig, on="term", how="left")

    return plot_df


def get_column_display_names():
    """Get the display names and tooltips for columns."""
    return {
        "minrank": {
            "display": "Min Rank",
            "tooltip": "Minimum rank across all statistical tests",
        },
        "term": {"display": "Term", "tooltip": "Term identifier"},
        "name": {"display": "Description", "tooltip": "Term description"},
        "domain": {"display": "Domain", "tooltip": "Domain categories"},
        "mwu_pvalue": {
            "display": "MWU P-value",
            "tooltip": "P-value for the Mann-Whitney U test",
        },
        "mcc_pvalue": {
            "display": "MCC P-value",
            "tooltip": "P-value for Matthews Correlation Coefficient",
        },
        "yjs_pvalue": {
            "display": "YJS P-value",
            "tooltip": "P-value for Youden J Statistic",
        },
        "lrp_pvalue": {
            "display": "LRP P-value",
            "tooltip": "P-value for Likelihood Ratio Positive",
        },
        "metric1_pvalue": {
            "display": "Metric1 P-value",
            "tooltip": "P-value for Metric 1",
        },
        "lrn_protective_pvalue": {
            "display": "LRN Protective P-value",
            "tooltip": "P-value for Likelihood Ratio Negative (Protective)",
        },
    }


# Called by phecode.html
@phecode_bp.route("/run-task/<string:phecode>", methods=["POST"])
def run_task(phecode):
    """Endpoint to run GWAS analysis."""
    try:
        results = run_gwas(phecode)
        formatted_results = format_gwas_results(results)
        return jsonify(
            {
                "status": "completed",
                "result": f"GWAS completed successfully with {len(formatted_results)} significant variants",
                "associations": formatted_results,
            }
        )
    except Exception as e:
        logger.exception(f"GWAS failed for {phecode}")
        error_message = f"GWAS analysis failed: {str(e)}"
        return jsonify({"status": "failed", "result": error_message}), 500


# Add a new route to update session settings
@phecode_bp.route("/update_settings/<string:phecode>", methods=["POST"])
def update_settings(phecode):
    """Update session settings and redirect back to phecode page."""
    # Get values from form
    run_version = request.form.get("run_version", "Run-v1")
    ancestry = request.form.get("ancestry", "EUR")

    print(f"Ancestry: {ancestry}")
    # Store in session
    session["run_version"] = run_version
    session["ancestry"] = ancestry

    # Redirect back to phecode page
    return redirect(url_for("phecode.show_phecode", phecode=phecode))


def main():
    """Main function for testing."""

    from config import Config

    services = ServiceRegistry.from_config(Config)

    stats_service: StatsService = services.stats_registry.get(
        run_version="Run-v1", ancestry="EUR"
    )

    # Some random tests of get_phecode_stats...
    print(stats_service.get_phecode_stats("332"))
    print(stats_service.get_phecode_stats("332", stats_types=["roc_stats_lrn_pvalue"]))
    print(
        stats_service.get_phecode_stats(
            "332", stats_types=["roc_stats_lrn_pvalue", "mwu_pvalue"]
        )
    )
    print(stats_service.get_phecode_stats("332", term="MP:0004957"))

    # Test of some of the weirdness...
    stats = stats_service.get_phecode_stats("332")
    prepare_nomaly_stats_response("332", stats)

    print("Done")

if __name__ == "__main__":
    main()
