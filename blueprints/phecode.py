import json
import logging

import numpy as np
from flask import (
    Blueprint,
    current_app,
    jsonify,
    redirect,
    render_template,
    request,
    session,
    url_for,
)

from blueprints.gwas import format_gwas_results, run_gwas
from blueprints.nomaly import make_qqplot
from data_services import (
    NomalyDataService,
    PhenotypeService,
    ServiceRegistry,
    StatsRegistry,
    StatsService,
)
from db import (
    get_all_phecodes,
    get_phecode_info,
    get_term_genes,
    get_term_names,
)

logger = logging.getLogger(__name__)
phecode_bp = Blueprint("phecode", __name__, template_folder="../templates")


@phecode_bp.route("/phecode/<string:phecode>", methods=["GET"])
def show_phecode(phecode):
    """Show the phecode page."""

    # Get Run version and ancestry from session
    run_version = session.get("run_version", "Run-v1")
    ancestry = session.get("ancestry", "EUR")

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
        phecode_data["flush"] = request.args.get("flush") == "1"

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

    data["population"] = ancestry
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
        error_msg = f"""
            No statistics available for {phecode} using ancestry {ancestry} in
            {run_version}. Use the feedback form if you think this is a mistake.
        """
        return jsonify({"error": True, "message": error_msg}), 500

    phecode_stats = extract_pvalue_columns(phecode_stats)
    phecode_stats = add_minrank_column(phecode_stats)
    phecode_stats = add_term_names(phecode_stats)

    phecode_stats = phecode_stats.sort_values("minrank")

    return prepare_nomaly_stats_response(phecode, phecode_stats)


def extract_pvalue_columns(phecode_stats):
    """Extract p-value columns from the phecode stats dataframe."""

    # Rename columns (for some reason)
    renames = {
        x: x.replace("roc_stats_", "")
        for x in phecode_stats.columns
        if x.startswith("roc_stats_")
    }
    phecode_stats.rename(columns=renames, inplace=True)

    # Pull out pvalue_columns into a new dataframe
    pvalue_columns = [
        col
        for col in phecode_stats.columns
        if col.endswith("_pvalue")
        # NOTE: There some pvalues we don't like...
        and not col.startswith("tti_pvalue")
        and not col.startswith("lrn_pvalue")
        and not col.startswith("lrp_protective_pvalue")
    ]
    return phecode_stats.loc[:, pvalue_columns]


def add_term_names(phecode_stats):
    """Add the term names (descriptions) to the phecode stats dataframe."""

    # NOTE: term = term_id = short code, term_name = term description
    terms = phecode_stats.index
    phecode_stats["term"] = terms

    term_name_dict = get_term_names(list(terms))

    phecode_stats["name"] = phecode_stats["term"].map(
        lambda x: term_name_dict.get(x, "-")
    )

    return phecode_stats


def add_minrank_column(phecode_stats):
    """Add the minrank column to the phecode stats dataframe

    This is calculated as the minimum rank from any of the pvalues.

    In case of ties (such as when all values in a column are 1),
    we use a tie-breaking mechanism to ensure unique ranks.
    """
    # First get the basic ranks with method="min"
    primary_ranks = phecode_stats.rank(method="min")

    # For tie-breaking, use average method on the same data
    avg_ranks = phecode_stats.rank(method="average")

    # Now create a composite rank score - this adds a small fraction to break ties
    # The 0.001 multiplier ensures the secondary ranking doesn't override the primary
    composite_ranks = primary_ranks + (0.001 * avg_ranks)

    # Find the minimum rank across all columns for each row
    # Keep all original ranks as integers but ensure they remain unique
    min_ranks = composite_ranks.min(axis=1)

    phecode_stats["minrank"] = min_ranks.round(0)

    return phecode_stats


def prepare_nomaly_stats_response(phecode, phecode_stats):
    """Prepare the JSON response for nomaly stats."""

    # Set the metric1_pvalue to na if it is 1
    # plot_df[plot_df["metric1_pvalue"].isna(), "metric1_pvalue"] = 1

    # Instead of HTML, send the plot data
    plot_data = get_qqplot_data(phecode_stats)

    # Format data for table display... why?
    plot_df = show_datatable_nomaly_stats(phecode_stats, phecode)

    # Add a link to the term page for each term
    # TODO: This should probably be done inside the template!
    plot_df["term"] = plot_df["term"].map(
        lambda x: f'<a href="{url_for("phecode_term.show_phecode_term", phecode=phecode, term=x)}" target="_blank">{x}</a>'
    )

    column_display_names = get_column_display_names()

    columns = [
        "minrank",
        "term",
        "name",
        "mwu_pvalue",
        "metric1_pvalue",
        "mcc_pvalue",
        "yjs_pvalue",
        "lrp_pvalue",
        "lrn_protective_pvalue",
    ]

    response = {
        "plotData": plot_data,  # Send JSON data instead of HTML
        # "affected": phecode_stats["num_rp"].values[0],
        # "control": phecode_stats["num_rn"].values[0],
        "data": plot_df.replace("nan", "1.00e+00").to_dict(orient="records"),
        "columns": columns,
        "columnNames": [column_display_names[col]["display"] for col in columns],
        "columnTooltips": [column_display_names[col]["tooltip"] for col in columns],
        "defaultColumns": [
            "minrank",
            "term",
            "name",
            # "domain",
            "mwu_pvalue",
            "metric1_pvalue",
            "mcc_pvalue",
        ],
        #     "numColumns": [
        #         "minrank",
        #         "mwu_pvalue",
        #         "metric1_pvalue",
        #         "mcc_pvalue",
        #     ],
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

    # limit to top 1000
    plot_df_filtered = plot_df.iloc[:1000]

    # get term names and domains
    # term_name_dict = get_term_names(plot_df_filtered["term"].tolist())
    # term_domain_dict = get_term_domains(plot_df_filtered["term"].tolist())

    # plot_df_filtered = plot_df_filtered.assign(
    #     name=plot_df_filtered["term"].map(lambda x: term_name_dict.get(x, "-"))
    # )

    # Remove duplicate terms??
    plot_df_filtered = plot_df_filtered.drop_duplicates(subset="term")

    # Create the domain column...

    # def map_term_to_domain(term):
    #    domains = term_domain_dict[term]
    #    if len(domains) < 10:
    #        return ", ".join(domains)
    #    else:
    #        return f"{len(domains)} domains"

    # plot_df_filtered["domain"] = plot_df_filtered["term"].map(map_term_to_domain)

    if addgene:
        plot_df_filtered = add_gene_info_to_DataTable(plot_df_filtered, phecode)

    # Round all P-values to scientific notation
    pval_columns = [col for col in plot_df_filtered.columns if col.endswith("_pvalue")]
    for col in pval_columns:
        plot_df_filtered[col] = plot_df_filtered[col].apply(
            lambda x: f"{float(x):0.2e}" if x != "None" else x
        )

    return plot_df_filtered.fillna("None")


def add_gene_info_to_DataTable(plot_df, phecode):
    term_gene_df = get_term_genes(plot_df["term"].tolist())

    # Get GWAS results for gene filtering
    # gwas_data = run_gwas(phecode, "EUR", "x")
    # sig_variants = format_gwas_results(gwas_data)
    # genefilter = set(v["Gene"] for v in sig_variants)

    # Filter and format gene information
    # term_gene_df_sig = term_gene_df[term_gene_df["gene"].isin(genefilter)].rename(
    #    columns={"gene": "sig gene"}
    # )

    # Group genes by term
    term_gene_df = (
        term_gene_df.groupby("term")["gene"]
        .apply(lambda x: ", ".join(x) if len(x) < 5 else f"{len(x)} genes")
        .reset_index()
    )

    # term_gene_df_sig = (
    #    term_gene_df_sig.groupby("term")["sig gene"]
    #    .apply(lambda x: ", ".join(x) if len(x) < 50 else f"{len(x)} genes")
    #    .reset_index()
    # )

    # Merge gene information with main dataframe
    plot_df = plot_df.merge(term_gene_df, on="term", how="left")
    # plot_df = plot_df.merge(term_gene_df_sig, on="term", how="left")

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

    # Get flush from POST data
    flush = request.form.get("flush", "False") == 1

    # Get ancestry from session
    ancestry = session.get("ancestry", "EUR")

    services: ServiceRegistry = current_app.extensions["nomaly_services"]
    phenotype_service: PhenotypeService = services.phenotype
    nomaly_data_service: NomalyDataService = services.nomaly_data

    try:
        results = run_gwas(
            phecode, ancestry, phenotype_service, nomaly_data_service, no_cache=flush
        )
        logger.info(f"GWAS completed successfully with {len(results)} variants")

        sig_p = 1
        formatted_results = format_gwas_results(results, significance_threshold=sig_p)
        return jsonify(
            {
                "status": "completed",
                "result": f"GWAS completed successfully with {len(formatted_results)} 'significant' variants (p < {sig_p})",
                "associations": formatted_results,
            }
        )
    except Exception as e:
        logger.exception(f"GWAS failed for {phecode}")
        error_message = f"GWAS analysis failed: {str(e)}"
        return jsonify({"status": "failed", "result": error_message}), 500


# Called by phecode.html
@phecode_bp.route("/update_settings/<string:phecode>", methods=["POST"])
def update_settings(phecode):
    """Update session settings and redirect back to phecode page."""
    # Get values from form
    run_version = request.form.get("run_version", "Run-v1")
    ancestry = request.form.get("ancestry", "EUR")

    # Store in session
    session["run_version"] = run_version
    session["ancestry"] = ancestry

    # Redirect back to phecode page
    return redirect(url_for("phecode.show_phecode", phecode=phecode))


@phecode_bp.route("/random_phecode", methods=["GET"])
def get_random_phecode():
    """Get a random phecode."""
    phecodes = get_all_phecodes()
    random_phecode = phecodes.sample(1).iloc[0]["phecode"]
    return redirect(url_for("phecode.show_phecode", phecode=random_phecode))


def main():
    """Main function for testing."""

    # Need to test this route in app context I think...
    from app import create_app

    app = create_app("development")
    with app.app_context():
        with app.test_request_context():
            # Set session variables that would normally come from the user's session
            session["run_version"] = "Run-v1"
            session["ancestry"] = "SAS"

            get_nomaly_stats("332")
            get_nomaly_stats("333.3")

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
    # stats = stats_service.get_phecode_stats("332")
    # Needs to be called afer 'get_nomaly_stats'
    # prepare_nomaly_stats_response("332", stats)

    print("Done")

if __name__ == "__main__":
    main()
