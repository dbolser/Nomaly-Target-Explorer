from flask import Blueprint, render_template, request, jsonify, url_for

from db import get_phecode_info, get_term_domains, get_term_names, get_term_genes
from blueprints.gwas import run_gwas, format_gwas_results
from blueprints.nomaly import nomaly_stats, nomaly_stats_v2, make_qqplot

import pandas as pd
import plotly.io as pio

import logging

logger = logging.getLogger(__name__)
phecode_bp = Blueprint("phecode", __name__, template_folder="../templates")


def get_stats_handler(version=1):
    """Get the appropriate stats handler based on version."""
    return nomaly_stats_v2 if version == 2 else nomaly_stats


@phecode_bp.route("/phecode/<string:phecode>", methods=["GET"])
def show_phecode(phecode):
    data = get_phecode_info(phecode)
    data["runbatch"] = "Run v1"
    data["show_gwas"] = request.args.get("gwas") == "1"
    return render_template("phecode.html", data=data)


@phecode_bp.route("/phecode2/<string:phecode>", methods=["GET"])
def show_phecode2(phecode):
    data = get_phecode_info(phecode)
    data["runbatch"] = "Run v2"
    data["show_gwas"] = request.args.get("gwas") == "1"
    return render_template("phecode.html", data=data)


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


def read_disease_stats_from_nomaly_statsHDF5(stats_handler, phecode):
    """Read disease stats from HDF5 file using the provided stats handler."""
    try:
        diseasestats = stats_handler.get_stats_by_disease(phecode)
    except Exception as e:
        logger.error(
            f"Failed to get Nomaly stats for Phecode {phecode}, exception was {e}",
            exc_info=True,
        )
        raise

    # rename columns
    for col in diseasestats.columns:
        if col.startswith("roc_stats_"):
            diseasestats = diseasestats.rename(
                columns={col: col.replace("roc_stats_", "")}
            )

    # select columns with pvalues
    pval_nondirect = ["mwu_pvalue", "mcc_pvalue", "yjs_pvalue", "lrp_pvalue"]
    pval_pos = ["metric1_pvalue"]
    pval_neg = ["lrn_protective_pvalue"]
    columns_pval = pval_nondirect + pval_pos + pval_neg
    plot_df_pval = diseasestats[columns_pval].copy()
    plot_df_pval.loc[:, "term"] = plot_df_pval.index

    # set metric1_pvalue to na if it is 1
    plot_df_pval.loc[:, "metric1_pvalue"] = plot_df_pval["metric1_pvalue"].map(
        lambda x: None if x == 1 else x
    )

    return diseasestats, plot_df_pval


def make_qqplot_html(plot_df_pval):
    fig = make_qqplot(plot_df_pval)
    return pio.to_html(fig, full_html=False)


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

    # TODO: what is this doing?
    plot_df_filtered["domain"] = plot_df_filtered["term"].map(
        lambda x: ", ".join(term_domain_dict[x])
        if len(term_domain_dict[x]) < 10
        else [f"{len(term_domain_dict[x])} domains"]
    )

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


def prepare_nomaly_stats_response(diseasestats, plot_df, phecode, version=1):
    """Prepare the JSON response for nomaly stats."""
    if diseasestats is None or plot_df is None:
        return jsonify({"error": "Failed to get Nomaly stats"}), 500

    graph_html = make_qqplot_html(plot_df)
    plot_df = show_datatable_nomaly_stats(plot_df, phecode)

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
        "qqplot": graph_html,
        "affected": diseasestats["num_rp"].values[0],
        "control": diseasestats["num_rn"].values[0],
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


@phecode_bp.route("/nomaly-stats/<string:phecode>", methods=["POST"])
def get_nomaly_stats(phecode):
    """Get nomaly stats for v1."""
    try:
        stats_handler = get_stats_handler(version=1)
        diseasestats, plot_df = read_disease_stats_from_nomaly_statsHDF5(
            stats_handler, phecode
        )
        return prepare_nomaly_stats_response(diseasestats, plot_df, phecode, version=1)
    except Exception as e:
        logger.error(f"Failed to get Nomaly stats for {phecode}: {e}")
        return jsonify({"error": "Failed to get Nomaly stats"}), 500


@phecode_bp.route("/nomaly-stats2/<string:phecode>", methods=["POST"])
def get_nomaly_stats2(phecode):
    """Get nomaly stats for v2."""
    try:
        stats_handler = get_stats_handler(version=2)
        diseasestats, plot_df = read_disease_stats_from_nomaly_statsHDF5(
            stats_handler, phecode
        )
        return prepare_nomaly_stats_response(diseasestats, plot_df, phecode, version=2)
    except Exception as e:
        logger.error(f"Failed to get Nomaly stats for {phecode}: {e}")
        return jsonify({"error": "Failed to get Nomaly stats"}), 500
