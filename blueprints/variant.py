import logging
import re
import threading
import traceback

from flask import Blueprint, current_app, jsonify, render_template, request, session

from blueprints.phewas import format_phewas_results, run_phewas_or_load_from_cache
from data_services import (
    GenotypeService,
    NomalyDataService,
    PhenotypeService,
    ServiceRegistry,
)
from db import get_all_phecodes

logger = logging.getLogger(__name__)

# Create the blueprint
variant_bp = Blueprint("variant", __name__, template_folder="../templates")

# Dictionary to store background task results
phewas_results = {}


@variant_bp.route("/variant/<string:variant>", methods=["POST", "GET"])
def show_variant(variant):
    """Route to render the Variant page"""
    try:
        # Get the ancestry parameter from request
        ancestry = session.get("ancestry", default="EUR")

        services: ServiceRegistry = current_app.extensions["nomaly_services"]
        genotype_service: GenotypeService = services.genotype
        nomaly_data_service: NomalyDataService = services.nomaly_data

        # Validate the variant format
        format_str = r"(\d+|MT|[XY])_\d+_[ACGT]+_[ACGT]+$"
        if not re.match(format_str, variant):
            return render_template(
                "error.html", error="Invalid variant ID format: " + variant
            )

        # Parse variant information from 'URL format' (e.g. "5_33951588_C_G")
        parts = variant.split("_")
        chrom = parts[0].replace("chr", "")
        pos = parts[1]
        allele1 = parts[2]
        allele2 = parts[3]

        # Nom'alize the variant format
        nomalized_variant_id = f"{chrom}_{pos}_{allele1}/{allele2}"
        plnkified_variant_id = f"{chrom}:{pos}_{allele1}/{allele2}"

        # I decided not to use the mapping file...
        variant_info = nomaly_data_service.get_variant_info_nomaly(nomalized_variant_id)

        if not variant_info:
            nomalized_variant_id = f"{chrom}_{pos}_{allele2}/{allele1}"
            variant_info = nomaly_data_service.get_variant_info_nomaly(
                nomalized_variant_id
            )
            if not variant_info:
                logger.info(f"Variant {chrom}:{pos} not found in Nomaly Variant data!")

        if not variant_info:
            gene = "TBD"
            rsid = "TBD"
        else:
            gene = variant_info["gene_id"]
            rsid = variant_info["RSID"]

        # I decided not to use the mapping file...
        if plnkified_variant_id not in genotype_service.plink_variant_ids:
            plnkified_variant_id = f"{chrom}:{pos}_{allele2}/{allele1}"
            if plnkified_variant_id not in genotype_service.plink_variant_ids:
                logger.info(f"Variant {chrom}:{pos} not found in genotype data!")
                return render_template(
                    "error.html",
                    error=f"Variant {chrom}:{pos} not found in genotype data!",
                )

        genotypeing_allele1 = plnkified_variant_id.split("_")[1][0]
        genotypeing_allele2 = plnkified_variant_id.split("_")[1][2]

        # Prepare data for the template
        variant_data = {
            # We probably only need to pass one of these...
            "nomalized_variant_id": nomalized_variant_id,
            "plnkified_variant_id": plnkified_variant_id,
            "rsid": rsid,
            "chromosome": chrom,
            "position": pos,
            "gene": gene,
            "genotyping_allele1": genotypeing_allele1,
            "genotyping_allele2": genotypeing_allele2,
            "ancestry": ancestry,
        }

        return render_template("variant.html", data=variant_data)

    except Exception as e:
        error_msg = f"Error processing variant {variant}: {str(e)}"
        print(error_msg)
        print(traceback.format_exc())
        return render_template("error.html", error=error_msg)


@variant_bp.route("/run-phewas/<string:variant>", methods=["POST"])
def run_phewas(variant):
    """Endpoint to trigger the PheWAS task"""

    # Get flush parameter from request
    flush = request.args.get("flush", default="0") == "1"

    # Get the ancestry parameter from request
    ancestry = session.get("ancestry", default="EUR")

    # Clear any existing results for this variant
    # Clear any existing results for this variant
    if variant in phewas_results:
        del phewas_results[variant]

    # Get services from the app context before starting the thread
    services = current_app.extensions["nomaly_services"]
    genotype_service = services.genotype
    phenotype_service = services.phenotype

    # Start the background task using threading
    task_thread = threading.Thread(
        target=background_task,
        args=(variant, genotype_service, phenotype_service, ancestry, flush),
    )
    task_thread.daemon = True
    task_thread.start()

    return jsonify({"status": "Task started"}), 202


def background_task(
    variant: str,
    genotype_service: GenotypeService,
    phenotype_service: PhenotypeService,
    ancestry: str = "EUR",
    flush: bool = False,
):
    """Background task to run the PheWAS"""
    try:
        phewas_df = run_phewas_or_load_from_cache(
            variant, genotype_service, phenotype_service, ancestry, no_cache=flush
        )

        # Get the phecode data (from the database)
        phecode_data = get_all_phecodes()

        # Merge the phewas results with the phecode data
        phewas_df = phewas_df.merge(phecode_data, on="phecode", how="inner")

        # Sanity check
        assert len(phecode_data) == len(phewas_df)

        # Handle missing phecode descriptions
        # phewas_df["description"] = phewas_df["description"].fillna("Unknown")
        # phewas_df["sex"] = phewas_df["sex"].fillna("Unknown")
        # phewas_df["phecode_group"] = phewas_df["phecode_group"].fillna("Unknown")

        phewas_df = phewas_df.sort_values(by="p_value", ascending=True)

        num_sig_1_star = (phewas_df["p_value"] < 0.05).sum()
        num_sig_2_star = (phewas_df["p_value"] < 0.01).sum()
        num_sig_3_star = (phewas_df["p_value"] < 0.001).sum()

        # Format the results for display
        phewas_dict = format_phewas_results(phewas_df)

        summary_text = (
            f"PheWAS identified {num_sig_1_star}*, {num_sig_2_star}**, {num_sig_3_star}***"
            "phecodes with association p<0.05, p<0.01, p<0.001 respectively."
        )
        phewas_results[variant] = {
            "result": summary_text,
            "associations": phewas_dict,
        }

    except Exception:
        error_msg = traceback.format_exc()
        logger.error(f"Error in background task: {error_msg}")
        phewas_results[variant] = {
            "result": f"Failed to get phecode-level stats for Variant {variant}, exception was <br> {error_msg}",
            "associations": [],
        }


@variant_bp.route("/phewas-result/<string:variant>", methods=["GET"])
def get_phewas_result(variant):
    """Endpoint to get the PheWAS results"""
    if variant not in phewas_results:
        return jsonify({"result": "Processing...", "associations": []})

    return jsonify(phewas_results[variant])
