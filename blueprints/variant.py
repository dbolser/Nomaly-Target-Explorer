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

logger = logging.getLogger(__name__)

# Create the blueprint
variant_bp = Blueprint("variant", __name__, template_folder="../templates")

# Dictionary to store background task results
phewas_results = {}


# Route to render the Variant page
@variant_bp.route("/variant/<string:variant>", methods=["POST", "GET"])
def show_variant(variant):
    try:
        # Get the ancestry parameter from request
        ancestry = session.get("ancestry", default="EUR")

        services: ServiceRegistry = current_app.extensions["nomaly_services"]
        genotype_service: GenotypeService = services.genotype
        nomaly_data_service: NomalyDataService = services.nomaly_data

        # Validate the variant format
        format_str = r"(\d+|MT|[XY])_\d+_[ACGT]+_[ACGT]+"
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

        # Check plink_id exists in the genotype service
        if plnkified_variant_id not in genotype_service.plink_variant_ids:
            logger.info(
                f"Variant {plnkified_variant_id} not found in genotype service, trying reverse alleles"
            )
            plnkified_variant_id = f"{chrom}:{pos}_{allele2}/{allele1}"

            if plnkified_variant_id not in genotype_service.plink_variant_ids:
                return render_template(
                    "error.html",
                    error=f"Variant at {chrom}:{pos} not found genotype data!",
                )

        # Use the 'nomaly_data service' to lookup variant information.
        variant_info = nomaly_data_service.get_variant_info_nomaly(nomalized_variant_id)

        if not variant_info:
            return render_template(
                "error.html",
                error=f"Variant {nomalized_variant_id} not found in the Nomaly Variant data!",
            )

        gene = variant_info["gene_id"]
        rsid = variant_info["RSID"]
        genotypeing_allele1 = plnkified_variant_id.split("_")[1][0]
        genotypeing_allele2 = plnkified_variant_id.split("_")[1][2]

        # Prepare data for the template
        variant_data = {
            "nomalized_variant_id": nomalized_variant_id,
            # We pass this through to the template so we can call phewas correctly
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


# Background task function
def background_task(
    variant: str,
    genotype_service: GenotypeService,
    phenotype_service: PhenotypeService,
    ancestry: str = "EUR",
    flush: bool = False,
):
    try:
        phewas_df = run_phewas_or_load_from_cache(
            variant, genotype_service, phenotype_service, ancestry, no_cache=flush
        )

        phewas_df = phewas_df.sort_values(by="p_value", ascending=True)

        num_sig_1_star = (phewas_df["p_value"] < 0.05).sum()
        num_sig_2_star = (phewas_df["p_value"] < 0.01).sum()
        num_sig_3_star = (phewas_df["p_value"] < 0.001).sum()

        # Format the results
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


# Endpoint to trigger the PheWAS task
@variant_bp.route("/run-phewas/<string:variant>", methods=["POST"])
def run_phewas(variant):

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


# Endpoint to get the PheWAS results
@variant_bp.route("/phewas-result/<string:variant>", methods=["GET"])
def get_phewas_result(variant):
    if variant not in phewas_results:
        return jsonify({"result": "Processing...", "associations": []})

    return jsonify(phewas_results[variant])
