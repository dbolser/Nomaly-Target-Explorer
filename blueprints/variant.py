from flask import Blueprint, render_template, jsonify, request, current_app
import threading
import traceback
import re

from blueprints.phewas import get_formatted_phewas_data


# Create the blueprint
variant_bp = Blueprint("variant", __name__, template_folder="../templates")

# Dictionary to store background task results
phewas_results = {}


# Route to render the Variant page
@variant_bp.route("/variant/<string:variant>", methods=["POST", "GET"])
def show_variant(variant):
    try:
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
        nomalized_variant = f"{chrom}_{pos}_{allele1}/{allele2}"
        plnkified_variant = f"{chrom}:{pos}_{allele1}/{allele2}"

        # THE KEY THING TO NOTE HERE IS THAT PLINK MAY HAVE FLIPPED THE ALLELES,
        # AND WE'RE TOO LAZY TO FIX IT!

        # Now we use the 'nomaly_data service' to lookup comprehensive variant
        # information from the 'mapping' hack.
        app = current_app
        services = app.extensions["nomaly_services"]

        variant_info = (
            False
            or services.nomaly_data.get_variant_info_nomaly(nomalized_variant)
            or services.nomaly_data.get_variant_info_plinky(plnkified_variant)
        )

        if not variant_info:
            return render_template(
                "error.html", error="Variant not found in the Nomaly mapping data."
            )

        # Finally get the data we want... easy eh?
        gene = variant_info["gene_id"]
        rsid = variant_info["RSID"]
        genotypeing_allele1 = variant_info["CHR_BP_A1_A2"].split("_")[1][0]
        genotypeing_allele2 = variant_info["CHR_BP_A1_A2"].split("_")[1][2]

        # Prepare data for the template
        variant_data = {
            "variant_id": nomalized_variant,
            "rsid": rsid,
            "chromosome": chrom,
            "position": pos,
            "gene": gene,
            "genotyping_allele1": genotypeing_allele1,
            "genotyping_allele2": genotypeing_allele2,
        }

        return render_template("variant.html", data=variant_data)

    except Exception as e:
        error_msg = f"Error processing variant {variant}: {str(e)}"
        print(error_msg)
        print(traceback.format_exc())
        return render_template("error.html", error=error_msg)


# Background task function
def background_task(variant: str, services, flush: bool = False):
    try:

        # Get formatted PheWAS results using injected services
        results = get_formatted_phewas_data(
            variant, services, no_cache=flush
        )  # None to get all phecodes

        if not results:
            phewas_results[variant] = {
                "result": f"No associations found for variant {variant}",
                "associations": [],
            }
        else:
            # Filter significant associations using raw p-value
            assoc_sig = [r for r in results if r["p_value"] < 0.05]

            # Drop the p_value key from the dictionary
            for r in assoc_sig:
                del r["p_value"]

            result = (
                f"PheWAS identified {len(assoc_sig)} phecodes with association p<0.05."
            )
            phewas_results[variant] = {"result": result, "associations": assoc_sig}

    except Exception:
        error_msg = traceback.format_exc()
        print(f"Error in background task: {error_msg}")
        phewas_results[variant] = {
            "result": f"Failed to get phecode-level stats for Variant {variant}, exception was <br> {error_msg}",
            "associations": [],
        }


# Endpoint to trigger the PheWAS task
@variant_bp.route("/run-phewas/<string:variant>", methods=["POST"])
def run_phewas(variant):
    from flask import current_app

    # Get flush parameter from request
    flush = request.args.get("flush", default="0") == "1"

    # Clear any existing results for this variant
    if variant in phewas_results:
        del phewas_results[variant]

    # Get services from the app context before starting the thread
    services = current_app.extensions["nomaly_services"]

    # Start the background task using threading
    task_thread = threading.Thread(
        target=background_task, args=(variant, services, flush)
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
