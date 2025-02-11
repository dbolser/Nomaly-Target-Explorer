from flask import Blueprint, render_template, jsonify, request
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
        # Parse variant information regardless of format
        if ":" in variant:
            # Format: "5:33951588_C/G"
            chrom = variant.split(":")[0].replace("chr", "")
            rest = variant.split(":")[1]
            pos = rest.split("_")[0]
            alleles = rest.split("_")[1].replace("/", "_")
        else:
            # Format: "5_33951588_C_G"
            parts = variant.split("_")
            chrom = parts[0].replace("chr", "")
            pos = parts[1]
            alleles = f"{parts[2]}_{parts[3]}"  # Reconstruct alleles

        # Normalize the variant format
        normalized_variant = f"{chrom}_{pos}_{alleles}"

        # Validate the normalized format
        format_str = r"(chr)?(\d+|MT|[XY])_\d+_[ACGT]+_[ACGT]+"
        if not re.match(format_str, normalized_variant):
            return render_template(
                "error.html", error="Invalid variant ID format: " + variant
            )

        # Prepare data for the template
        variant_data = {
            "variant_id": f"{chrom}:{pos}_{alleles.replace('_', '/')}",  # Display format
            "chromosome": chrom,
            "position": pos,
            "allele1": alleles.split("_")[0],
            "allele2": alleles.split("_")[1],
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
            variant, None, services
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
