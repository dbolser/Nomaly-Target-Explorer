from flask import Blueprint, render_template, jsonify, request
import threading
import pandas as pd
import os
import traceback

from blueprints.phewas import phecode_level_assoc, PHEWAS_PHENO_DIR

import re

# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #
# GWAS_PHENO_DIR = '/data/clu/ukbb/by_pheno/'

# ----------------------------------------------------- #
# Variant Blueprint
# ----------------------------------------------------- #

# Create the blueprint
variant_bp = Blueprint("variant", __name__, template_folder="../templates")


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

    # ----------------------------------------------------- #
    # Run pheWAS on the variant
    # ----------------------------------------------------- #
    # GenotypeHDF5, ICD10HDF5 and phecodeHDF5 are needed


# ----------------------------------------------------- #
# Background task to run PheWAS
# ----------------------------------------------------- #

# Dictionary to store background task results
phewas_results = {}


# Background task function
def background_task(variant, flush=False):
    try:
        result = run_phewas_if_not_done(variant, flush)
        if result is None:
            phewas_results[variant] = {
                "result": f"No associations found for variant {variant}",
                "associations": []
            }
        else:
            # result already contains the formatted data from run_phewas_if_not_done
            phewas_results[variant] = result
            
    except Exception:
        error_msg = traceback.format_exc()
        print(f"Error in background task: {error_msg}")  # Add logging
        phewas_results[variant] = {
            "result": f"Failed to get phecode-level stats for Variant {variant}, exception was <br> {error_msg}",
            "associations": []
        }


# Background task function for PheWAS
def run_phewas_if_not_done(variant, flush=False):
    """
    Check if PheWAS results exist for this variant.
    If not (or if flush=True), run PheWAS.
    Returns the processed DataFrame of results.
    """
    output_prefix = f"variant_{variant}"
    phewas_path = f"{PHEWAS_PHENO_DIR}{output_prefix}.assoc_nomaly.tsv"

    if not flush and os.path.exists(phewas_path):
        # Read existing results
        assoc = pd.read_csv(phewas_path, sep="\t")
    else:
        # Run PheWAS (which now handles file saving internally)
        variant_colon = variant.replace("_", ":")
        assoc = phecode_level_assoc(variant_colon)

    if assoc is None or assoc.empty:
        return None

    # Filter significant associations
    assoc_sig = assoc[assoc["p_value"] < 0.05]

    # Format the results
    result = f"PheWAS identified {assoc_sig.shape[0]} phecodes with association p<0.05."

    # Format associations for the table
    associations = []
    for _, row in assoc_sig.iterrows():
        associations.append(
            {
                "Phecode": row["phecode"],
                "Sex": row.get("sex", "TBD"),
                "Description": row.get("description", "TBD"),
                "Group": row.get("phecode_group", "TBD"),
                "Counts": f"{row['n_cases']}<br/>{row['n_controls']}",
                "RefAF": f"{row['ref_allele_freq_cases']:.5f}<br/>{row['ref_allele_freq_controls']:.5f}",
                "AltAF": f"{row['alt_allele_freq_cases']:.5f}<br/>{row['alt_allele_freq_controls']:.5f}",
                "Ref_HMOZ": f"{row.get('homozygous_ref_cases', 0)}<br/>{row.get('homozygous_ref_controls', 0)}",
                "Alt_HMOZ": f"{row.get('homozygous_alt_cases', 0)}<br/>{row.get('homozygous_alt_controls', 0)}",
                "HTRZ": f"{row.get('heterozygous_cases', 0)}<br/>{row.get('heterozygous_controls', 0)}",
                "P": f"{row['p_value']:.2e}<br/>",
                "OR": f"{row['odds_ratio']:.2f}<br/>",
            }
        )

    return {"result": result, "associations": associations}


# Endpoint to trigger the PheWAS task
@variant_bp.route("/run-phewas/<string:variant>", methods=["POST"])
def run_phewas(variant):
    # Get flush parameter from request
    flush = request.args.get("flush", default="0") == "1"
    
    # Clear any existing results for this variant
    if variant in phewas_results:
        del phewas_results[variant]
        
    # Start the background task using threading
    task_thread = threading.Thread(target=background_task, args=(variant, flush))
    task_thread.daemon = True  # Make thread daemon so it doesn't block shutdown
    task_thread.start()
    return jsonify({"status": "Task started"}), 202


# Endpoint to get the PheWAS results
@variant_bp.route("/phewas-result/<string:variant>", methods=["GET"])
def get_phewas_result(variant):
    if variant not in phewas_results:
        return jsonify({"result": "Processing...", "associations": []})
    return jsonify(phewas_results[variant])
