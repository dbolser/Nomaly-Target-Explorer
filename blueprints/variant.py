from flask import Blueprint, render_template, jsonify
import threading
import pandas as pd
import os
import traceback

from blueprints.gwas import variant_level_assoc, GWAS_PHENO_DIR
from db import get_term_domains, get_term_names, get_term_genes, get_phecode_info
from blueprints.nomaly import nomaly_stats, make_qqplot
import plotly.express as px
import plotly.io as pio

import re

# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #
# GWAS_PHENO_DIR = '/data/clu/ukbb/by_pheno/'

# ----------------------------------------------------- #
# Phecode Blueprint
# ----------------------------------------------------- #

# Create the blueprint
variant_bp = Blueprint('variant', __name__, template_folder='../templates')

# Dictionary to store background task results
phewas_results = {}

# Background task function for PheWAS
def phewas_background_task(variant):
    try:
        # Placeholder for actual PheWAS computation
        result = f"PheWAS analysis completed for variant {variant}. Found 15 phenotype associations with P < 0.05"
        
        # More realistic dummy data
        associations = [
            {
                "Phecode": "250.2",
                "Description": "Type 2 diabetes",
                "Cases": 12453,
                "Controls": 85674,
                "OR": 1.45,
                "P": 2.3e-8
            },
            {
                "Phecode": "401.1",
                "Description": "Essential hypertension",
                "Cases": 28941,
                "Controls": 69186,
                "OR": 1.22,
                "P": 3.4e-6
            },
            {
                "Phecode": "272.1",
                "Description": "Hyperlipidemia",
                "Cases": 15678,
                "Controls": 82449,
                "OR": 1.18,
                "P": 4.7e-5
            },
            {
                "Phecode": "278.1",
                "Description": "Obesity",
                "Cases": 9876,
                "Controls": 88251,
                "OR": 1.31,
                "P": 8.9e-5
            },
            {
                "Phecode": "411.4",
                "Description": "Coronary atherosclerosis",
                "Cases": 7654,
                "Controls": 90473,
                "OR": 1.25,
                "P": 1.2e-4
            },
            {
                "Phecode": "362.2",
                "Description": "Retinopathy",
                "Cases": 3421,
                "Controls": 94706,
                "OR": 1.42,
                "P": 3.5e-4
            },
            {
                "Phecode": "585.3",
                "Description": "Chronic kidney disease",
                "Cases": 5632,
                "Controls": 92495,
                "OR": 1.29,
                "P": 4.8e-4
            },
            {
                "Phecode": "530.1",
                "Description": "Esophageal reflux",
                "Cases": 11234,
                "Controls": 86893,
                "OR": 0.82,
                "P": 8.3e-4
            },
            {
                "Phecode": "495.2",
                "Description": "Asthma",
                "Cases": 8765,
                "Controls": 89362,
                "OR": 0.88,
                "P": 1.7e-3
            },
            {
                "Phecode": "274.1",
                "Description": "Gout",
                "Cases": 2345,
                "Controls": 95782,
                "OR": 1.35,
                "P": 2.1e-3
            },
            {
                "Phecode": "696.4",
                "Description": "Psoriasis",
                "Cases": 1987,
                "Controls": 96140,
                "OR": 1.41,
                "P": 3.2e-3
            },
            {
                "Phecode": "244.1",
                "Description": "Hypothyroidism",
                "Cases": 6543,
                "Controls": 91584,
                "OR": 0.91,
                "P": 8.7e-3
            },
            {
                "Phecode": "714.1",
                "Description": "Rheumatoid arthritis",
                "Cases": 2198,
                "Controls": 95929,
                "OR": 1.28,
                "P": 1.2e-2
            },
            {
                "Phecode": "571.5",
                "Description": "Fatty liver disease",
                "Cases": 4321,
                "Controls": 93806,
                "OR": 1.19,
                "P": 2.8e-2
            },
            {
                "Phecode": "365.1",
                "Description": "Glaucoma",
                "Cases": 3456,
                "Controls": 94671,
                "OR": 0.89,
                "P": 4.1e-2
            }
        ]
        
        phewas_results[variant] = {
            "result": result,
            "associations": associations
        }
    except Exception as e:
        error_msg = f"Failed to run PheWAS for variant {variant}: {str(e)}"
        print(error_msg)
        print(traceback.format_exc())
        phewas_results[variant] = {"result": error_msg, "associations": []}

# Route to render the Variant page
@variant_bp.route('/variant/<string:variant>', methods=['POST', 'GET'])
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
        format_str = r"(chr)?\d+_\d+_[ACGT]+_[ACGT]+"
        if not re.match(format_str, normalized_variant):
            return render_template(
                "error.html", error="Invalid variant ID format: " + variant
            )

        # Prepare data for the template
        variant_data = {
            "variant_id": f"{chrom}:{pos}_{alleles.replace('_', '/')}",  # Display format
            "chromosome": chrom,
            "position": pos,
            "alleles": alleles.replace("_", "/"),  # Convert back to C/G for display
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

# Endpoint to trigger the PheWAS task
@variant_bp.route('/run-phewas/<string:variant>', methods=['POST'])
def run_phewas(variant):
    # Start the background task using threading
    task_thread = threading.Thread(target=phewas_background_task, args=(variant,))
    task_thread.start()
    return jsonify({"status": "Task started"}), 202

# Endpoint to get the PheWAS results
@variant_bp.route('/phewas-result/<string:variant>', methods=['GET'])
def get_phewas_result(variant):
    result = phewas_results.get(variant, {"result": "Processing...", "associations": []})
    return jsonify(result)
