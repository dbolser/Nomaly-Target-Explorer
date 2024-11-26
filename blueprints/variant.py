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
