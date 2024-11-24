from flask import Blueprint, render_template, jsonify
import threading
import pandas as pd
import os
import traceback

from blueprints.gwas import variant_level_assoc, GWAS_PHENO_DIR
# from db import get_db_connection
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
    # ----------------------------------------------------- #
    # make sure variant id is in the correct format, and in our list
    # ----------------------------------------------------- #
    # variant = 'rs1234567'
    format_str = 'chr\\d+_\\d+_[ACGT]+_[ACGT]+'
    if not re.match(format_str, variant):
        return render_template('error.html', error='Invalid variant ID format: ' + variant)
    
    # # all_variants_in_nomaly = 
    # if variant not in all_variants_in_nomaly:
    #     return render_template('error.html', error='Variant ID not found in database: ' + variant)

    # ----------------------------------------------------- #
    # Run pheWAS on the variant
    # ----------------------------------------------------- #
    # GenotypeHDF5, ICD10HDF5 and phecodeHDF5 are needed