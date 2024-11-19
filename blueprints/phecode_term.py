from flask import Blueprint, render_template, request, jsonify, session
import threading
import time
import pandas as pd
import os
import traceback

from blueprints.gwas import variant_level_assoc, GWAS_PHENO_DIR
# from db import get_db_connection
from db import get_term_domains, get_term_names, get_term_genes, get_phecode_info

from blueprints.nomaly import nomaly_stats, make_qqplot
import plotly.express as px
import plotly.io as pio


# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #
# GWAS_PHENO_DIR = '/data/clu/ukbb/by_pheno/'

# # load ontology_info.tsv
# ontology_info = pd.read_csv('/home/snow/datadir/1.3/GRCh38/ontology_info.tsv', sep='\t')
# # keep columns 0, 2
# ontology_info = ontology_info.iloc[:, [0, 2]]
# ontology_info.columns = ['term', 'name']

# ----------------------------------------------------- #
# Phecode Blueprint
# ----------------------------------------------------- #

# Create the blueprint
phecode_term_bp = Blueprint('phecode', __name__, template_folder='../templates')

# Route to render the Phecode page
@phecode_term_bp.route('/phecode/<string:phecode>/term/<string:term>', methods=['GET'])
def show_phecode_term(phecode, term):

    data = get_phecode_info(phecode)

    # # Example data (replace with actual data retrieval logic)
    # data = {
    #     "phecode": phecode,
    #     "sex": "Both",
    #     "affected": 150,
    #     "excluded": 10,
    #     "phecode_exclude": "None",
    #     "phecode_group": "Cardiovascular"
    # }

    # get term names
    term_name_dict = get_term_names([term])
    term_gene_df = get_term_genes([term])

    # term scores

    # term variants matrix


    return render_template('phecode_term.html', data=data)
