from flask import Blueprint, render_template, request, jsonify, session
import threading
import time
import pandas as pd
import os
import traceback

from blueprints.gwas import variant_level_assoc, GWAS_PHENO_DIR
# from db import get_db_connection
from db import get_term_domains, get_term_names, get_term_genes, get_phecode_info

import plotly.express as px
import plotly.io as pio

from blueprints.nomaly import nomaly_scores, nomaly_stats, nomaly_genotype

# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #
# GWAS_PHENO_DIR = '/data/clu/ukbb/by_pheno/'

# ----------------------------------------------------- #
# Phecode Blueprint
# ----------------------------------------------------- #

# Create the blueprint
phecode_term_bp = Blueprint('phecode_term', __name__, template_folder='../templates')

# Route to render the Phecode page
@phecode_term_bp.route('/phecode/<string:phecode>/term/<string:term>', methods=['POST', 'GET'])
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

    stats_dict = nomaly_stats.get_stats_by_term_disease(term, phecode)
    # dict_keys(['num_rp', 'num_rn', 'mwu_pvalue', 'tti_pvalue', 'metric1_pvalue', 'roc_stats_mcc_value', 'roc_stats_mcc_pvalue', 'roc_stats_yjs_value', 'roc_stats_yjs_pvalue', 'roc_stats_lrp_value', 'roc_stats_lrp_pvalue', 'roc_stats_lrp_protective_value', 'roc_stats_lrp_protective_pvalue', 'roc_stats_lrn_value', 'roc_stats_lrn_pvalue', 'roc_stats_lrn_protective_value', 'roc_stats_lrn_protective_pvalue'])

    # get term names
    term_name = get_term_names([term])[term]
    term_domains = get_term_domains([term])[term]
    # print(term_domains)
    term_gene_df = get_term_genes([term])

    data['term'] = term
    data['termname'] = term_name
    data['domainlen'] = len(term_domains)
    # data['domains'] = ', '.join(term_domains)
    data['genelen'] = term_gene_df['gene'].nunique()
    data['genes'] = ', '.join(term_gene_df['gene'].unique())

    #   <span>Term: {{ data.term }}</span> |
    #             <span>Meaning: {{ data.termname }}</span> |
    #             <span> {{ data.domainlen }} Domains associated: {{ data.domains }}</span> |
    #             <span> {{data.genelen}} Genes associated: {{ data.genes }}</span> |

    # term scores
    

    # term variants matrix

    return render_template('phecode_term.html', phecode=phecode, term=term, data=data)

from blueprints.nomaly import pharos, pp
from db import get_term_domain_genes_variant, get_term_domain_genes

# Route to render the Phecode - term variants table
@phecode_term_bp.route('/phecode/<string:phecode>/term/<string:term>/tableVar', methods=['POST'])
def show_phecode_term_tableVar(phecode, term):

    # get term variants and other info including gene, domain
    term_df = get_term_domain_genes_variant(term)
    print(term_df.shape)

    # add pharos 
    term_df = term_df.merge(pharos, on='gene', how='left')
    print(term_df.shape)

    # add pp
    term_df = term_df.merge(pp, on='gene', how='left')
    print(term_df.shape)

    # ['variant_id', 'gene', 'aa', 'hmm', 'hmm_pos', 'sf', 'fa', 'tdl', 'tbio', 'pubmed_scores', 'generifs', 'antibodies', 'go_terms', 'tchem', 'tclin', 'drug_list', 'tdark', 'gwas', 'gwas_hits', 'disease_associations1', 'disease_associations2', 'classification', 'notes', 'indication_mesh_term', 'catall', 'succ_1_2', 'succ_2_3', 'succ_3_a', 'year_launch']

    print(term_df.columns)  
    
    # replace na
    term_df = term_df.fillna('None')

    nomalyResults = {
        'data': term_df.to_dict(orient='records'),
        'columns': ['variant_id', 'gene', 'tdl', 'tbio','classification', 'pubmed_scores', 'generifs', 'antibodies', 'go_terms', 'tchem', 'tclin', 'drug_list', 'tdark', 'drug_program_indication', 'gwas', 'gwas_hits', 'disease_associations1', 'disease_associations2', 'notes'],
        'defaultColumns': ['variant_id', 'gene', 'tdl', 'tbio','classification','gwas_hits', 'disease_associations1', 'drug_program_indication'],
        'numColumns': [],
    }

    return nomalyResults


@phecode_term_bp.route('/phecode/<string:phecode>/term/<string:term>/tableGene', methods=['POST'])
def show_phecode_term_tableGene(phecode, term):

    # get term variants and other info including gene, domain
    term_df = get_term_domain_genes(term)
    print(term_df.shape)

    # add pharos 
    term_df = pharos[pharos['gene'].isin(term_df['gene'])]
    print(term_df.shape)

    # add pp
    term_df = term_df.merge(pp, on='gene', how='left')
    print(term_df.shape)

    # ['variant_id', 'gene', 'aa', 'hmm', 'hmm_pos', 'sf', 'fa', 'tdl', 'tbio', 'pubmed_scores', 'generifs', 'antibodies', 'go_terms', 'tchem', 'tclin', 'drug_list', 'tdark', 'gwas', 'gwas_hits', 'disease_associations1', 'disease_associations2', 'classification', 'notes', 'indication_mesh_term', 'catall', 'succ_1_2', 'succ_2_3', 'succ_3_a', 'year_launch']

    print(term_df.columns)  
    
    # replace na
    term_df = term_df.fillna('None')

    nomalyResults = {
        'data': term_df.to_dict(orient='records'),
        'columns': ['gene', 'tdl', 'tbio','classification', 'pubmed_scores', 'generifs', 'antibodies', 'go_terms', 'tchem', 'tclin', 'drug_list', 'tdark', 'drug_program_indication', 'gwas', 'gwas_hits', 'disease_associations1', 'disease_associations2', 'notes'],
        'defaultColumns': ['gene', 'tdl', 'tbio','classification','gwas_hits', 'disease_associations1', 'drug_program_indication'],
        'numColumns': [],
    }

    return nomalyResults