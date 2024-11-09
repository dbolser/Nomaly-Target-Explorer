from flask import Blueprint, render_template
import pandas as pd

# Create the blueprint
page1_bp = Blueprint('page1', __name__)

# Sample data for each table
datadir = '/data/general/UKBB/Phenotypes/'

individual_phenotypes_data = pd.read_csv(datadir + 'individual_phenotypes.tsv', sep='\t')[[
    'Sex', 'suppop', 'icd10_any', 'phecodes'
]].rename(columns={'Sex':'sex', 'icd10_any': 'icd10', 'phecodes': 'phecode', 'suppop': 'population'}).to_dict(orient='records')

phecode_data = pd.read_csv(datadir + 'phecodes.tsv', sep='\t', dtype={'phecode': str})[[
    'phecode', 'description', 'sex', 'group', 'phecode_exclude_range', 'cases', 'exclude'
]].rename(columns={'exclude': 'excluded', 'phecode_exclude_range': 'phecode_exclude'}).to_dict(orient='records')

icd10_coding_data = pd.read_csv(datadir + "createUKBphenome/data/coding19.tsv", sep='\t')[[
    'coding', 'meaning'
]].to_dict(orient='records')

# Route for page1 with tables
@page1_bp.route('/page1')
def show_tables():
    return render_template('page1.html', 
                           individual_phenotypes=individual_phenotypes_data, 
                           phecode=phecode_data, 
                           icd10_coding=icd10_coding_data)