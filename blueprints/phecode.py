from flask import Blueprint, render_template, request, jsonify, session
import threading
import time
import pandas as pd
import os
import traceback

from blueprints.gwas import variant_level_assoc, GWAS_PHENO_DIR
from db import get_db_connection

# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #
# GWAS_PHENO_DIR = '/data/clu/ukbb/by_pheno/'


# Create the blueprint
phecode_bp = Blueprint('phecode', __name__, template_folder='../templates')

# Dictionary to store background task results (in-memory for simplicity)
task_results = {}

# Route to render the Phecode page
@phecode_bp.route('/phecode/<string:phecode>', methods=['GET'])
def show_phecode(phecode):

    conn = get_db_connection()
    cur = conn.cursor()

    # filter the phecodes and the ICD10 table based on the query
    cur.execute(
        f"""
        SELECT p.phecode, p.description, p.sex, p.phecode_group, p.phecode_exclude, p.affected, p.excluded, icd10.icd10, icd10.meaning, icd10.icd10_count
        FROM phecode_definition p
        JOIN icd10_phecode ip ON p.phecode = ip.phecode
        JOIN icd10_coding icd10 ON ip.icd10 = icd10.icd10
        WHERE p.phecode = '{phecode}';
        """
    )
    results = cur.fetchall()

    cur.close()
    conn.close()

    # 4. Get the column names
    columns = [desc[0] for desc in cur.description]

    # 7. Convert the rows and columns into a Pandas DataFrame
    filtered_df = pd.DataFrame(results, columns=columns)

    # Add icd10_count to meaning
    filtered_df['meaning'] = filtered_df['meaning'] + ' | ' + filtered_df['icd10_count'].astype(str)

    # Group by 'description' and aggregate meanings into a list
    grouped = filtered_df.groupby(['description', 'phecode', 'sex', 'affected', 'excluded', 'phecode_exclude', 'phecode_group'
    ]).agg({
        'meaning': lambda x: list(x)  # Group meanings into a list
    }).reset_index()

    # Convert results to a list of dictionaries
    results = grouped.to_dict(orient='records')

    # results only has one row
    data = results[0]

    # # Example data (replace with actual data retrieval logic)
    # data = {
    #     "phecode": phecode,
    #     "sex": "Both",
    #     "affected": 150,
    #     "excluded": 10,
    #     "phecode_exclude": "None",
    #     "phecode_group": "Cardiovascular"
    # }
    return render_template('phecode.html', data=data)

# Background task function
def background_task(phecode):
    # ----------------------------------------------------- #
    # Task is to run GWAS and read the results
    # ----------------------------------------------------- #
    try:
        run_gwas_if_not_done(phecode)
    except Exception:
        task_results['result'] = f"Failed to get variant-level stats for Phecode {phecode}, exception was <br> {traceback.format_exc()}"
    

def run_gwas_if_not_done(phecode):
    # ----------------------------------------------------- #
    # GWAS results file path
    # ----------------------------------------------------- #
    output_prefix = f'phecode_{phecode}'
    gwas_path = f'{GWAS_PHENO_DIR}{output_prefix}.assoc_nomaly.tsv'

    # ----------------------------------------------------- #
    # Check if GWAS has been run for this Phecode
    # ----------------------------------------------------- #
    if not os.path.exists(gwas_path):
        task_results[phecode] = 'No GWAS data found for this Phecode. Processing...'
        assoc = variant_level_assoc(pheno_type = 'PheCode', code = phecode)
    else:
        assoc = pd.read_csv(gwas_path, sep='\t')
    
    # ----------------------------------------------------- #
    # If assoc is empty, return error message
    # ----------------------------------------------------- #
    if assoc is None:
        result = f"Failed to get variant-level stats for Phecode {phecode}, ask admin to check logs."

    # ----------------------------------------------------- #
    # If assoc is not empty, process the data
    # ----------------------------------------------------- #
    # remove any na in the columns
    assoc = assoc.dropna()

    assoc_sig = assoc[assoc['P']<0.05]
    result = f"GWAS identified {assoc_sig.shape[0]} variants in {assoc_sig['gene_id'].nunique()} unique genes has association p<0.05."
    
    # Store the result in the task_results dictionary
    task_results[phecode] = result    


def read_gwas(phecode):
    # ----------------------------------------------------- #
    # GWAS results file path
    # ----------------------------------------------------- #
    output_prefix = f'phecode_{phecode}'
    gwas_path = f'{GWAS_PHENO_DIR}{output_prefix}.assoc_nomaly.tsv'

    if not os.path.exists(gwas_path):
        return []
    
    assoc = pd.read_csv(gwas_path, sep='\t')

    # remove any na in the columns
    assoc = assoc.dropna()

    # assoc columns are:  Index(['nomaly_variant', 'gene_id', 'RSID', 'CHR_BP_A1_A2', 'F_A', 'F_U', 'OR', 'P'], dtype='object')
    # rename columns DataFrame
    assoc = assoc.rename(columns={
        'nomaly_variant': 'Variant',
        'gene_id': 'Gene',
        # 'RSID': 'RSID',
        # 'F_A': 'F_A',
        # 'F_U': 'F_U',
        # 'OR': 'OR',
        # 'P': 'P'
    })
    # columns to keep
    columns_to_keep = ['Variant', 'Gene', 'P', 'OR', 'F_A', 'F_U', 'RSID']
    assoc = assoc[columns_to_keep]

    assoc_sig = assoc[assoc['P']<0.05]

    return assoc_sig.to_dict(orient='records')


# Endpoint to trigger the background task
@phecode_bp.route('/run-task/<string:phecode>', methods=['POST'])
def run_task(phecode):
    # Start the background task using threading
    task_thread = threading.Thread(target=background_task, args=(phecode,))
    task_thread.start()
    return jsonify({"status": "Task started"}), 202

# Endpoint to get the result of the background task
@phecode_bp.route('/task-result/<string:phecode>', methods=['GET'])
def get_task_result(phecode):
    result = task_results.get(phecode, "Processing...")
    associations = read_gwas(phecode)
    return jsonify({"result": result,
                    "associations": associations})