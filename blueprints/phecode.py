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
phecode_bp = Blueprint('phecode', __name__, template_folder='../templates')

# Route to render the Phecode page
@phecode_bp.route('/phecode/<string:phecode>', methods=['GET'])
def show_phecode(phecode):

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

    return render_template('phecode.html', data=data)

# ----------------------------------------------------- #
# Nomaly Stats
# ----------------------------------------------------- #
@phecode_bp.route('/nomaly-stats/<string:phecode>', methods=['POST'])
def get_nomaly_stats(phecode):
    # Get the Nomaly stats for the Phecode
    # rows: term, columns: phecode, 3rd dim: statstype
    # num_rp num_rn mwu_pvalue tti_pvalue metric1_pvalue roc_stats_mcc_value roc_stats_mcc_pvalue roc_stats_yjs_value roc_stats_yjs_pvalue roc_stats_lrp_value roc_stats_lrp_pvalue roc_stats_lrp_protective_value roc_stats_lrp_protective_pvalue roc_stats_lrn_value roc_stats_lrn_pvalue roc_stats_lrn_protective_value roc_stats_lrn_protective_pvalue
    
    # ----------------------------------------------------- #
    # qq plot
    # ----------------------------------------------------- #
    try:
        diseasestats = nomaly_stats.get_stats_by_disease(phecode)
    except:
        return jsonify({"error": f"Failed to get Nomaly stats for Phecode {phecode}, ask admin to check logs."})
    # rename colums
    for col in diseasestats.columns:
        if '_p_value' in col:
            diseasestats = diseasestats.rename(columns={col: col.replace('_p_value', '_pvalue')})
            col = col.replace('_p_value', '_pvalue')
        if col.startswith('roc_stats_'):
            diseasestats = diseasestats.rename(columns={col: col.replace('roc_stats_', '')})

    # ['mwu_pvalue', 'mcc_pvalue', 'metric1_pvalue', 'yjs_pvalue', 'lrp_pvalue', 'lrn_protective_pvalue']

    # get columns with pvalues
    pval_nondirect = ['mwu_pvalue', 'mcc_pvalue', 'yjs_pvalue', 'lrp_pvalue']
    pval_pos = ['metric1_pvalue']
    pval_neg = ['lrn_protective_pvalue']
    columns_pval = pval_nondirect + pval_pos + pval_neg

    plot_df = diseasestats[columns_pval]
    plot_df['term'] = plot_df.index

    # set metric1_pvalue to na if it is 1
    plot_df['metric1_pvalue'] = plot_df['metric1_pvalue'].map(lambda x: pd.NA if x == 1 else x)

    # Generate an interactive plot using Plotly
    fig = make_qqplot(plot_df)

    # Convert the plot to HTML
    graph_html = pio.to_html(fig, full_html=False)

    # ----------------------------------------------------- #
    # Plot done, dataTable now
    # ----------------------------------------------------- #

    TERM_THRESHOLD = 0.000001
    # set metric1 threshold by total number of terms non-zero
    total_metric1 = sum(plot_df['metric1_pvalue']<1)
    TERM_THRESHOLD_metric1 = 1/total_metric1

    # speed up by requiring at least one pvalue< threshold or metric1 < threshold
    plot_df = plot_df[
        (plot_df[pval_nondirect + pval_pos + pval_neg] < TERM_THRESHOLD).any(axis=1) | (plot_df['metric1_pvalue'] < TERM_THRESHOLD_metric1)
    ]
    # plot_df = plot_df[(plot_df[pval_nondirect + pval_pos + pval_neg] < TERM_THRESHOLD).any(axis=1)]

    # add minimum rank from any of the pvalues
    columns_rank =[]
    for col in pval_nondirect + pval_pos:
        plot_df[col.replace('_pvalue', '_minrank')] = plot_df[col].rank(method='min')
        columns_rank.append(col.replace('_pvalue', '_minrank'))
        
    plot_df['minrank'] = plot_df[columns_rank].min(axis=1)
    plot_df.drop(columns = columns_rank, inplace=True)
    plot_df.sort_values('minrank', inplace=True)

    # limit to top 50
    plot_df = plot_df.iloc[:50]

    # get term names
    term_name_dict = get_term_names(plot_df['term'].tolist())

    plot_df['name'] = plot_df['term'].map(lambda x: term_name_dict.get(x, '-'))

    # get term to gene mapping
    print('getting term to gene mapping', flush=True)
    term_gene_df = get_term_genes(plot_df['term'].tolist())
    print('got term to gene mapping', flush=True)

    # filter gene by assoc var
    var_assoc_sig = read_gwas(phecode)
    genefilter = set([x['Gene'] for x in var_assoc_sig])
    term_gene_df_sig = term_gene_df[term_gene_df['gene'].isin(genefilter)].rename(columns={'gene': 'sig gene'})
    # term_gene_df_other = term_gene_df[~term_gene_df['gene'].isin(genefilter)]

    # group by term, no significance filter (to 
    term_gene_df = term_gene_df.groupby('term')['gene'].apply(
        lambda x: ', '.join(x) if len(x) < 5 else f"{len(x)} genes"
        ).reset_index()
    
    # group by term, use sig filter, uncomment above)
    term_gene_df_sig = term_gene_df_sig.groupby('term')['sig gene'].apply(
        lambda x: ', '.join(x) if len(x) < 50 else f"{len(x)} genes"
        ).reset_index()
    
    # # fill term_gene_df_sig NA with term_gene_df_other 
    # term_gene_df_other = term_gene_df.groupby('term')['gene'].apply(
    #     lambda x: ', '.join(x)
    #     ).reset_index().set_index('term')
    
    # term_gene_df_sig['gene'] = term_gene_df_sig['gene'].fillna(
    #     term_gene_df_sig['term'].map(lambda x: f"None ({term_gene_df_other.loc[x, 'gene']})")
    # )

    plot_df = plot_df.merge(term_gene_df, on='term', how='left')
    plot_df = plot_df.merge(term_gene_df_sig, on='term', how='left')

    # print(plot_df, flush=True)

    # get term to domain mapping
    term_domain_dict = get_term_domains(plot_df['term'].tolist())

    plot_df['domain'] = plot_df['term'].map(
        lambda x: ', '.join(term_domain_dict[x]) if len(term_domain_dict[x]) < 10 else [f"{len(term_domain_dict[x])} domains"]
    )

    # Replace NaN values with None (valid JSON format)
    plot_df = plot_df.fillna('None')

    # term to link
    plot_df['term'] = plot_df['term'].map(lambda x: f'<a href="/phecode/{phecode}/term/{x}">{x}</a>')

    nomalyResults = {
        'qqplot': graph_html,
        'data': plot_df.to_dict(orient='records'),
        'columns': ['minrank', 'term', 'name', 'gene', 'sig gene', 'domain'] + columns_pval,
        'defaultColumns': ['minrank','name', 'sig gene', 'domain'],
        #   + ['mwu_pvalue', 'metric1_pvalue', 'yjs_pvalue','mcc_pvalue'],
        'numColumns': columns_pval,
    }
    return nomalyResults

# ----------------------------------------------------- #
# Background task to run GWAS
# ----------------------------------------------------- #

# Dictionary to store background task results (in-memory for simplicity)
task_results = {}

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
    result = f"GWAS identified {assoc_sig.shape[0]} missense variants in {assoc_sig['gene_id'].nunique()} unique genes has association p<0.05."
    
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
        'CHR_BP_A1_A2': 'Variant',
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

    # Only show significant associations
    assoc_sig = assoc[assoc['P']<0.05]

    # change RSID to link
    assoc_sig['RSID'] = assoc_sig['RSID'].map(lambda x: f'<a href="https://www.ncbi.nlm.nih.gov/snp/{x}">v</a>')
    # assoc_sig['Variant'] = assoc_sig['Variant'].map(lambda x: f'<a href="https://www.ncbi.nlm.nih.gov/snp/?term={x.split('_')[0]}">{x}</a>')
    # assoc_sig['Variant'] = assoc_sig['Variant'] + assoc_sig['RSID']

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