from flask import Blueprint, render_template, jsonify

import pandas as pd
import os
import traceback

from blueprints.gwas import variant_level_assoc, GWAS_PHENO_DIR
from blueprints.phewas import phecode_level_assoc, PHEWAS_PHENO_DIR

# from db import get_db_connection
from db import get_term_domains, get_term_names, get_term_genes, get_phecode_info
from db import get_db_connection
from db import get_term_domain_genes_variant, get_term_domain_genes

from blueprints.nomaly import pharos, pp
from blueprints.nomaly import nomaly_stats

from errors import DataNotFoundError, GWASError
import logging


# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #
# GWAS_PHENO_DIR = '/data/clu/ukbb/by_pheno/'

# ----------------------------------------------------- #
# Phecode Blueprint
# ----------------------------------------------------- #

# Create the blueprint
phecode_term_bp = Blueprint("phecode_term", __name__, template_folder="../templates")

logger = logging.getLogger(__name__)


# Route to render the Phecode page
@phecode_term_bp.route(
    "/phecode/<string:phecode>/term/<string:term>", methods=["POST", "GET"]
)
def show_phecode_term(phecode, term):
    try:
        data = get_phecode_info(phecode)
        stats_dict = nomaly_stats.get_stats_by_term_disease(term, phecode)

        term_name = get_term_names([term])[term]
        term_domains = get_term_domains([term])[term]
        term_gene_df = get_term_genes([term])

        # Update data dictionary
        data.update(
            {
                "term": term,
                "termname": term_name,
                "domainlen": len(term_domains),
                "genelen": term_gene_df["gene"].nunique(),
                "genes": ", ".join(term_gene_df["gene"].unique()),
            }
        )

        return render_template(
            "phecode_term.html", phecode=phecode, term=term, data=data
        )

    except DataNotFoundError as e:
        logger.warning(f"Data not found: {str(e)}")
        return render_template("error.html", error=str(e)), 404
    except Exception as e:
        logger.error(f"Error in show_phecode_term: {str(e)}", exc_info=True)
        return render_template("error.html", error="An unexpected error occurred"), 500


# Route to render the Phecode - term variants table
@phecode_term_bp.route('/phecode/<string:phecode>/term/<string:term>/tableVar', methods=['POST'])
def show_phecode_term_tableVar(phecode, term):

    # get term variants and other info including gene, domain
    term_df = get_term_domain_genes_variant(term)
    print(f"Initial data shape: {term_df.shape}")
    print(f"Initial genes: {term_df['gene'].nunique()}")

    # add pharos
    print(f"Pharos genes: {len(pharos['gene'].unique())}")
    term_df = term_df.merge(pharos, on='gene', how='left')
    print(f"After pharos merge shape: {term_df.shape}")
    print(f"Genes after pharos merge: {term_df['gene'].nunique()}")

    # add pp
    print(f"PP genes: {len(pp['gene'].unique())}")
    term_df = term_df.merge(pp, on='gene', how='left')
    print(f"Final shape: {term_df.shape}")
    print(f"Final genes: {term_df['gene'].nunique()}")

    # Sample some data
    print("\nSample of final data:")
    print(term_df.head())

    # ['variant_id', 'gene', 'aa', 'hmm', 'hmm_pos', 'sf', 'fa', 'tdl', 'tbio', 'pubmed_scores', 'generifs', 'antibodies', 'go_terms', 'tchem', 'tclin', 'drug_list', 'tdark', 'gwas', 'gwas_hits', 'disease_associations1', 'disease_associations2', 'classification', 'notes', 'indication_mesh_term', 'catall', 'succ_1_2', 'succ_2_3', 'succ_3_a', 'year_launch']

    print(term_df.columns)  
    
    # replace na
    term_df = term_df.fillna('None')

    nomalyResults = {
        'data': term_df.to_dict(orient='records'),
        'columns': ['variant_id', 'gene', 'tdl', 'tbio','classification', 'drug_list', 'drug_program_indication', 'disease_associations1', 'disease_associations2', 'notes'],
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


def ensure_gwas(phecode):
    """
    Check if GWAS results for this phecode exist.
    If not, run GWAS.
    """
    output_prefix = f"phecode_{phecode}"
    gwas_path = f"{GWAS_PHENO_DIR}{output_prefix}.assoc_nomaly.tsv"
    if not os.path.exists(gwas_path):
        # Run GWAS
        assoc = variant_level_assoc(pheno_type="PheCode", code=phecode)
        if assoc is not None and not assoc.empty:
            # assoc is already saved by variant_level_assoc function call or we can save here
            pass
        else:
            print(f"GWAS failed or returned no results for {phecode}")


def ensure_phewas(variant_id):
    """
    Check if PheWAS results for this variant exist.
    If not, run PheWAS.
    """
    output_prefix = f"variant_{variant_id}"
    phewas_path = f"{PHEWAS_PHENO_DIR}{output_prefix}.assoc_nomaly.tsv"

    if os.path.exists(phewas_path):
        # If file exists, just return
        return

    try:
        assoc = phecode_level_assoc(variant_id)
        if assoc is None or assoc.empty:
            print(f"PheWAS returned no results for variant {variant_id}")
    except Exception as e:
        print(f"Error running PheWAS for variant {variant_id}: {e}")


def load_gwas_data(phecode, variant_id):
    """
    Load GWAS data for a given phecode and return row for this variant.
    Returns a dict with keys: GWAS_P, GWAS_OR, etc.
    If not found, returns None.
    """
    ensure_gwas(phecode)
    gwas_file = f"{GWAS_PHENO_DIR}phecode_{phecode}.assoc_nomaly.tsv"
    if not os.path.exists(gwas_file):
        return None
    gwas_df = pd.read_csv(gwas_file, sep="\t")

    # Debug: Print column names
    print(f"GWAS columns for {phecode}: {gwas_df.columns.tolist()}")

    # variant_id is like: 8_6870776_C/T
    # GWAS variant format: "CHR:BP_REF/ALT"
    # Convert variant_id to GWAS format
    # variant_id split: chrom_pos_ref/alt (already given)
    # Just replace first underscore with colon to get CHR:POS_REF/ALT
    # variant_gwas_format = variant_id.replace("_", ":", 1)

    variant_col = "nomaly_variant"

    print(variant_id)
    print(gwas_df[variant_col])

    row = gwas_df[gwas_df[variant_col] == variant_id]
    if row.empty:
        print(f"Variant {variant_id} not found in GWAS file")
        return None

    r = row.iloc[0]

    # Map expected column names to possible alternatives
    column_map = {
        "P": ["P", "p", "pvalue", "p_value", "P_value"],
        "OR": ["OR", "or", "odds_ratio"],
        "F_A": ["F_A", "freq_a", "freq_alt"],
        "F_U": ["F_U", "freq_u", "freq_ref"],
        "RSID": ["RSID", "rsid", "rs", "SNP"],
    }

    result = {}
    for key, alternatives in column_map.items():
        found_col = next((col for col in alternatives if col in gwas_df.columns), None)
        result[f"GWAS_{key}"] = r[found_col] if found_col else None

    return result


def load_phewas_data(variant_id, phecode):
    """
    Load PheWAS data for a given variant and filter by phecode.
    Returns dict with allele freq, OR, p-value, counts, etc.
    If not found, returns None.
    """
    ensure_phewas(variant_id)
    phewas_file = f"{PHEWAS_PHENO_DIR}variant_{variant_id}.assoc_nomaly.tsv"
    if not os.path.exists(phewas_file):
        return None

    phewas_df = pd.read_csv(phewas_file, sep="\t")
    phewas_df["phecode"] = phewas_df["phecode"].astype(str)
    phecode = str(phecode)
    phewas_row = phewas_df[phewas_df["phecode"] == phecode]
    if phewas_row.empty:
        return None
    r = phewas_row.iloc[0]
    return {
        "PheWAS_P": r["p_value"],
        "PheWAS_OR": r["odds_ratio"],
        "Ref_AF_Cases": r["ref_allele_freq_cases"],
        "Alt_AF_Cases": r["alt_allele_freq_cases"],
        "Ref_AF_Controls": r["ref_allele_freq_controls"],
        "Alt_AF_Controls": r["alt_allele_freq_controls"],
        "Cases_Count": r["n_cases"],
        "Controls_Count": r["n_controls"],
    }


def get_term_variants(term: str) -> pd.DataFrame:
    """
    Return a DataFrame of variant_id, gene, aa, hmm_score for the given term.
    """
    conn = get_db_connection()
    cur = conn.cursor()

    query = """
    SELECT term, variant_id, gene, aa, ABS(wild - mutant) AS hmm_score
    FROM variants_consequences
    INNER JOIN terms2snps USING (variant_id)
    WHERE term = %s
    """
    print(f"\nExecuting query for term: {term}")
    print(f"Query: {query}")

    cur.execute(query, (term,))
    results = cur.fetchall()
    print(f"Number of results from DB: {len(results)}")

    columns = [desc[0] for desc in cur.description]
    cur.close()
    conn.close()

    df = pd.DataFrame(results, columns=columns) if results else pd.DataFrame()
    print(f"Created DataFrame with shape: {df.shape}")

    # Convert numpy types to Python native types
    if not df.empty:
        df["hmm_score"] = df["hmm_score"].astype(float)
        print("Sample of data:")
        print(df.head())

    return df


@phecode_term_bp.route(
    "/phecode/<string:phecode>/term/<string:term>/tableVariantDetail",
    methods=["GET", "POST"],
)
def show_phecode_term_variant_detail(phecode, term):
    try:
        # Get variants from DB
        print(f"\nFetching variants for term: {term}")
        term_df = get_term_variants(term)
        print(f"Initial term_df shape: {term_df.shape}")

        if term_df.empty:
            print("No variants found for this term!")
            return jsonify(
                {"data": [], "columns": [], "defaultColumns": [], "numColumns": []}
            )

        data_records = []
        for idx, row in term_df.iterrows():
            variant_id = row["variant_id"]
            print(f"\nProcessing variant {idx+1}/{len(term_df)}: {variant_id}")
            standard_variant_id = variant_id.replace("/", "_")

            # Load PheWAS data
            phewas_data = load_phewas_data(standard_variant_id, phecode)
            print(f"PheWAS data found: {phewas_data is not None}")

            # Load GWAS data
            gwas_data = load_gwas_data(phecode, variant_id)
            print(f"GWAS data found: {gwas_data is not None}")

            # Construct record
            record = {
                "Variant": variant_id,
                "Gene": row["gene"],
                "AA_Change": row["aa"],
                "HMM_Score": row["hmm_score"],
            }

            # Add PheWAS data if available
            if phewas_data:
                record.update(phewas_data)

            # Add GWAS data if available
            if gwas_data:
                record.update(gwas_data)

            data_records.append(record)

        print(f"\nFinal number of records: {len(data_records)}")
        if len(data_records) > 0:
            print("Sample record:", data_records[0])

        # Define columns based on the data
        columns = ["Variant", "Gene", "AA_Change", "HMM_Score"]
        if any(phewas_data for r in data_records):
            columns.extend(["PheWAS_P", "PheWAS_OR", "Cases_Count", "Controls_Count"])
        if any(gwas_data for r in data_records):
            columns.extend(["GWAS_P", "GWAS_OR"])

        default_columns = ["Variant", "Gene", "AA_Change", "PheWAS_P", "GWAS_P"]
        num_columns = ["HMM_Score", "PheWAS_P", "PheWAS_OR", "GWAS_P", "GWAS_OR"]

        return jsonify(
            {
                "data": data_records,
                "columns": columns,
                "defaultColumns": default_columns,
                "numColumns": num_columns,
            }
        )

    except Exception as e:
        print("Error in tableVariantDetail:", traceback.format_exc())
        return jsonify({"error": str(e)}), 500


def get_term_domain_genes_variant(term):
    conn = get_db_connection()
    cur = conn.cursor()

    # Print the query and term value
    query = """YOUR QUERY HERE"""
    print(f"Executing query with term: {term}")
    print(f"Query: {query}")

    cur.execute(query, (term,))
    results = cur.fetchall()
    print(f"Raw results count: {len(results)}")

    # Rest of the function...
