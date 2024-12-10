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

    if not os.path.exists(phewas_path):
        # Run PheWAS
        # Convert variant_id to colon format for pheWAS
        variant_colon = variant_id.replace("_", ":").replace("/", ":")
        print(variant_colon)
        try:
            # Attempt to run PheWAS
            assoc = phecode_level_assoc(variant_colon)
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
    variant_gwas_format = variant_id.replace("_", ":", 1)

    variant_col = "nomaly_variant"

    row = gwas_df[gwas_df[variant_col] == variant_gwas_format]
    if row.empty:
        print(f"Variant {variant_gwas_format} not found in GWAS file")
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
    row = phewas_df[phewas_df["phecode"] == phecode]
    if row.empty:
        return None
    r = row.iloc[0]
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
    cur.execute(query, (term,))
    results = cur.fetchall()
    columns = [desc[0] for desc in cur.description]
    cur.close()
    conn.close()

    return pd.DataFrame(results, columns=columns) if results else pd.DataFrame()


@phecode_term_bp.route(
    "/phecode/<string:phecode>/term/<string:term>/tableVariantDetail",
    methods=["GET", "POST"],
)
def show_phecode_term_variant_detail(phecode, term):
    """
    Create a detailed view showing variants associated with the term,
    including allele frequencies, counts from PheWAS, PheWAS p-values, GWAS p-values, etc.
    """
    try:
        # Get variants from DB
        term_df = get_term_variants(term)
        if term_df.empty:
            return jsonify(
                {"data": [], "columns": [], "defaultColumns": [], "numColumns": []}
            )

        # For each variant, get PheWAS and GWAS results
        data_records = []
        for _, row in term_df.iterrows():
            variant_id = row["variant_id"]
            gene = row["gene"]
            aa = row["aa"]
            hmm_score = row["hmm_score"]

            # Extract alleles from variant_id. Format: CHR_POS_REF/ALT
            # Example: "8_6870776_C/T"
            # alleles: REF = C, ALT = T
            parts = variant_id.split("_")  # [CHR, POS, REF/ALT]
            if len(parts) != 3:
                print(f"Unexpected variant_id format: {variant_id}")
                continue
            chrom = parts[0]
            pos = parts[1]
            refalt = parts[2]
            if "/" in refalt:
                ref, alt = refalt.split("/")
            else:
                print(f"Unexpected ref/alt format in variant_id: {variant_id}")
                continue

            # Load PheWAS data for this variant and phecode
            phewas_data = load_phewas_data(variant_id, phecode)
            if phewas_data is None:
                # no data found from PheWAS
                phewas_data = {
                    "PheWAS_P": None,
                    "PheWAS_OR": None,
                    "Ref_AF_Cases": None,
                    "Alt_AF_Cases": None,
                    "Ref_AF_Controls": None,
                    "Alt_AF_Controls": None,
                    "Cases_Count": None,
                    "Controls_Count": None,
                }

            # Load GWAS data for this phecode and variant
            gwas_data = load_gwas_data(phecode, variant_id)
            if gwas_data is None:
                gwas_data = {
                    "GWAS_P": None,
                    "GWAS_OR": None,
                    "GWAS_F_A": None,
                    "GWAS_F_U": None,
                    "GWAS_RSID": None,
                }

            record = {
                "variant_id": variant_id,
                "gene": gene,
                "aa": aa,
                "hmm_score": f"{hmm_score:.4f}" if hmm_score is not None else None,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "phewas_p": f"{phewas_data['PheWAS_P']:.2e}"
                if phewas_data["PheWAS_P"]
                else None,
                "phewas_or": f"{phewas_data['PheWAS_OR']:.2f}"
                if phewas_data["PheWAS_OR"]
                else None,
                "cases_count": phewas_data["Cases_Count"],
                "controls_count": phewas_data["Controls_Count"],
                "ref_af_cases": f"{phewas_data['Ref_AF_Cases']:.4f}"
                if phewas_data["Ref_AF_Cases"]
                else None,
                "alt_af_cases": f"{phewas_data['Alt_AF_Cases']:.4f}"
                if phewas_data["Alt_AF_Cases"]
                else None,
                "ref_af_controls": f"{phewas_data['Ref_AF_Controls']:.4f}"
                if phewas_data["Ref_AF_Controls"]
                else None,
                "alt_af_controls": f"{phewas_data['Alt_AF_Controls']:.4f}"
                if phewas_data["Alt_AF_Controls"]
                else None,
                "gwas_p": f"{gwas_data['GWAS_P']:.2e}" if gwas_data["GWAS_P"] else None,
                "gwas_or": f"{gwas_data['GWAS_OR']:.2f}"
                if gwas_data["GWAS_OR"]
                else None,
                "gwas_f_a": f"{gwas_data['GWAS_F_A']:.4f}"
                if gwas_data["GWAS_F_A"]
                else None,
                "gwas_f_u": f"{gwas_data['GWAS_F_U']:.4f}"
                if gwas_data["GWAS_F_U"]
                else None,
                "gwas_rsid": gwas_data["GWAS_RSID"],
            }

            data_records.append(record)

        # Define columns and defaults
        columns = [
            "variant_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "gene",
            "aa",
            "hmm_score",
            "cases_count",
            "controls_count",
            "ref_af_cases",
            "alt_af_cases",
            "ref_af_controls",
            "alt_af_controls",
            "phewas_p",
            "phewas_or",
            "gwas_p",
            "gwas_or",
            "gwas_f_a",
            "gwas_f_u",
            "gwas_rsid",
        ]
        default_columns = [
            "variant_id",
            "gene",
            "aa",
            "hmm_score",
            "cases_count",
            "controls_count",
            "alt_af_cases",
            "alt_af_controls",
            "phewas_p",
            "phewas_or",
            "gwas_p",
            "gwas_or",
        ]
        num_columns = [
            "hmm_score",
            "alt_af_cases",
            "alt_af_controls",
            "ref_af_cases",
            "ref_af_controls",
            "phewas_p",
            "phewas_or",
            "gwas_p",
            "gwas_or",
            "gwas_f_a",
            "gwas_f_u",
        ]

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
