import math
from flask import Blueprint, render_template, jsonify, request

import pandas as pd
import os
import traceback

from blueprints.gwas import variant_level_assoc, GWAS_PHENO_DIR

from db import get_term_domains, get_term_names, get_term_genes, get_phecode_info
from db import get_term_domain_genes_variant, get_term_domain_genes
from db import get_term_variants

from blueprints.nomaly import pharos, pp
from blueprints.nomaly import nomaly_stats

from blueprints.nomaly import nomaly_genotype

from errors import DataNotFoundError, GWASError
import logging

# from blueprints.phewas import get_formatted_phewas_data

from blueprints.phecode_term_helper import load_cached_results, save_results


# Create the blueprint
phecode_term_bp = Blueprint("phecode_term", __name__, template_folder="../templates")

logger = logging.getLogger(__name__)


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


@phecode_term_bp.route(
    "/phecode/<string:phecode>/term/<string:term>/tableVariantDetail",
    methods=["GET", "POST"],
)
def show_phecode_term_variant_detail(phecode: str, term: str, flush: bool = False):
    # Define the columns
    columns = [
        "Variant",
        "Gene",
        "AA_Change",
        "HMM_Score",
        "TDL",
        "TBIO",
        "Classification",
        "Drug_Program_Indication",
        "GWAS_P",
        "GWAS_OR",
        "GWAS_F_A",
        "GWAS_F_U",
        "GWAS_RSID",
        "vs00",
        "vs01",
        "vs11",
        "hmoz_alt",
        "hmoz_ref",
        "htrz",
    ]

    default_columns = [
        "Variant",
        "Gene",
        "Classification",
        "AA_Change",
        "HMM_Score",
        "GWAS_P",
        "GWAS_OR",
        "GWAS_RSID",
        "vs00",
        "vs01",
        "vs11",
    ]

    numeric_columns = ["HMM_Score", "GWAS_P", "GWAS_OR", "vs00", "vs01", "vs11"]

    try:
        # Get flush parameter from POST body or URL query parameter
        flush = request.args.get("flush", "false").lower() == "true"  # Handle URL parameter
        if request.is_json:
            flush = request.get_json().get("flush", flush)  # POST body overrides URL parameter if present
        logger.info(f"Flush parameter received: {flush}")

        # First check cache
        cached_data = load_cached_results(phecode, term, flush)
        
        if cached_data is not None and "data" in cached_data:
            logger.info(f"Using cached data for phecode {phecode}, term {term}")
            return jsonify({
                "data": cached_data["data"],
                "columns": columns,
                "defaultColumns": default_columns,
                "numColumns": numeric_columns,
            })
        
        logger.warning(f"Cache miss or flush requested for phecode {phecode}, term {term}")

        # If no cache, get variants from DB
        print(f"Fetching variants for term: {term}")
        term_df = get_term_variants(term)
        print(f"Initial term_df shape: {term_df.shape}")

        # Fill NA values before merging
        term_df = term_df.fillna("None")

        # Add pharos data
        term_df = term_df.merge(pharos, on="gene", how="left")
        print(f"After pharos merge shape: {term_df.shape}")

        term_df = term_df.fillna("None")  # Fill NAs after first merge

        # Add pp data
        term_df = term_df.merge(pp, on="gene", how="left")
        print(f"After pp merge shape: {term_df.shape}")
        term_df = term_df.fillna("None")  # Fill NAs after second merge

        if term_df.empty:
            print("No variants found for this term!")
            return jsonify(
                {"data": [], "columns": [], "defaultColumns": [], "numColumns": []}
            )

        # Load GWAS data for all variants
        gwas_data = load_gwas_data(phecode)
        print(f"GWAS data found: {gwas_data is not None}")

        # If P is NaN, set to 1
        gwas_data["P"] = gwas_data["P"].fillna(1)
        gwas_data["OR"] = gwas_data["OR"].fillna(1)

        # Load genotype counts for all variants
        genotype_counts = nomaly_genotype.get_variant_counts()

        data_records = []

        # Calculate 'the score'...
        for idx, row in term_df.iterrows():
            variant_id = row["variant_id"]
            variant_id = variant_id.replace("_", ":")
            variant_id = variant_id.replace("/", ":")

            if variant_id not in genotype_counts.index:
                # Swap the alleles (the last two parts of the variant_id)
                alleles = variant_id.split(":")
                alleles[-1], alleles[-2] = alleles[-2], alleles[-1]
                variant_id = ":".join(alleles)
                if variant_id not in genotype_counts.index:
                    print(f"Variant {variant_id} not found in genotype counts")
                    continue

            hmoz_alt = genotype_counts.loc[variant_id, "homozygous_alt"]
            htrz = genotype_counts.loc[variant_id, "heterozygous"]
            hmoz_ref = genotype_counts.loc[variant_id, "homozygous_ref"]
            total = hmoz_alt + hmoz_ref + htrz

            f00 = hmoz_alt / total
            f01 = htrz / total
            f11 = hmoz_ref / total

            HMM_score = row["hmm_score"]

            term_df.at[idx, "vs00"] = HMM_score * f01 + (HMM_score**2 * f11)
            term_df.at[idx, "vs01"] = HMM_score**2 * (f00 + f01)
            term_df.at[idx, "vs11"] = HMM_score * f01 + (HMM_score**2 * f00)
            term_df.at[idx, "hmoz_alt"] = hmoz_alt
            term_df.at[idx, "hmoz_ref"] = hmoz_ref
            term_df.at[idx, "htrz"] = htrz

        for idx, row in term_df.iterrows():
            variant_id = row["variant_id"]
            print(f"\nProcessing variant {idx+1}/{len(term_df)}: {variant_id}")
            # standard_variant_id = variant_id.replace("/", "_")

            # Construct record
            record = {
                "Variant": variant_id,
                "Gene": row["gene"],
                "AA_Change": str(row["aa"]),
                "HMM_Score": f"{float(row['hmm_score']):.2f}",
                "TDL": str(row.get("tdl", "None")),
                "TBIO": str(row.get("tbio", "None")),
                "Classification": str(row.get("classification", "None")),
                "Drug_Program_Indication": str(
                    row.get("drug_program_indication", "None")
                ),
                "vs00": f"{float(row['vs00']):.2f}",
                "vs01": f"{float(row['vs01']):.2f}",
                "vs11": f"{float(row['vs11']):.2f}",
                "hmoz_alt": row["hmoz_alt"],
                "hmoz_ref": row["hmoz_ref"],
                "htrz": row["htrz"],
            }

            # Add GWAS data if available
            row = gwas_data[gwas_data["nomaly_variant"] == variant_id]

            if row.empty:
                print(f"Variant {variant_id} not found in GWAS file")

            else:
                r = row.iloc[0]
                record.update(
                    {
                        "GWAS_P": f"{r.get('P', 1):.2e}",
                        "GWAS_OR": f"{r.get('OR', 1):.2f}",
                        "GWAS_F_A": f"{r.get('F_A', 1):.5f}",
                        "GWAS_F_U": f"{r.get('F_U', 1):.5f}",
                        "GWAS_RSID": r.get("RSID"),
                    }
                )

            data_records.append(record)

        # Cache the results
        save_results(phecode, term, data_records)

        print(f"\nFinal number of records: {len(data_records)}")
        if len(data_records) > 0:
            print("Sample record:", data_records[0])

        result = {
            "data": data_records,
            "columns": columns,
            "defaultColumns": default_columns,
            "numColumns": numeric_columns,
        }

        if request:
            return jsonify(result)
        else:
            return pd.DataFrame(result["data"])

    except Exception as e:
        error_msg = f"Error processing phecode {phecode}, term {term}: {str(e)}"
        print(error_msg)
        print(traceback.format_exc())
        return jsonify({"error": error_msg}), 500


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


def load_gwas_data(phecode):
    """
    Load GWAS data for a given phecode and return row for this variant.
    Returns a dict with keys: GWAS_P, GWAS_OR, etc.
    If not found, returns None.
    """
    ensure_gwas(phecode)
    gwas_file = f"{GWAS_PHENO_DIR}phecode_{phecode}.assoc_nomaly.tsv"
    if not os.path.exists(gwas_file):
        print(f"GWAS file not found for phecode {phecode}")
        return None
    gwas_df = pd.read_csv(gwas_file, sep="\t")

    return gwas_df


def main():
    print("Hello")

    term = "MP:0004957"
    term_df = get_term_variants(term)
    print(term_df.head())

    phecode = "282"
    # Load GWAS data for all variants
    gwas_data = load_gwas_data(phecode)
    print(f"GWAS data found: {gwas_data is not None}")

    result = show_phecode_term_variant_detail(phecode, term, flush=True)
    print(result)

    result.iloc[0]

    print(f"Hello")


if __name__ == "__main__":
    main()
