from flask import Blueprint, render_template, jsonify, request

import pandas as pd
import traceback

from blueprints.gwas import run_gwas, format_gwas_results

from db import (
    get_term_variants,
    get_term_domains,
    get_term_names,
    get_term_genes,
    get_phecode_info,
)

from blueprints.nomaly import pharos, pp

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
    logger.info(
        f"Starting variant detail processing for phecode {phecode}, term {term}"
    )
    logger.debug(f"Flush parameter: {flush}")

    # Define the columns
    columns = [
        "Variant",
        "Gene",
        "AA_Change",
        "HMM_Score",
        "Classification",
        "Drug_Program_Indication",
        "TDL",
        "TBIO",
        "vs00",
        "vs01",
        "vs11",
        "hmoz_alt",
        "hmoz_ref",
        "htrz",
        "GWAS_P",
        "GWAS_OR",
        "GWAS_RSID",
    ]
    logger.debug(f"Defined {len(columns)} columns")

    default_columns = [
        "Variant",
        "Gene",
        "AA_Change",
        "HMM_Score",
        "Classification",
        "vs00",
        "vs01",
        "vs11",
        "GWAS_P",
        "GWAS_OR",
        "GWAS_RSID",
    ]

    numeric_columns = ["HMM_Score", "vs00", "vs01", "vs11", "GWAS_P", "GWAS_OR"]

    try:
        # Get flush parameter from POST body or URL query parameter
        flush = (
            request.args.get("flush", "false").lower() == "true"
        )  # Handle URL parameter
        if request.is_json:
            flush = request.get_json().get(
                "flush", flush
            )  # POST body overrides URL parameter if present
        logger.info(f"Flush parameter received: {flush}")

        # Check cache
        cached_data = load_cached_results(phecode, term, flush)

        if cached_data is not None and "data" in cached_data:
            logger.info(f"Using cached data for phecode {phecode}, term {term}")
            return jsonify(
                {
                    "data": cached_data["data"],
                    "columns": columns,
                    "defaultColumns": default_columns,
                    "numColumns": numeric_columns,
                }
            )

        logger.warning(
            f"Cache miss or flush requested for phecode {phecode}, term {term}"
        )

        # Get variants from DB
        print(f"Fetching variants for term: {term}")
        term_df = get_term_variants(term)
        print(f"Initial term_df shape: {term_df.shape}")

        # Fill NA values and merge with pharos and pp data
        term_df = term_df.fillna("None")
        term_df = term_df.merge(pharos, on="gene", how="left").fillna("None")
        term_df = term_df.merge(pp, on="gene", how="left").fillna("None")
        print(f"After merges shape: {term_df.shape}")

        if term_df.empty:
            print("No variants found for this term!")
            return jsonify(
                {
                    "data": [],
                    "columns": columns,
                    "defaultColumns": default_columns,
                    "numColumns": numeric_columns,
                }
            )

        # Load GWAS data
        gwas_data = run_gwas(phecode)
        logger.info(f"GWAS data found: {not gwas_data.empty}")
        logger.info(f"GWAS data columns: {gwas_data.columns.tolist()}")
        logger.info(f"First few rows of GWAS data:\n{gwas_data.head()}")
        print(f"GWAS data found: {not gwas_data.empty}", flush=True)
        print(f"GWAS data columns: {gwas_data.columns.tolist()}", flush=True)
        print(f"First few rows of GWAS data:\n{gwas_data.head()}", flush=True)

        # Fill missing values in GWAS data
        gwas_data["P"] = gwas_data["P"].fillna(1)
        gwas_data["OR"] = gwas_data["OR"].fillna(1)

        # Load genotype counts
        genotype_counts = nomaly_genotype.get_variant_counts()

        data_records = []
        for idx, row in term_df.iterrows():
            variant_id = row["variant_id"].replace("_", ":").replace("/", ":")
            logger.info(f"\nProcessing variant: {variant_id}")
            print(f"\nProcessing variant: {variant_id}", flush=True)

            # Try both allele orientations for genotype counts
            if variant_id not in genotype_counts.index:
                alleles = variant_id.split(":")
                alleles[-1], alleles[-2] = alleles[-2], alleles[-1]
                variant_id = ":".join(alleles)
                if variant_id not in genotype_counts.index:
                    print(f"Variant {variant_id} not found in genotype counts")
                    continue

            # Calculate genotype frequencies
            counts = genotype_counts.loc[variant_id]
            total = int(counts.sum())
            f00 = float(counts["homozygous_alt"]) / total
            f01 = float(counts["heterozygous"]) / total
            f11 = float(counts["homozygous_ref"]) / total

            # Calculate variant scores
            HMM_score = float(row["hmm_score"])
            vs00 = HMM_score * f01 + (HMM_score**2 * f11)
            vs01 = HMM_score**2 * (f00 + f01)
            vs11 = HMM_score * f01 + (HMM_score**2 * f00)

            # Build record with explicit type conversion
            record = {
                "Variant": variant_id,
                "Gene": str(row["gene"]),
                "AA_Change": str(row["aa"]),
                "HMM_Score": f"{float(row['hmm_score']):.2f}",
                "TDL": str(row.get("tdl", "None")),
                "TBIO": str(row.get("tbio", "None")),
                "Classification": str(row.get("classification", "None")),
                "Drug_Program_Indication": str(
                    row.get("drug_program_indication", "None")
                ),
                "vs00": f"{vs00:.2f}",
                "vs01": f"{vs01:.2f}",
                "vs11": f"{vs11:.2f}",
                "hmoz_alt": int(counts["homozygous_alt"]),
                "hmoz_ref": int(counts["homozygous_ref"]),
                "htrz": int(counts["heterozygous"]),
                # Initialize GWAS fields with default values
                "GWAS_P": "",
                "GWAS_OR": "",
                "GWAS_F_A": "",
                "GWAS_F_U": "",
                "GWAS_RSID": "",
            }

            # Add GWAS data if available
            # Convert variant format from 1:23519074:A:G to 1_23519074_A/G
            gwas_variant_id = variant_id.replace(
                ":", "_", 2
            )  # Replace first two colons with underscore
            gwas_variant_id = gwas_variant_id.replace(
                ":", "/"
            )  # Replace last colon with slash
            print(f"\nLooking for variant: {variant_id}")
            print(f"Converted to GWAS format: {gwas_variant_id}")

            gwas_row = gwas_data[gwas_data["nomaly_variant"] == gwas_variant_id]

            print(f"GWAS data variants: {gwas_data['nomaly_variant'].head()}")
            print(f"Found match: {not gwas_row.empty}")

            logger.info(f"Found GWAS data for variant: {not gwas_row.empty}")
            logger.info(
                f"GWAS row data:\n{gwas_row.to_dict('records')}"
                if not gwas_row.empty
                else "No GWAS data found"
            )
            print(f"Found GWAS data for variant: {not gwas_row.empty}", flush=True)
            print(
                f"GWAS row data:\n{gwas_row.to_dict('records')}"
                if not gwas_row.empty
                else "No GWAS data found",
                flush=True,
            )

            record.update(
                {
                    "GWAS_P": f"{float(gwas_row.iloc[0].get('P', 1)):.2e}"
                    if not gwas_row.empty
                    else "",
                    "GWAS_OR": f"{float(gwas_row.iloc[0].get('OR', 1)):.2f}"
                    if not gwas_row.empty
                    else "",
                    "GWAS_F_A": f"{float(gwas_row.iloc[0].get('F_A', 1)):.5f}"
                    if not gwas_row.empty
                    else "",
                    "GWAS_F_U": f"{float(gwas_row.iloc[0].get('F_U', 1)):.5f}"
                    if not gwas_row.empty
                    else "",
                    "GWAS_RSID": str(gwas_row.iloc[0].get("RSID", ""))
                    if not gwas_row.empty
                    else "",
                }
            )
            logger.info(f"Final record GWAS fields: {record}")
            print(f"Final record GWAS fields: {record}", flush=True)

            data_records.append(record)

        # Cache and return results
        save_results(phecode, term, data_records)
        print(f"\nFinal number of records: {len(data_records)}")

        result = {
            "data": data_records,
            "columns": columns,
            "defaultColumns": default_columns,
            "numColumns": numeric_columns,
        }

        print(f"Hello")  # Test print at end

        return jsonify(result) if request else pd.DataFrame(result["data"])

    except Exception as e:
        error_msg = f"Error processing phecode {phecode}, term {term}: {str(e)}"
        print(error_msg)
        print(traceback.format_exc())
        return jsonify({"error": error_msg}), 500


def main():
    print("Hello")

    term = "MP:0004957"
    term_df = get_term_variants(term)
    print(term_df.head())

    phecode = "282"
    # Load GWAS data for all variants
    gwas_data = run_gwas(phecode)
    print(f"GWAS data found: {gwas_data is not None}")

    result = show_phecode_term_variant_detail(phecode, term, flush=True)
    print(result)

    result.iloc[0]

    print(f"Hello")


if __name__ == "__main__":
    main()
