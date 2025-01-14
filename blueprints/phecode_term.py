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

from errors import DataNotFoundError
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
        formatted_gwas = pd.DataFrame(
            format_gwas_results(gwas_data, significance_threshold=0.1)
        )
        if not formatted_gwas.empty:
            logger.info(f"Formatted GWAS data shape: {formatted_gwas.shape}")
            logger.info(f"Formatted GWAS columns: {formatted_gwas.columns.tolist()}")

        # Load genotype counts
        genotype_counts = nomaly_genotype.get_variant_counts()

        data_records = []
        for _, row in term_df.iterrows():
            variant_id = str(row["variant_id"])
            genotype_variant_id = variant_id.replace("_", ":").replace("/", ":")
            logger.info(f"\nProcessing variant: {variant_id} ({genotype_variant_id})")

            # Try both allele orientations for genotype counts
            if genotype_variant_id not in genotype_counts.index:
                alleles = genotype_variant_id.split(":")
                alleles[-1], alleles[-2] = alleles[-2], alleles[-1]
                genotype_variant_id = ":".join(alleles)
                if genotype_variant_id not in genotype_counts.index:
                    logger.warning(f"Variant {variant_id} not found in genotype counts")
                    continue

            try:
                # Calculate genotype frequencies
                counts = genotype_counts.loc[genotype_variant_id]
                total = int(float(counts.sum()))
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
                    "hmoz_alt": int(float(counts["homozygous_alt"])),
                    "hmoz_ref": int(float(counts["homozygous_ref"])),
                    "htrz": int(float(counts["heterozygous"])),
                    # Initialize GWAS fields with default values
                    "GWAS_P": "",
                    "GWAS_OR": "",
                    "GWAS_F_A": "",
                    "GWAS_F_U": "",
                    "GWAS_RSID": "",
                }

                # Add GWAS data if available

                if not formatted_gwas.empty:
                    gwas_row = formatted_gwas[
                        formatted_gwas["nomaly_variant"] == variant_id
                    ]
                    if not gwas_row.empty:
                        gwas_data = gwas_row.iloc[0].to_dict()
                        record.update(
                            {
                                "GWAS_P": f"{float(gwas_data.get('P', 1.0)):.2e}",
                                "GWAS_OR": f"{float(gwas_data.get('OR', 1.0)):.2f}",
                                "GWAS_F_A": f"{float(gwas_data.get('F_A', 0.0)):.5f}",
                                "GWAS_F_U": f"{float(gwas_data.get('F_U', 0.0)):.5f}",
                                "GWAS_RSID": str(gwas_data.get("RSID", "")),
                            }
                        )

                logger.debug(f"Final record GWAS fields: {record}")
                data_records.append(record)
            except (ValueError, TypeError) as e:
                logger.error(f"Error processing variant {variant_id}: {e}")
                continue

        # Cache and return results
        save_results(phecode, term, data_records)
        logger.info(f"\nFinal number of records: {len(data_records)}")

        result = {
            "data": data_records,
            "columns": columns,
            "defaultColumns": default_columns,
            "numColumns": numeric_columns,
        }

        return jsonify(result)

    except Exception as e:
        error_msg = f"Error processing phecode {phecode}, term {term}: {str(e)}"
        logger.error(error_msg)
        logger.error(traceback.format_exc())
        return jsonify({"error": error_msg}), 500


def main():
    phecode = "561"
    phecode = "564.1"

    term = "MP:0004957"
    term = "HP:0000789"
    term = "KW:0544"

    gwas_data = run_gwas(phecode)
    term_data = get_term_variants(term)

    print(f"GWAS data found: {gwas_data is not None}")

    result = show_phecode_term_variant_detail(phecode, term, flush=True)

    print(result)


if __name__ == "__main__":
    main()
