import logging
import traceback
from typing import Optional

import pandas as pd
from flask import Blueprint, current_app, jsonify, render_template, request

from blueprints.gwas import format_gwas_results, run_gwas
from blueprints.nomaly import pharos, pp
from blueprints.phecode_term_helper import load_cached_results, save_results
from db import (
    get_phecode_info,
    get_term_domains,
    get_term_genes,
    get_term_names,
    get_term_variants,
)
from errors import DataNotFoundError

# Create a 'dummy' profile decorator if we don't have line_profiler installed
try:
    from line_profiler import profile
except ImportError:

    def profile(func):
        return func


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
def show_phecode_term_variant_detail(
    phecode: str,
    term: str,
    sex: Optional[str] = None,
    ancestry: Optional[str] = None,
    flush: bool = False,
):
    flush = request.args.get("flush", "false").lower() == "true"  # Handle URL parameter
    if request.is_json:
        flush = request.get_json().get(
            "flush", flush
        )  # POST body overrides URL parameter if present

    try:
        result = calculate_phecode_term_variant_detail(
            phecode, term, sex, ancestry, flush
        )
    except Exception as e:
        logger.error(f"Error in show_phecode_term_variant_detail: {str(e)}")
        return jsonify({"error": str(e)}), 500

    return jsonify(result)


# @profile
def calculate_phecode_term_variant_detail(
    phecode: str,
    term: str,
    sex: Optional[str] = None,
    ancestry: Optional[str] = None,
    flush: bool = False,
) -> dict:
    services = current_app.extensions["nomaly_services"]
    genotype_service = services.genotype
    assert genotype_service is not None

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
        logger.info(f"Flush parameter received: {flush}")

        # Check cache
        cached_data = load_cached_results(phecode, term, flush)

        if cached_data is not None and "data" in cached_data:
            logger.info(f"Using cached data for phecode {phecode}, term {term}")
            return {
                "data": cached_data["data"],
                "columns": columns,
                "defaultColumns": default_columns,
                "numColumns": numeric_columns,
            }

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

        # Load GWAS data
        gwas_data = run_gwas(phecode)
        formatted_gwas = pd.DataFrame(
            format_gwas_results(gwas_data, significance_threshold=0.1)
        )
        if not formatted_gwas.empty:
            logger.info(f"Formatted GWAS data shape: {formatted_gwas.shape}")
            logger.info(f"Formatted GWAS columns: {formatted_gwas.columns.tolist()}")

        data_records = []
        for _, row in term_df.iterrows():
            nomaly_variant_id = str(row["variant_id"])
            logger.info(f"\nProcessing variant: {nomaly_variant_id}")

            try:
                # Calculate genotype frequencies
                counts = genotype_service._hdf.get_variant_counts(
                    nomaly_variant_id=nomaly_variant_id
                )
                total = counts["total"]
                f00 = float(counts["homozygous_ref"]) / total
                f11 = float(counts["homozygous_alt"]) / total
                f01 = float(counts["heterozygous"]) / total

                # Calculate variant scores
                HMM_score = float(row["hmm_score"])
                hmm2 = HMM_score * HMM_score

                vs00 = hmm2 * f01 + hmm2 * 4 * f11
                vs11 = hmm2 * f01 + hmm2 * 4 * f00
                vs01 = hmm2 * (f00 + f11)

                # Build record with explicit type conversion
                record = {
                    "Variant": nomaly_variant_id,
                    "Gene": str(row["gene"]),
                    "AA_Change": str(row["aa"]),
                    "HMM_Score": f"{float(row['hmm_score']):.2f}",
                    "TDL": str(row.get("tdl", "None")),
                    "TBIO": str(row.get("tbio", "None")),
                    "Classification": str(row.get("classification", "None")),
                    "Drug_Program_Indication": str(
                        row.get("drug_program_indication", "None")
                    ),
                    "vs00": f"{vs00:.4f}",
                    "vs01": f"{vs01:.4f}",
                    "vs11": f"{vs11:.4f}",
                    "hmoz_ref": counts["homozygous_ref"],
                    "hmoz_alt": counts["homozygous_alt"],
                    "htrz": counts["heterozygous"],
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
                        formatted_gwas["nomaly_variant"] == nomaly_variant_id
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
                logger.error(f"Error processing variant {nomaly_variant_id}: {e}")
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

        return result

    except Exception as e:
        error_msg = f"Error processing phecode {phecode}, term {term}: {str(e)}"
        logger.error(error_msg)
        logger.error(traceback.format_exc())
        raise e


def main():
    phecode = "561"
    phecode = "564.1"
    phecode = "338"

    term = "MP:0004957"
    term = "HP:0000789"
    term = "KW:0544"
    term = "MP:0000948"

    from app import create_app

    app = create_app("development")

    with app.app_context():
        gwas_data = run_gwas(phecode)
        print(gwas_data)

        term_data = get_term_variants(term)
        print(term_data)

        result = calculate_phecode_term_variant_detail(phecode, term, flush=True)
        print(result)


if __name__ == "__main__":
    main()
