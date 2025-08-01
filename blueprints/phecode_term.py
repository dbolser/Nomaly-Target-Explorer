import logging

import numpy as np
from flask import Blueprint, current_app, jsonify, render_template, request, session

from blueprints.gwas import format_gwas_results, run_gwas

# TODO: MAKE A DATA SERVICE FOR THIS!!!!
from blueprints.nomaly import pharos, pp
from blueprints.phecode import get_phecode_data
from blueprints.phecode_term_helper import (
    CacheMiss,
    load_cache_results,
    save_cache_results,
)
from data_services import (
    GenotypeService,
    NomalyDataService,
    PhenotypeService,
    ServiceRegistry,
)
from db import (
    get_phecode_info,
    get_term_domains,
    get_term_genes,
    get_term_names,
    get_term_variants,
)

# # Create a 'dummy' profile decorator if we don't have line_profiler installed
# try:
#     from line_profiler import profile  # type: ignore
# except ImportError:

#     def profile(func):
#         return func


logger = logging.getLogger(__name__)
phecode_term_bp = Blueprint("phecode_term", __name__, template_folder="../templates")


@phecode_term_bp.route(
    "/phecode/<string:phecode>/term/<string:term>", methods=["POST", "GET"]
)
def show_phecode_term(phecode, term):
    # Get ancestry from session
    ancestry = session.get("ancestry", "EUR")

    try:
        # Get services while we're in 'app context'
        services: ServiceRegistry = current_app.extensions["nomaly_services"]
        phenotype_service = services.phenotype

        # NOTE: We inject the appropriate services into 'backend' functions
        # (dependency injection)

        phecode_term_data = get_phecode_info(phecode)

        case_counts = phenotype_service.get_case_counts_for_phecode(phecode, ancestry)

        phecode_term_data["population"] = ancestry
        phecode_term_data["affected"] = case_counts["affected"]
        phecode_term_data["excluded"] = case_counts["excluded"]
        phecode_term_data["control"] = case_counts["control"]

        phecode_term_data = get_phecode_data(phecode, services.phenotype, ancestry)

        term_name = get_term_names([term])[term]
        term_domains = get_term_domains([term])[term]
        term_gene_df = get_term_genes([term])

        # Update data dictionary
        phecode_term_data.update(
            {
                "term": term,
                "termname": term_name,
                "domainlen": len(term_domains),
                "genelen": term_gene_df["gene"].nunique(),
                "genes": ", ".join(term_gene_df["gene"].unique()),
            }
        )

        return render_template(
            "phecode_term.html", phecode=phecode, term=term, data=phecode_term_data
        )
    except Exception as e:
        logger.error(f"Failed to get Phecode Term data for {phecode} / {term}: {e}")
        return jsonify({"error": "Failed to get Phecode / Term data"}), 500


# Called by phecode_term.html
@phecode_term_bp.route(
    "/phecode/<string:phecode>/term/<string:term>/tableVariantDetail",
    methods=["GET", "POST"],
)
def show_phecode_term_variant_detail(
    phecode: str,
    term: str,
):
    flush = request.args.get("flush", "0") == "1"
    logger.info(f"Flush: {flush}")

    # Get ancestry from session
    ancestry = session.get("ancestry", "EUR")

    # Get services while we're in 'app context'
    services: ServiceRegistry = current_app.extensions["nomaly_services"]

    try:
        result = calculate_phecode_term_variant_detail(
            phecode,
            term,
            services.genotype,
            services.phenotype,
            services.nomaly_data,
            ancestry,
            flush,
        )
    except Exception as e:
        logger.error(f"Error in show_phecode_term_variant_detail: {str(e)}")
        return jsonify({"error": str(e)}), 500

    return jsonify(result)


# @profile
def calculate_phecode_term_variant_detail(
    phecode: str,
    term: str,
    genotype_service: GenotypeService,
    phenotype_service: PhenotypeService,
    nomaly_data_service: NomalyDataService,
    ancestry: str = "EUR",
    no_cache: bool = False,
) -> dict:

    logger.info(
        f"Getting variant details for phecode/term/ancestry: {phecode}/{term}/{ancestry} (flush: {no_cache})"
    )

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
        if not no_cache:
            try:
                cached_data = load_cache_results(phecode, term, ancestry)
                return {
                    "data": cached_data["data"],
                    "columns": columns,
                    "defaultColumns": default_columns,
                    "numColumns": numeric_columns,
                }

            except CacheMiss:
                logger.warning(
                    f"Cache miss for phecode {phecode}, term {term}, ancestry {ancestry}"
                )
                # continue...

        # Get variants from DB
        logger.info(f"Fetching variants for term: {term}")
        # NOTE: We're calling a function to get term -> variant mappings from
        # the Nomaly DB. Naturally, the IDs we get back are nomaly_variant_ids!
        term_df = get_term_variants(term)
        logger.debug(f"Initial term_df shape: {term_df.shape}")

        # Fill NA values and
        # term_df = term_df.fillna("None")

        # merge with pharos and pp data
        term_df = term_df.merge(pharos, on="gene", how="left").fillna("None")
        term_df = term_df.merge(pp, on="gene", how="left").fillna("None")
        logger.debug(f"After merges shape: {term_df.shape}")

        # Load GWAS data
        gwas_df = run_gwas(
            phecode, ancestry, phenotype_service, nomaly_data_service, no_cache
        )
        gwas_df = format_gwas_results(gwas_df)

        # TODO: JUST MERGE GWAS DF HERE!

        logger.info(f"Formatted GWAS data shape: {gwas_df.shape}")
        logger.info(f"Formatted GWAS columns: {gwas_df.columns.tolist()}")

        # TODO: CREATE A TEST SO WE CAN BE SURE THE FOLLOWING IS CORRECT...

        # counts = genotype_service.get_genotype_counts_and_freqs(
        #     vids=np.array(term_df["variant_id"].unique()),
        #     nomaly_ids=True,
        #     ancestry=ancestry,
        # )

        # merged_df = term_df.merge(counts, on="variant_id")
        # assert len(merged_df) == len(term_df)

        # merged_df["hmm_score_squared"] = merged_df["hmm_score"] ** 2
        # merged_df["vs00"] = (
        #     merged_df["hmm_score_squared"] * merged_df["het_freq"]
        #     + merged_df["hmm_score_squared"] * 4 * merged_df["alt_freq"]
        # )
        # merged_df["vs01"] = (
        #     merged_df["hmm_score_squared"] * merged_df["ref_freq"]
        #     + merged_df["hmm_score_squared"] * 4 * merged_df["alt_freq"]
        # )
        # merged_df["vs11"] = (
        #     merged_df["hmm_score_squared"] * merged_df["ref_freq"]
        #     + merged_df["hmm_score_squared"] * 4 * merged_df["alt_freq"]
        # )

        # # For reference to the above:
        # # hmm_squared = float(row["hmm_score"]) ** 2
        # # vs00 = hmm_squared * f01 + hmm_squared * 4 * f11
        # # vs11 = hmm_squared * f01 + hmm_squared * 4 * f00
        # # vs01 = hmm_squared * (f00 + f11)

        data_records = []
        for _, row in term_df.iterrows():
            nomaly_variant_id = str(row["variant_id"])
            logger.debug(f"\nProcessing variant: {nomaly_variant_id}")

            try:
                # Calculate genotype frequencies

                # TODO: USE THE NEW PHEWAS FUNCTION TO DO THIS IN ONE GO???
                # Actually, I think we can use the genotype_counts function in the
                # genotype service to do this in one go...
                # Yes we can! Just get counts outside this loop for all variants

                counts = genotype_service.get_genotype_counts_and_freqs(
                    vids=np.array([nomaly_variant_id]),
                    nomaly_ids=True,
                    ancestry=ancestry,
                )
                counts = counts.iloc[0].to_dict()

                f00 = counts["ref_freq"]
                f11 = counts["alt_freq"]
                f01 = counts["het_freq"]

                # TODO: Write tests for this!

                # Calculate variant scores...
                hmm_squared = float(row["hmm_score"]) ** 2
                vs00 = hmm_squared * f01 + hmm_squared * 4 * f11
                vs11 = hmm_squared * f01 + hmm_squared * 4 * f00
                vs01 = hmm_squared * (f00 + f11)

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
                    "hmoz_ref": counts["ref_count"],
                    "hmoz_alt": counts["alt_count"],
                    "htrz": counts["het_count"],
                    # Initialize GWAS fields with default values
                    "GWAS_P": "",
                    "GWAS_OR": "",
                    "GWAS_F_A": "",
                    "GWAS_F_U": "",
                    "GWAS_RSID": "",
                }

                # Add GWAS data
                # TODO: Merge this above (with the genotype counts)!!!

                if nomaly_variant_id in gwas_df["nomaly_variant"].values:
                    gwas_row = (
                        gwas_df[gwas_df["nomaly_variant"] == nomaly_variant_id]
                        .iloc[0]
                        .to_dict()
                    )

                    # TODO: Write tests for this!
                    record.update(
                        {
                            "GWAS_P": f"{gwas_row['P']:.2e}",
                            "GWAS_OR": f"{gwas_row['OR']:.2f}",
                            "GWAS_F_A": f"{gwas_row['F_A']:.5f}",
                            "GWAS_F_U": f"{gwas_row['F_U']:.5f}",
                            "GWAS_RSID": str(gwas_row["RSID"]),
                        }
                    )

                logger.debug(f"Final record GWAS fields: {record}")
                data_records.append(record)
            except (ValueError, TypeError) as e:
                logger.error(f"Error processing variant {nomaly_variant_id}: {e}")
                continue

        # Cache and return results
        # TODO: ANCESTRY
        save_cache_results(phecode, term, ancestry, data_records)
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
        raise e


def main():
    """Main function for testing."""

    from config import Config

    services = ServiceRegistry.from_config(Config)

    # phecode = "561"
    # phecode = "564.1"
    # phecode = "338"
    # phecode = "332"

    # term = "MP:0004957"
    # term = "HP:0000789"
    # term = "KW:0544"
    # term = "MP:0000948"
    # term = "MP:0004986"
    # term = "GO:0045244"

    phecode = "722.7"
    term = "GO:0045717"

    # gwas_data = run_gwas(
    #     phecode, "EUR", services.phenotype, services.nomaly_data, no_cache=True
    # )
    # print(gwas_data)

    term_data = get_term_variants(term)
    print(term_data)

    result = calculate_phecode_term_variant_detail(
        phecode,
        term,
        services.genotype,
        services.phenotype,
        services.nomaly_data,
        ancestry="EUR",
        no_cache=True,
    )
    print(result)

    exit(0)

    import numpy as np

    # Run some random GWAS for fun...
    ancestries = np.array(["AFR", "EAS", "EUR", "SAS"])
    phecodes = services.phenotype.phecodes

    done = set()
    for _ in range(10000):
        phecode = np.random.choice(phecodes)
        ancestry = np.random.choice(ancestries)

        terms = services.stats_registry.get("Run-v1", ancestry)._hdf.terms
        term = np.random.choice(terms)

        if (phecode, ancestry, term) in done:
            continue
        done.add((phecode, ancestry, term))

        print(f"Running VP for {phecode} / {term} ({ancestry})")
        try:
            calculate_phecode_term_variant_detail(
                phecode,
                term,
                services.genotype,
                services.phenotype,
                services.nomaly_data,
                ancestry,
                no_cache=True,
            )
            print(f"Successfully ran VP for {phecode} / {term} ({ancestry})")
        except Exception as e:
            print(f"Error running VP for {phecode} / {term} ({ancestry}): {e}")
            continue


if __name__ == "__main__":
    main()
