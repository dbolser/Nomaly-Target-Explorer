"""
This module contains refactored functions for prioritizing variants by Nomaly scores.
The original mega-function has been broken down into smaller, more focused functions.
"""

import json
import logging
from typing import Any, Dict, Optional, Sequence, Tuple, cast

import numpy as np
import pandas as pd


from blueprints.prioritisation_by_nomaly_scores import (
    convert_numpy_types,
    load_cached_results,
    save_results_to_cache,
    read_nomaly_filtered_genotypes_new,
    process_individual_variants,
    variant_scores,
    log_and_stream,
)

from db import get_term_variants


# Create a logger
logger = logging.getLogger(__name__)


def term_variant_prioritisation(
    vids, eids, ancestry, genotype_service, stream_logger=None
) -> pd.DataFrame:
    """Prioritise a set of variants from a set of eids and their genotypes.

    Typically we're looking at the set of variants for a given term
    (term-variants) and a set of 'high scoring' *cases* (True Positives).

    Genotypes are 'scored' using background data from Nomaly (HMM Score and
    Nomaly Score for the given eids).
    """

    # Remove the duplicate aa column but keep the largest hmm_score
    vids = (
        vids.drop(columns=["aa"])
        .sort_values(by="hmm_score", ascending=False)
        .drop_duplicates(subset=["variant_id", "gene"], keep="first")
        .reset_index(drop=True)
    )

    # Group by variant_id, colapse duplicate genes into a list and select the largest hmm_score
    term_variants_genes = vids.groupby("variant_id")["gene"].apply(list).reset_index()
    term_variants_hmmsc = vids.groupby("variant_id")["hmm_score"].max().reset_index()
    vids = term_variants_genes.merge(term_variants_hmmsc, on="variant_id")

    log_and_stream(f"Reading genotypes for {len(vids)} variants", stream_logger)

    sel_genotypes = read_nomaly_filtered_genotypes_new(
        eids, vids["variant_id"], ancestry, genotype_service
    )

    log_and_stream(
        f"Genotypes for {len(sel_genotypes['row_eids'])} individuals and {len(sel_genotypes['col_variants'])} variants are read.",
        stream_logger,
    )
    if len(sel_genotypes["error_variants"]) > 0:
        log_and_stream(
            f"Failed to get genotype data for {len(sel_genotypes['error_variants'])} variants.",
            stream_logger,
            level="warning",
        )

    variant_scores_table = variant_scores.loc[sel_genotypes["col_variants"]]
    term_variant_scores = variant_scores_table[["vs00", "vs01", "vs11"]]

    print(
        f"We read {len(term_variant_scores)} variant scores for our {len(vids)} variants"
    )

    # TODO:This is now the bottleneck!
    ind_top_variants_df = process_individual_variants(
        sel_genotypes, term_variant_scores
    )

    print(
        f"For some reason, we get just {len(ind_top_variants_df)} top variants for {len(sel_genotypes['row_eids'])} individuals"
    )

    top_variants = vids.join(ind_top_variants_df, on="variant_id", how="right")
    top_variants = top_variants.join(term_variant_scores, on="variant_id")

    # Convert gene lists to strings to match the original output format
    top_variants["gene"] = top_variants["gene"].map(lambda x: ", ".join(x))

    log_and_stream(f"Term variants: {len(vids)}", stream_logger)
    log_and_stream(f"Top variants: {len(top_variants)}", stream_logger)
    log_and_stream(
        f"Individual top variants: {len(ind_top_variants_df)}", stream_logger
    )

    return top_variants


def fetch_phenotype_data(
    phecode: str,
    phenotype_service,
    population: str = "EUR",
    stream_logger=None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Fetch phenotype data and separate cases, controls, and excludes.

    Args:
        phecode: The disease code
        phenotype_service: Service to fetch phenotype data
        population: Optional population filter, defaults to "EUR"
        stream_logger: Optional logger for streaming progress

    Returns:
        Tuple containing (case_eids, control_eids, exclude_eids, all_eids, phenotypes)
    """
    # Get phenotype data for the given phecode
    eids, phenotypes = phenotype_service.get_cases_for_phecode(phecode, population)

    log_and_stream(
        f"Got {len(eids)} samples for phecode {phecode} with population {population}",
        stream_logger,
    )

    # Separate cases, controls, and excludes
    case_eids = eids[phenotypes == 1]
    control_eids = eids[phenotypes == 0]
    exclude_eids = eids[phenotypes == 9]

    # Verify that all eids are accounted for
    assert len(eids) == len(case_eids) + len(control_eids) + len(exclude_eids)

    log_and_stream(
        f"Got {len(case_eids)} cases for {phecode} with population {population}",
        stream_logger,
    )

    # Assert that the case eids are sorted
    assert np.all(np.diff(case_eids) > 0), "case_eids are not sorted!"

    return case_eids, control_eids, exclude_eids, eids, phenotypes


def fetch_nomaly_scores(
    case_eids: np.ndarray, control_eids: np.ndarray, nomaly_scores_service, term: str
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fetch Nomaly scores for cases and controls.

    Args:
        case_eids: EIDs for cases
        control_eids: EIDs for controls
        nomaly_scores_service: Service to fetch Nomaly scores
        term: The term to get scores for

    Returns:
        Tuple containing (case_scores, control_scores)
    """
    case_scores = nomaly_scores_service.get_scores_by_eids_unsorted(
        case_eids, terms=[term]
    )
    assert len(case_eids) == len(case_scores)

    control_scores = nomaly_scores_service.get_scores_by_eids_unsorted(
        control_eids, terms=[term]
    )
    assert len(control_eids) == len(control_scores)

    return case_scores, control_scores


def fetch_stats_data(term: str, phecode: str, stats_service) -> Dict[str, Any]:
    """
    Fetch statistics for the given term and phecode.

    Args:
        term: Term to get stats for
        phecode: Disease code to get stats for
        stats_service: Service to fetch stats

    Returns:
        Dictionary of statistics
    """
    stats = stats_service.get_stats_by_term_phecode(
        term=term,
        phecode=phecode,
        statstype=[
            "num_rp",
            "num_rn",
            "metric1_pvalue",
            "roc_stats_mcc_pvalue",
            "roc_stats_mcc_or",
            "roc_stats_mcc_threshold",
            "roc_stats_mcc_all_index",
            "roc_stats_mcc_tp",
            "roc_stats_mcc_fp",
            "roc_stats_mcc_fn",
            "roc_stats_mcc_tn",
            "roc_stats_yjs_pvalue",
            "roc_stats_yjs_or",
            "roc_stats_yjs_threshold",
            "roc_stats_yjs_tp",
            "roc_stats_yjs_fp",
            "roc_stats_yjs_fn",
            "roc_stats_yjs_tn",
            "roc_stats_lrp_pvalue",
            "roc_stats_lrp_or",
            "roc_stats_lrp_threshold",
            "roc_stats_lrp_tp",
            "roc_stats_lrp_fp",
            "roc_stats_lrp_fn",
            "roc_stats_lrp_tn",
            "roc_stats_lrn_protective_pvalue",
            "roc_stats_lrn_protective_or",
            "roc_stats_lrn_protective_threshold",
            "roc_stats_lrn_protective_tp",
            "roc_stats_lrn_protective_fp",
            "roc_stats_lrn_protective_fn",
            "roc_stats_lrn_protective_tn",
        ],
    )

    return clean_stats_data(stats)


def clean_stats_data(stats: Dict[str, Any]) -> Dict[str, Any]:
    """
    Clean statistics data by handling NaN and Inf values.

    Args:
        stats: Statistics dictionary

    Returns:
        Cleaned statistics dictionary
    """
    # Make a copy to avoid modifying the original
    stats = stats.copy()

    # Handle NaN in pvalue
    if pd.isna(stats.get("metric1_pvalue")):
        stats["metric1_pvalue"] = 1

    # Clean up other values
    for key in stats.keys():
        if key.endswith("_threshold"):
            if stats[key] == float("inf"):
                stats[key] = 0
        if key.endswith("_or"):
            if np.isnan(stats[key]):
                stats[key] = 1

    return stats


def compute_derived_stats(
    stats: Dict[str, Any],
    case_eids: np.ndarray,
    control_eids: np.ndarray,
    case_scores: np.ndarray,
    control_scores: np.ndarray,
    phecode: str,
) -> Dict[str, Any]:
    """
    Compute derived statistics from base stats.

    Args:
        stats: Base statistics
        case_eids: EIDs for cases
        control_eids: EIDs for controls
        case_scores: Nomaly scores for cases
        control_scores: Nomaly scores for controls
        phecode: Disease code

    Returns:
        Updated statistics dictionary
    """
    # Make a copy to avoid modifying the original
    stats = stats.copy()

    # Check for mismatches between phenotype and stats services
    if len(case_eids) != stats["num_rp"] or len(control_eids) != stats["num_rn"]:
        logger.warning(
            f"Phecode {phecode} has {len(case_eids)} / {len(control_eids)} cases and controls"
            + f" from phenotype service but {stats['num_rp']} / {stats['num_rn']} cases and controls"
            + " from stats service!"
        )

    # Fill in the 'missing' values for metric1
    stats["metric1_threshold"] = 0.022

    # Compute true/false positives/negatives for metric1
    stats["metric1_tp"] = np.sum(case_scores >= stats["metric1_threshold"])
    stats["metric1_fp"] = np.sum(control_scores >= stats["metric1_threshold"])
    stats["metric1_fn"] = len(case_eids) - stats["metric1_tp"]
    stats["metric1_tn"] = len(control_eids) - stats["metric1_fp"]

    return stats


def process_statistic(
    stat: str,
    stats: Dict[str, Any],
    case_eids: np.ndarray,
    control_eids: np.ndarray,
    case_scores: np.ndarray,
    control_scores: np.ndarray,
    term_variants: pd.DataFrame,  # Changed from List[str] to match actual usage
    genotype_service,
    phecode: str,
    ancestry: str,
    stream_logger=None,
) -> Dict[str, Any]:
    """
    Process a single statistic (e.g., metric1, roc_stats_mcc).

    Args:
        stat: Name of statistic to process
        stats: Statistics dictionary
        case_eids: EIDs for cases
        control_eids: EIDs for controls
        case_scores: Nomaly scores for cases
        control_scores: Nomaly scores for controls
        term_variants: DataFrame of variants for the term
        genotype_service: Service to fetch genotype data
        phecode: Disease code
        stream_logger: Optional logger for streaming progress

    Returns:
        Updated statistics dictionary with computed values for this statistic
    """
    # Make a copy to avoid modifying the original
    stats = stats.copy()

    # Handle edge case where threshold is 0 and tp is 0
    if stats[f"{stat}_threshold"] <= 0 and stats[f"{stat}_tp"] == 0:
        stats[f"{stat}_threshold"] = 1e308

    # Get case EIDs above threshold
    if stat.endswith("protective"):
        eids_above_threshold = control_eids[
            control_scores >= stats[f"{stat}_threshold"]
        ]
    else:
        eids_above_threshold = case_eids[case_scores >= stats[f"{stat}_threshold"]]

    if stream_logger:
        stream_logger.info(
            f"Processing statistic {stat} with {len(eids_above_threshold)} individuals"
        )
    print(f"Selected for statistic {stat}: {len(eids_above_threshold)}")

    # Check for mismatches between computed TP and stats service TP
    if len(eids_above_threshold) != stats[f"{stat}_tp"]:
        logger.warning(
            f"Phecode {phecode} has {len(eids_above_threshold)} True Positives for {stat} from "
            + f"nomaly scores but {stats[f'{stat}_tp']} True Positives from stats service!"
        )

    # Calculate derived metrics
    stats[f"{stat}_tpr"] = stats[f"{stat}_tp"] / stats["num_rp"]
    stats[f"{stat}_fpr"] = stats[f"{stat}_fp"] / stats["num_rn"]
    stats[f"{stat}_lrp"] = stats[f"{stat}_tpr"] / (stats[f"{stat}_fpr"] + 1e-10)

    # Prioritize variants
    top_variants = term_variant_prioritisation(
        term_variants,
        eids_above_threshold,
        ancestry,
        genotype_service,
        stream_logger,
    )

    # Generate gene-level statistics
    gene_data = process_gene_level_stats(top_variants)

    # Add results to stats
    stats[f"{stat}_top_variants"] = top_variants.to_dict(orient="records")
    stats[f"{stat}_top_gene_set"] = gene_data

    return stats


def process_gene_level_stats(top_variants: pd.DataFrame) -> Sequence[Dict[str, Any]]:
    """
    Process gene-level statistics from top variants.

    Args:
        top_variants: DataFrame of top variants

    Returns:
        Sequence of dictionaries with gene-level statistics
    """
    # Explode genes (one variant might be associated with multiple genes)
    # First, convert any string genes to lists for consistent processing
    if top_variants["gene"].dtype == object and isinstance(
        top_variants["gene"].iloc[0], str
    ):
        # If genes are already strings (e.g., "GENE1, GENE2"), split them
        top_variants_copy = top_variants.copy()
        top_variants_copy["gene"] = top_variants_copy["gene"].str.split(", ")
        top_variants_per_gene = top_variants_copy.explode("gene")
    else:
        # If genes are already lists, just explode them
        top_variants_per_gene = top_variants.explode("gene")

    # Aggregate by gene
    gene_variants = top_variants_per_gene.groupby("gene")["variant_id"].agg(list)
    gene_hmm_score = top_variants_per_gene.groupby("gene")["hmm_score"].sum()
    gene_total_vs = top_variants_per_gene.groupby("gene")["vs"].sum()
    gene_num_individual = top_variants_per_gene.groupby("gene")["num_individuals"].sum()
    gene_num_variants = top_variants_per_gene.groupby("gene").size()

    # Create gene-level dataframe
    top_gene_set = pd.DataFrame(
        {
            "variant_id": gene_variants.map(lambda x: ", ".join(x)),
            "hmm_score": gene_hmm_score,
            "total_vs": gene_total_vs,
            "variant_num": gene_num_variants,
            "num_individuals": gene_num_individual,
        }
    )

    # Sort by hmm_score
    top_gene_set = top_gene_set.sort_values("hmm_score", ascending=False)

    # Convert to records - using Sequence type to fix type error
    return cast(
        Sequence[Dict[str, Any]],
        top_gene_set.reset_index()
        .rename(columns={"index": "gene"})
        .to_dict(orient="records"),
    )


def get_top_variants_refactored(
    phecode: str,
    term: str,
    phenotype_service,
    genotype_service,
    score_service,
    stats_service,
    run_version: str,
    ancestry: str,
    stream_logger=None,
    protective: bool = False,
    no_cache: bool = False,
) -> Dict[str, Any]:
    """
    Get the top variants for the given disease code and term.

    This refactored version breaks down the original function into smaller,
    more focused functions with clear responsibilities.

    Args:
        phecode: Disease code to analyze
        term: Term to analyze
        phenotype_service: Service to fetch phenotype data
        genotype_service: Service to fetch genotype data
        nomaly_scores_service: Service to fetch Nomaly scores
        stats_service: Service to fetch statistics
        population: Optional population filter
        biological_sex: Optional biological sex filter
        stream_logger: Optional logger for streaming progress
        protective: Whether to include protective variants
        no_cache: Whether to skip cache lookup

    Returns:
        Dictionary with analysis results
    """
    # Try to load from cache first (unless no_cache is True)
    if not no_cache:
        cached_results = load_cached_results(phecode, term, run_version, ancestry)
        if cached_results:
            if stream_logger:
                stream_logger.info("Loaded results from cache")
            return cached_results

    # Compute results if not cached...
    log_and_stream(
        f"Computing results (not found in cache, flush={no_cache}, protective={protective})",
        stream_logger,
    )

    try:
        # Step 1: Fetch phenotype data
        case_eids, control_eids, exclude_eids, all_eids, phenotypes = (
            fetch_phenotype_data(phecode, phenotype_service, ancestry, stream_logger)
        )

        # Step 2: Get term variants
        term_variants = get_term_variants(term)
        log_and_stream(
            f"Got {len(term_variants)} variants for term {term}", stream_logger
        )

        # Step 3: Fetch Nomaly scores
        case_scores, control_scores = fetch_nomaly_scores(
            case_eids, control_eids, score_service, term
        )

        # Step 4: Fetch stats data
        stats = fetch_stats_data(term, phecode, stats_service)

        # Step 5: Compute derived statistics
        stats = compute_derived_stats(
            stats, case_eids, control_eids, case_scores, control_scores, phecode
        )

        # Step 6: Determine which statistics to process
        stats_to_process = [
            "metric1",
            "roc_stats_mcc",
            "roc_stats_yjs",
            "roc_stats_lrp",
        ]
        if protective:
            stats_to_process.append("roc_stats_lrn_protective")

        # Step 7: Process each statistic
        for stat in stats_to_process:
            stats = process_statistic(
                stat,
                stats,
                case_eids,
                control_eids,
                case_scores,
                control_scores,
                term_variants,
                genotype_service,
                phecode,
                ancestry,
                stream_logger,
            )

        # Step 8: Save to cache
        save_results_to_cache(phecode, term, stats, run_version, ancestry)

        return stats

    except Exception as e:
        logger.error(f"Error in get_top_variants: {str(e)}", exc_info=True)
        raise


def main():
    """This code exists for debugging purposes."""

    # disease_code = "290.11"
    # term = "GO:0016861"

    # disease_code = "290.11"
    # term = "CD:MESH:D009139"

    # disease_code = "334.2"
    # term = "GO:0044559"

    # disease_code = "334.2"
    # term = "GO:0044559"

    # disease_code = "324.1"
    # term = "GO:0009225"

    phecode = "332"
    term = "GO:0030800"
    term = "MP:0004986"

    # disease_code = "300.13"
    # term = "GO:1901136"

    # disease_code = "334.2"
    # term = "CD:MESH:D009139"

    # disease_code = "705"
    # term = "GO:0034081"
    # term = "GO:0003960"

    phecode = "256"
    term = "MP:0004819"

    from app import create_app

    # TODO: Avoid using app context when we should just be able to 'inject' services...
    app = create_app("development")
    with app.app_context():
        services = app.extensions["nomaly_services"]

        run_version = "Run-v1"
        ancestry = "EUR"

        data = get_top_variants_refactored(
            phecode,
            term,
            services.phenotype,
            services.genotype,
            services.nomaly_score,
            services.stats,
            run_version=run_version,
            ancestry=ancestry,
            no_cache=True,
            # protective=True,
        )

        # Search for values that will cause problems in JavaScript JSON parsing
        def check_json_safety(obj):
            if isinstance(obj, dict):
                for key, value in obj.items():
                    if value is None:
                        print(f"Found None value for key: {key}")
                        obj[key] = ""  # Replace None with empty string
                    else:
                        check_json_safety(value)
            elif isinstance(obj, list):
                for i, item in enumerate(obj):
                    if item is None:
                        print(f"Found None value for index: {i}")
                        obj[i] = ""  # Replace None with empty string
                    else:
                        check_json_safety(item)
            return obj

        # Clean up any None/null values that would break JavaScript
        data = check_json_safety(data)

        # print(json.dumps(data, indent=2, sort_keys=True))

        for key, value in data.items():
            if "top_variants" in key or "top_gene_set" in key:
                data[key] = len(data[key])

        data = convert_numpy_types(data)

        print(json.dumps(data, indent=2, sort_keys=True))

    print("OK")


if __name__ == "__main__":
    main()
