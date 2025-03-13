"""
This module contains refactored functions for prioritizing variants by Nomaly scores.
The original mega-function has been broken down into smaller, more focused functions.
"""

import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

# Create a logger
logger = logging.getLogger(__name__)

# Import profile decorator (keeping the same behavior)
try:
    from line_profiler import profile  # type: ignore
except ImportError:

    def profile(func):
        """Dummy profile decorator when line_profiler is not available."""
        return func


class StreamLogger:
    """Helper class to capture and stream progress messages."""

    def __init__(self, message_queue):
        self.message_queue = message_queue

    def info(self, message):
        """Send a progress message immediately to the queue."""
        msg = {"type": "progress", "data": message}
        self.message_queue.put(msg)


def get_cache_path(disease_code: str, term: str) -> Path:
    """Get the path to the cache file for the given disease code and term."""
    # Implementation as in original code
    pass


def load_cached_results(disease_code: str, term: str) -> dict:
    """Load cached results for the given disease code and term."""
    # Implementation as in original code
    pass


def save_results_to_cache(disease_code: str, term: str, data: dict):
    """Save results to cache for the given disease code and term."""
    # Implementation as in original code
    pass


def get_term_variants(term: str) -> List[str]:
    """Get list of variant IDs for a given term."""
    # Implementation as in original code
    pass


def term_variant_prioritisation(vids, eids, genotype_service, stream_logger=None):
    """Process variants using the thread pool."""
    # Implementation as in original code
    pass


def fetch_phenotype_data(
    phecode: str,
    phenotype_service,
    population: Optional[str] = None,
    biological_sex: Optional[str] = None,
    stream_logger=None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Fetch phenotype data and separate cases, controls, and excludes.

    Args:
        phecode: The disease code
        phenotype_service: Service to fetch phenotype data
        population: Optional population filter
        biological_sex: Optional biological sex filter
        stream_logger: Optional logger for streaming progress

    Returns:
        Tuple containing (case_eids, control_eids, exclude_eids, all_eids, phenotypes)
    """
    # Get phenotype data for the given phecode
    eids, phenotypes = phenotype_service.get_cases_for_phecode(
        phecode, population, biological_sex
    )

    log_message = f"Got {len(eids)} samples for phecode {phecode} with population {population} and sex {biological_sex}"
    if stream_logger:
        stream_logger.info(log_message)
    else:
        logger.info(log_message)

    # Separate cases, controls, and excludes
    case_eids = eids[phenotypes == 1]
    control_eids = eids[phenotypes == 0]
    exclude_eids = eids[phenotypes == 9]

    # Verify that all eids are accounted for
    assert len(eids) == len(case_eids) + len(control_eids) + len(exclude_eids)

    log_message = f"Got {len(case_eids)} cases for {phecode} with population {population} and sex {biological_sex}"
    if stream_logger:
        stream_logger.info(log_message)
    else:
        logger.info(log_message)

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
    stats["metric1_fp"] = np.sum(control_scores < stats["metric1_threshold"])
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
    term_variant_ids: List[str],
    genotype_service,
    phecode: str,
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
        term_variant_ids: Variant IDs for the term
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
        term_variant_ids,
        eids_above_threshold,
        genotype_service,
        stream_logger,
    )

    # Generate gene-level statistics
    gene_data = process_gene_level_stats(top_variants)

    # Add results to stats
    stats[f"{stat}_top_variants"] = top_variants.to_dict(orient="records")
    stats[f"{stat}_top_gene_set"] = gene_data

    return stats


def process_gene_level_stats(top_variants: pd.DataFrame) -> List[Dict[str, Any]]:
    """
    Process gene-level statistics from top variants.

    Args:
        top_variants: DataFrame of top variants

    Returns:
        List of dictionaries with gene-level statistics
    """
    # Explode genes (one variant might be associated with multiple genes)
    top_variants_per_gene = top_variants.explode("gene")

    # Format gene lists for display
    top_variants_copy = top_variants.copy()
    top_variants_copy["gene"] = top_variants_copy["gene"].map(lambda x: ", ".join(x))

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

    # Convert to records
    return (
        top_gene_set.reset_index()
        .rename(columns={"index": "gene"})
        .to_dict(orient="records")
    )


@profile
def get_top_variants(
    phecode: str,
    term: str,
    phenotype_service,
    genotype_service,
    nomaly_scores_service,
    stats_service,
    population: Optional[str] = None,
    biological_sex: Optional[str] = None,
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
        cached_results = load_cached_results(phecode, term)
        if cached_results:
            if stream_logger:
                stream_logger.info("Loaded results from cache")
            return cached_results

    # Log that we're computing results
    if stream_logger:
        stream_logger.info(
            f"Computing results (not found in cache, flush={no_cache}, protective={protective})"
        )

    try:
        # Step 1: Fetch phenotype data
        case_eids, control_eids, exclude_eids, all_eids, phenotypes = (
            fetch_phenotype_data(
                phecode, phenotype_service, population, biological_sex, stream_logger
            )
        )

        # Step 2: Get term variants
        term_variant_ids = get_term_variants(term)
        if stream_logger:
            stream_logger.info(f"Got {len(term_variant_ids)} variants for term {term}")

        # Step 3: Fetch Nomaly scores
        case_scores, control_scores = fetch_nomaly_scores(
            case_eids, control_eids, nomaly_scores_service, term
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
                term_variant_ids,
                genotype_service,
                phecode,
                stream_logger,
            )

        # Step 8: Save to cache
        save_results_to_cache(phecode, term, stats)

        return stats

    except Exception as e:
        logger.error(f"Error in get_top_variants: {str(e)}", exc_info=True)
        if stream_logger:
            stream_logger.info(f"Error in analysis: {str(e)}")
        raise
