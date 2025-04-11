"""
Prioritisation of individual variants by Nomaly prediction scores and disease status.
per disease and term.

Genes are prioritised by the sum of Nomaly scores of their variants.
"""

import json
import logging
import traceback
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from queue import Queue
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from flask import (
    Blueprint,
    Response,
    current_app,
    render_template,
    request,
    session,
    stream_with_context,
)

from blueprints.phecode import get_phecode_data
from config import Config
from data_services import (
    GenotypeService,
    NomalyScoreService,
    PhenotypeService,
    ServiceRegistry,
    StatsService,
)
from db import get_term_domains, get_term_genes, get_term_names, get_term_variants

logger = logging.getLogger(__name__)


from config import Config

# # Create a 'dummy' profile decorator if we don't have line_profiler installed
# try:
#     from line_profiler import profile  # type: ignore
# except ImportError:

#     def profile(func):
#         return func


# TODO: Move to a new variant scores (or 'nomaly_data') data service?

variant_scores = pd.read_csv(
    Config.NOMALY_VARIANT_SCORES_PATH,
    sep="\t",
    index_col="nomaly_variant_id",
)

variant2gene = pd.read_csv(
    Config.NOMALY_VARIANT2GENE_PATH, sep="\t", header=None, index_col=0
)

# Global thread pool for variant processing
# Limit concurrent processing to avoid overwhelming the server
MAX_WORKERS = 20  # Adjust based on server capacity
variant_processor = ThreadPoolExecutor(
    max_workers=MAX_WORKERS, thread_name_prefix="variant_processor"
)


# Convert NumPy types to native Python types
def convert_numpy_types(obj: Any) -> Any:
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    return obj


class StreamLogger:
    """Helper class to capture and stream progress messages.

    This is basically a Flask thing as far as I know.
    """

    def __init__(self, message_queue: Queue):
        self.message_queue = message_queue

    def log(self, message, level="info"):
        """Send a message immediately to the queue (if defined)."""

        if self.message_queue:
            self.message_queue.put({"type": "progress", "data": message})
        else:
            # Not sure how we could end up on this path...
            getattr(logger, level)(message)


# @profile
def individual_variant_prioritisation(row, term_variant_scores) -> pd.DataFrame:
    """Return numpy array of variant scores for the selected variants."""
    variants = term_variant_scores["variant_id"].values
    sel_vs = term_variant_scores[["vs00", "vs01", "vs11"]].to_numpy()

    geno_matrix = np.zeros((len(row), 3))
    for i, val in enumerate(row):
        if val != -9:
            geno_matrix[i, int(val)] = 1

    scores = np.dot(geno_matrix, sel_vs.T)
    diagonal_scores = np.diagonal(scores)

    sorted_indices = np.argsort(diagonal_scores)[::-1]
    sorted_scores = diagonal_scores[sorted_indices]
    sorted_variants = variants[sorted_indices]

    top_variants = pd.DataFrame({"variant_id": sorted_variants, "vs": sorted_scores})
    top_variants = top_variants.assign(rank=range(1, len(top_variants) + 1))
    top_variants = top_variants[(top_variants["vs"] > 1) | (top_variants["rank"] <= 5)]

    return top_variants


# @profile
def process_individual_variants(genotypes, term_variant_scores):
    """Process variants for all individuals more efficiently...

    More efficiently than what?
    """
    ind_top_variants = defaultdict(int)
    ind_top_variants_scores = defaultdict(float)

    genotype_data = genotypes.T

    for individual in genotype_data:
        top_variants = individual_variant_prioritisation(
            individual, term_variant_scores
        )
        for _, row in top_variants.iterrows():
            variant = row["variant_id"]
            vs = row["vs"]
            ind_top_variants[variant] += 1
            ind_top_variants_scores[variant] += vs

    return pd.DataFrame(
        {
            "variant_id": list(ind_top_variants.keys()),
            "num_individuals": list(ind_top_variants.values()),
            "vs": list(ind_top_variants_scores.values()),
        }
    )


# @profile
def term_variant_prioritisation(
    term_variants: pd.DataFrame,
    case_eids_above_threshold: np.ndarray,
    genotype_service: GenotypeService,
    stream_logger: StreamLogger,
) -> pd.DataFrame:
    """Prioritise a set of variants from a set of eids and their genotypes.

    Typically we're looking at the set of variants for a given term
    (term-variants) and a set of 'high scoring' *cases* (True Positives).

    Genotypes are 'scored' using background data from Nomaly (HMM Score and
    Nomaly Score for the given eids).
    """

    # TODO: All this term_variant dataframe hacking can be done back in
    # 'get_top_variants' outside the 'for each stat in stats' loop.

    # Remove the duplicate aa column but keep the largest hmm_score
    term_variants = (
        term_variants.drop(columns=["aa"])
        .sort_values(by="hmm_score", ascending=False)
        .drop_duplicates(subset=["variant_id", "gene"], keep="first")
        .reset_index(drop=True)
    )

    # Group by variant_id, colapse duplicate genes into a list and select the largest hmm_score
    term_variants_genes = (
        term_variants.groupby("variant_id")["gene"].apply(list).reset_index()
    )
    term_variants_hmmsc = (
        term_variants.groupby("variant_id")["hmm_score"].max().reset_index()
    )
    term_variants = term_variants_genes.merge(term_variants_hmmsc, on="variant_id")

    stream_logger.log(
        f"Reading genotypes for {len(term_variants)} variants across {len(case_eids_above_threshold)} individuals",
    )

    vids_array = term_variants["variant_id"].to_numpy()

    genotypes = genotype_service.get_genotypes(
        eids=case_eids_above_threshold, vids=vids_array, nomaly_ids=True
    )

    genotype_counts = genotype_service.get_genotype_counts_and_freqs(
        eids=case_eids_above_threshold,
        vids=vids_array,
        nomaly_ids=True,
    )

    # This step just serves to align the variant_id columns
    genotype_counts_and_hmm_scores = genotype_counts.merge(
        term_variants_hmmsc, on="variant_id"
    )

    f00 = genotype_counts_and_hmm_scores["ref_freq"]
    f11 = genotype_counts_and_hmm_scores["alt_freq"]
    f01 = genotype_counts_and_hmm_scores["het_freq"]

    hmm2 = genotype_counts_and_hmm_scores["hmm_score"] ** 2

    vs00 = hmm2 * f01 + hmm2 * 4 * f11
    vs11 = hmm2 * f01 + hmm2 * 4 * f00
    vs01 = hmm2 * (f00 + f11)

    # Create term_variant_scores DataFrame with variant_id as a column, not index
    term_variant_scores = pd.DataFrame(
        {
            "variant_id": genotype_counts_and_hmm_scores["variant_id"],
            "vs00": vs00,
            "vs01": vs01,
            "vs11": vs11,
        }
    )

    # Process individual variants with the calculated scores
    sel_genotypes = {
        "row_eids": case_eids_above_threshold,
        "col_variants": vids_array,
        "data": genotypes.T,
        "error_variants": [],
    }

    ind_top_variants_df = process_individual_variants(genotypes, term_variant_scores)

    stream_logger.log(
        f"Found {len(ind_top_variants_df)} top variants for {len(case_eids_above_threshold)} individuals",
    )

    # Use merge instead of join since ind_top_variants_df now has variant_id as a column
    top_variants = term_variants.merge(
        ind_top_variants_df, on="variant_id", how="right"
    )
    # Join term_variant_scores as well, but using merge since it now has variant_id as a column
    top_variants = top_variants.merge(term_variant_scores, on="variant_id", how="left")

    stream_logger.log(
        f"Term variants: {len(term_variants)}"
        f"Top variants: {len(top_variants)}"
        f"Individual top variants: {len(ind_top_variants_df)}"
    )

    return top_variants


def get_cache_path(
    disease_code: str, term: str, run_version: str, ancestry: str
) -> Path:
    """Get the cache file path for this disease/term combination."""
    cache_dir = Path(Config.VARIANT_SCORES_DIR)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return (
        cache_dir
        / f"variant_prioritization_{disease_code}_{term}_{run_version}_{ancestry}.json"
    )


def load_cached_results(
    disease_code: str, term: str, run_version: str, ancestry: str
) -> Optional[Dict[str, Any]]:
    """Load cached results if they exist."""
    cache_path = get_cache_path(disease_code, term, run_version, ancestry)
    if not cache_path.exists():
        return None

    try:
        # Load the cached JSON data
        with open(cache_path) as f:
            data = json.load(f)

        return data
    except Exception as e:
        logger.error(f"Error loading cached results: {e}")
        return None


def save_results_to_cache(
    disease_code: str, term: str, data: dict, run_version: str, ancestry: str
) -> Optional[Dict[str, Any]]:
    """Save results to cache file."""
    try:
        cache_path = get_cache_path(disease_code, term, run_version, ancestry)

        # Convert data before saving
        converted_data = convert_numpy_types(data)

        # Save as JSON
        with open(cache_path, "w") as f:
            json.dump(converted_data, f, indent=2)  # indent for readability
        # Return the loaded data
        return load_cached_results(disease_code, term, run_version, ancestry)
    except Exception as e:
        logger.error(f"Error saving results to cache: {e}")
        return None


# @profile
def get_top_variants(
    phecode: str,
    term: str,
    phenotype_service: PhenotypeService,
    genotype_service: GenotypeService,
    nomaly_scores_service: NomalyScoreService,
    stats_service: StatsService,
    run_version: str,
    ancestry: str,
    stream_logger: StreamLogger,
    protective: bool = False,
    no_cache: bool = False,
) -> dict:
    """Get the top variants for the disease and term."""

    if not no_cache:
        try:
            cached_results = load_cached_results(phecode, term, run_version, ancestry)
            if cached_results:
                stream_logger.log("Loaded results from cache")
                return cached_results
        except Exception as e:
            logger.error(f"Error loading cached results: {e}")

    stream_logger.log(
        f"Computing results (not found in cache, flush={no_cache}, protective={protective})",
    )

    # Get phenotype data for the given phecode...
    phenotypes = phenotype_service.get_cases_for_phecode(phecode, ancestry)

    stream_logger.log(
        f"Got {len(phenotypes)} samples for phecode {phecode} with ancestry {ancestry}",
    )

    case_eids = phenotypes[phenotypes["phenotype"] == 1]["eid"].to_numpy()
    ctrl_eids = phenotypes[phenotypes["phenotype"] == 0]["eid"].to_numpy()
    excl_eids = phenotypes[phenotypes["phenotype"] == 9]["eid"]

    assert len(phenotypes) == len(case_eids) + len(ctrl_eids) + len(excl_eids)

    stream_logger.log(
        f"Got {len(case_eids)} cases for {phecode} with ancestry {ancestry}",
    )

    # TODO: This should be conditioned on the run_version
    case_scores = nomaly_scores_service.get_scores_by_eids_unsorted(
        case_eids, terms=np.array([term])
    )

    assert len(case_eids) == len(case_scores)

    stream_logger.log(
        f"Got {len(ctrl_eids)} controls for {phecode} with ancestry {ancestry}",
    )

    # TODO: This should be conditioned on the run_version
    ctrl_scores = nomaly_scores_service.get_scores_by_eids_unsorted(
        ctrl_eids, terms=np.array([term])
    )

    assert len(ctrl_eids) == len(ctrl_scores)

    stream_logger.log(
        f"Got {len(excl_eids)} excluded samples for {phecode} with ancestry {ancestry}",
    )

    try:
        # NOTE: The stats service is specific to run_version and ancestry
        stats = stats_service.get_term_stats(term=term, phecode=phecode)
    except Exception as e:
        logger.error(f"Error getting stats for {term} and {phecode}: {e}")
        raise

    # Quick cleanup of stats...

    # Just ignore these!
    if np.any(stats[stats.filter(like="_threshold").columns] == float("inf")):
        stats[stats.filter(like="_threshold").columns] = stats.filter(
            like="_threshold"
        ).replace(float("inf"), 0)
    if np.any(stats[stats.filter(like="_or").columns].isna()):
        stats[stats.filter(like="_or").columns] = stats.filter(like="_or").fillna(1)

    # Convert stats to dict..
    stats = stats.to_dict(orient="records")[0]

    # Err..
    if pd.isna(stats["metric1_pvalue"]):
        stats["metric1_pvalue"] = 1

    # TODO: Something is wrong somewhere...
    if len(case_eids) != stats["num_rp"]:
        logger.warning(
            f"""Phecode {phecode} has {len(case_eids)} cases
            from phenotype service but {stats["num_rp"]} cases
            from stats service!"""
        )

    if len(ctrl_eids) != stats["num_rn"]:
        logger.info(
            f"""Phecode {phecode} has {len(ctrl_eids)} controls
            from phenotype service but {stats["num_rn"]} controls
            from stats service!
            I THINK THIS IS TODO WITH EXCLUDED SEX-MATCHED CONTROLS!"""
        )

    term_variants = get_term_variants(term)

    stream_logger.log(f"Got {len(term_variants)} variants for term {term}")

    # Fill in the 'missing' values for metric1

    stats["metric1_threshold"] = 0.022

    # It's hard being this dumb...

    # True positives are cases above the threshold
    stats["metric1_tp"] = np.sum(case_scores >= stats["metric1_threshold"])

    # False positives are controls above the threshold
    stats["metric1_fp"] = np.sum(ctrl_scores >= stats["metric1_threshold"])

    # False negatives are cases below the threshold (or, all the other cases)
    stats["metric1_fn"] = len(case_eids) - stats["metric1_tp"]

    # True negatives are controls below the threshold (or, all the other controls)
    stats["metric1_tn"] = len(ctrl_eids) - stats["metric1_fp"]

    stats_to_process = [
        "metric1",
        "roc_stats_mcc",
        "roc_stats_yjs",
        "roc_stats_lrp",
    ]

    if protective:
        stats_to_process.append("roc_stats_lrn_protective")

    for stat in stats_to_process:
        stream_logger.log(f"Processing statistic {stat}")

        # How did this happen?
        if stats[f"{stat}_threshold"] <= 0 and stats[f"{stat}_tp"] == 0:
            # I think this should have been set this way?
            stats[f"{stat}_threshold"] = 1e308

        case_eids_above_threshold = case_eids[case_scores >= stats[f"{stat}_threshold"]]

        if stat.endswith("protective"):
            case_eids_above_threshold = ctrl_eids[
                ctrl_scores >= stats[f"{stat}_threshold"]
            ]

        stream_logger.log(
            f"Selected for statistic {stat}: {len(case_eids_above_threshold)}",
        )

        if len(case_eids_above_threshold) != stats[f"{stat}_tp"]:
            # position_in_index = str(stats.get(f"{stat}_all_index", "Unknown"))
            logger.warning(
                f"Phecode {phecode} has {len(case_eids_above_threshold)} True Positives for {stat} from "
                + f"nomaly scores but {stats[f'{stat}_tp']} True Positives from stats service!"
                # TODO: The index looks wrong.
                # + f"Note that the index is {position_in_index}."
            )

        # Calculate True Positive Rate (Sensitivity (TP/(TP+FN) or TP/RP))
        # Calculate False Positive Rate (1 - Specificity (FP/(FP+TN) or FP/RN))
        stats[f"{stat}_tpr"] = stats[f"{stat}_tp"] / stats["num_rp"]
        stats[f"{stat}_fpr"] = stats[f"{stat}_fp"] / stats["num_rn"]

        # Calculate Positive Likelihood Ratio (TPR/FPR) (TPR/(1-TNR))
        stats[f"{stat}_lrp"] = stats[f"{stat}_tpr"] / (stats[f"{stat}_fpr"] + 1e-10)

        top_variants = term_variant_prioritisation(
            term_variants,
            case_eids_above_threshold,
            genotype_service,
            stream_logger,
        )

        # Tables

        # So do we explode gene here first or what?
        top_variants_per_gene = top_variants.explode("gene")

        # Only now can we format these fuckers
        top_variants["gene"] = top_variants["gene"].map(lambda x: ", ".join(x))

        gene_variants = top_variants_per_gene.groupby("gene")["variant_id"].agg(list)
        gene_hmm_score = top_variants_per_gene.groupby("gene")["hmm_score"].sum()
        gene_total_vs = top_variants_per_gene.groupby("gene")["vs"].sum()
        gene_num_individual = top_variants_per_gene.groupby("gene")[
            "num_individuals"
        ].sum()
        gene_num_variants = top_variants_per_gene.groupby("gene").size()

        top_gene_set = pd.DataFrame(
            {
                "variant_id": gene_variants.map(lambda x: ", ".join(x)),
                "hmm_score": gene_hmm_score,
                "total_vs": gene_total_vs,
                "variant_num": gene_num_variants,
                "num_individuals": gene_num_individual,
            }
        )
        # top_gene_set["avg_individuals_per_variant"] = top_gene_set["num_individuals"] / top_gene_set["variant_num"]
        top_gene_set = top_gene_set.sort_values("hmm_score", ascending=False)

        # Create gene set data without tuple wrapping
        gene_set_data = (
            top_gene_set.reset_index()
            .rename(columns={"index": "gene"})
            .to_dict(orient="records")
        )

        # Do this last
        stats[f"{stat}_top_variants"] = top_variants.to_dict(orient="records")
        stats[f"{stat}_top_gene_set"] = gene_set_data

    # Save to cache
    save_results_to_cache(phecode, term, stats, run_version, ancestry)

    return stats


prioritisation_bp = Blueprint("prioritisation", __name__)


@prioritisation_bp.route("/variant_scores/<disease_code>/<term>")
def show_variant_scores(
    disease_code: str,
    term: str,
):
    """Show the variant scores page"""
    flush = bool(request.args.get("flush", False))
    protective = bool(request.args.get("protective", False))
    print(f"Flush: {flush}, Protective: {protective}")

    services = current_app.extensions["nomaly_services"]
    phenotype_service = services.phenotype

    # Get run version and ancestry from the session
    # TODO: Implement run version!
    run_version = session.get("run_version", "Run-v1")
    ancestry = session.get("ancestry", "EUR")

    # Get phecode data
    phecode_data = get_phecode_data(disease_code, phenotype_service, ancestry)

    # Get term data
    term_names = get_term_names([term])
    term_domains = get_term_domains([term])
    term_genes = get_term_genes([term])

    # Add term details to 'page data'
    phecode_data["term"] = term
    phecode_data["termname"] = term_names.get(term, "")
    phecode_data["domainlen"] = len(term_domains.get(term, []))

    genes = term_genes[term_genes["term"] == term]["gene"].tolist()
    phecode_data["genelen"] = len(genes)
    phecode_data["genes"] = (
        ", ".join(genes) if len(genes) < 5 else f"{len(genes)} genes"
    )

    return render_template(
        "variant_scores.html",
        disease_code=disease_code,
        term=term,
        data=phecode_data,
        flush=flush,  # Pass to template
        protective=protective,
        top_variants=pd.DataFrame(),
        top_gene_set=pd.DataFrame(),
    )


# Called by variant_scores.html?
@prioritisation_bp.route("/stream_progress/<disease_code>/<term>")
def stream_progress(disease_code: str, term: str, ancestry: str = "EUR"):
    """Stream progress updates and final results."""
    no_cache = bool(request.args.get("flush", False))
    protective = bool(request.args.get("protective", False))
    print(f"Flush: {no_cache}, Protective: {protective}")

    def process_variants(
        disease_code: str,
        term: str,
        phenotype_service: PhenotypeService,
        genotype_service: GenotypeService,
        nomaly_scores_service: NomalyScoreService,
        stats_service: StatsService,
        run_version: str,
        ancestry: str,
        message_queue: Queue,
        protective: bool = False,
        no_cache: bool = False,
    ):
        """Process variants using the thread pool."""
        try:
            stream_logger = StreamLogger(message_queue)
            data = get_top_variants(
                disease_code,
                term,
                phenotype_service,
                genotype_service,
                nomaly_scores_service,
                stats_service,
                run_version,
                ancestry,
                stream_logger,
                protective,
                no_cache,
            )

            # Send all stats data
            message_data = {
                "type": "results",
                "data": {
                    "metric1": {
                        "top_variants": data["metric1_top_variants"],
                        "top_gene_set": data["metric1_top_gene_set"],
                        "stats": {
                            "num_rp": int(data["num_rp"]),
                            "num_rn": int(data["num_rn"]),
                            "pvalue": float(data["metric1_pvalue"]),
                            "tpr": float(data["metric1_tpr"]),
                            "fpr": float(data["metric1_fpr"]),
                            "lrp": float(data["metric1_lrp"]),
                            "tp": int(data["metric1_tp"]),
                            "tn": int(data["metric1_tn"]),
                            "fp": int(data["metric1_fp"]),
                            "fn": int(data["metric1_fn"]),
                            "threshold": float(data["metric1_threshold"]),
                            "meaning": "True Positives, True Negatives, False Positives and False Negatives determined using a specific, empirically determinedNomaly Score threshold (0.022).",
                        },
                    },
                    "roc_stats_mcc": {
                        "top_variants": data["roc_stats_mcc_top_variants"],
                        "top_gene_set": data["roc_stats_mcc_top_gene_set"],
                        "stats": {
                            "num_rp": int(data["num_rp"]),
                            "num_rn": int(data["num_rn"]),
                            "pvalue": float(data["roc_stats_mcc_pvalue"]),
                            "tpr": float(data["roc_stats_mcc_tpr"]),
                            "fpr": float(data["roc_stats_mcc_fpr"]),
                            "lrp": float(data["roc_stats_mcc_lrp"]),
                            "tp": int(data["roc_stats_mcc_tp"]),
                            "tn": int(data["roc_stats_mcc_tn"]),
                            "fp": int(data["roc_stats_mcc_fp"]),
                            "fn": int(data["roc_stats_mcc_fn"]),
                            "threshold": float(data["roc_stats_mcc_threshold"]),
                            "meaning": "Matthews Correlation Coefficient (MCC) selects the threshold that maximises the chi-square statistic, minimising the contingency table's p-value.",
                        },
                    },
                    "roc_stats_yjs": {
                        "top_variants": data["roc_stats_yjs_top_variants"],
                        "top_gene_set": data["roc_stats_yjs_top_gene_set"],
                        "stats": {
                            "num_rp": int(data["num_rp"]),
                            "num_rn": int(data["num_rn"]),
                            "pvalue": float(data["roc_stats_yjs_pvalue"]),
                            "tpr": float(data["roc_stats_yjs_tpr"]),
                            "fpr": float(data["roc_stats_yjs_fpr"]),
                            "lrp": float(data["roc_stats_yjs_lrp"]),
                            "tp": int(data["roc_stats_yjs_tp"]),
                            "tn": int(data["roc_stats_yjs_tn"]),
                            "fp": int(data["roc_stats_yjs_fp"]),
                            "fn": int(data["roc_stats_yjs_fn"]),
                            "threshold": float(data["roc_stats_yjs_threshold"]),
                            "meaning": "Youden's J statistic (YJS) selects the threshold that maximises the difference between the true positive rate (TPR) and false positive rate FPR (TPR - FPR).",
                        },
                    },
                    "roc_stats_lrp": {
                        "top_variants": data["roc_stats_lrp_top_variants"],
                        "top_gene_set": data["roc_stats_lrp_top_gene_set"],
                        "stats": {
                            "num_rp": int(data["num_rp"]),
                            "num_rn": int(data["num_rn"]),
                            "pvalue": float(data["roc_stats_lrp_pvalue"]),
                            "tpr": float(data["roc_stats_lrp_tpr"]),
                            "fpr": float(data["roc_stats_lrp_fpr"]),
                            "lrp": float(data["roc_stats_lrp_lrp"]),
                            "tp": int(data["roc_stats_lrp_tp"]),
                            "tn": int(data["roc_stats_lrp_tn"]),
                            "fp": int(data["roc_stats_lrp_fp"]),
                            "fn": int(data["roc_stats_lrp_fn"]),
                            "threshold": float(data["roc_stats_lrp_threshold"]),
                            "meaning": "Likelihood Ratio for being Positive (LRP) selects a threshold that maximises the ratio of true positive rate over false positive rate (TPR/FPR).",
                        },
                    },
                },
            }

            if protective:
                try:
                    message_data["data"]["roc_stats_lrn_protective"] = {
                        "top_variants": data["roc_stats_lrn_protective_top_variants"],
                        "top_gene_set": data["roc_stats_lrn_protective_top_gene_set"],
                        "stats": {
                            "num_rp": int(data["num_rp"]),
                            "num_rn": int(data["num_rn"]),
                            "pvalue": float(data["roc_stats_lrn_protective_pvalue"]),
                            "tpr": float(data["roc_stats_lrn_protective_tpr"]),
                            "fpr": float(data["roc_stats_lrn_protective_fpr"]),
                            "lrp": float(data["roc_stats_lrn_protective_lrp"]),
                            "tp": int(data["roc_stats_lrn_protective_tp"]),
                            "tn": int(data["roc_stats_lrn_protective_tn"]),
                            "fp": int(data["roc_stats_lrn_protective_fp"]),
                            "fn": int(data["roc_stats_lrn_protective_fn"]),
                            "threshold": float(
                                data["roc_stats_lrn_protective_threshold"]
                            ),
                            "meaning": "Threshold by maximum fraction of cases over threshold to that of controls.",
                        },
                    }
                except KeyError as e:
                    logger.error(f"KeyError in protective message: {e}")
                    if not no_cache:
                        # Try again but with no cache set to true...
                        process_variants(
                            disease_code,
                            term,
                            phenotype_service,
                            genotype_service,
                            nomaly_scores_service,
                            stats_service,
                            run_version,
                            ancestry,
                            message_queue,
                            protective,
                            no_cache=True,
                        )

            message_queue.put(message_data)
            print("Results message sent")

        except Exception as e:
            print(f"Error in process_variants: {str(e)}")
            print(traceback.format_exc())
            logger.error(f"Error processing variants: {e}")
            message_queue.put({"type": "error", "data": str(e)})
        finally:
            print("Sending done message...")
            message_queue.put({"type": "done"})
            print("Done message sent")

    def generate():
        message_queue = Queue()
        services: ServiceRegistry = current_app.extensions["nomaly_services"]

        # Get run version and ancestry from the session
        run_version = session.get("run_version", "Run-v1")
        ancestry = session.get("ancestry", "EUR")
        stats_service = services.stats_registry.get(run_version, ancestry)

        try:
            # Submit to thread pool with the service
            print("Submitting to thread pool...")
            _ = variant_processor.submit(
                process_variants,
                disease_code,
                term,
                services.phenotype,
                services.genotype,
                services.nomaly_score,
                stats_service,
                run_version,
                ancestry,
                message_queue,
                protective,
                no_cache,
            )
            print("Submitted to thread pool")

            # Stream messages as they arrive
            while True:
                try:
                    print("Waiting for message...")
                    msg = message_queue.get(timeout=600)  # 10 minute timeout
                    print(f"Got message of type: {msg.get('type')}")
                    yield f"data: {json.dumps(msg)}\n\n"
                    if msg["type"] == "done":
                        print("Got done message, breaking")
                        break
                except Exception as e:
                    print(f"Error in message loop: {str(e)}")
                    logger.error(f"Error in stream: {e}")
                    yield f"data: {json.dumps({'type': 'error', 'data': str(e)})}\n\n"
                    yield f"data: {json.dumps({'type': 'done'})}\n\n"
                    break

        except Exception as e:
            logger.error(f"Error in stream: {e}")
            yield f"data: {json.dumps({'type': 'error', 'data': str(e)})}\n\n"
            yield f"data: {json.dumps({'type': 'done'})}\n\n"

    return Response(stream_with_context(generate()), mimetype="text/event-stream")


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


def main():
    """This code exists for debugging purposes.

    NOTES:

    calls:
    get_top_variants
        (phecode, term)

        1a) Gets the case and ctrl EIDs for the phecode (x ancestry)
        1b) Gets the scores for the case and ctrl EIDs (x run_version)

        2a) Gets the stats for the term
        2b) Gets the case EIDs 'above threshold' for each stat

        3) Gets the variants for the term

        calls:
        term_variant_prioritisation
            (term_variants, case_eids_above_threshold)

            1) Gets the *genotypes* for the case EIDs above the threshold
            2a) Gets the 'nomaly variant scores' for the term variants

            calls:
            process_individual_variants
                (eids, variants, genotypes)


                calls:
                individual_variant_prioritisation
                    (row, term_variant_scores)

    """
    run = "Run-v1"
    ancestry = "EUR"

    phecode = "722.7"
    term = "KW:1167"
    term = "GO:1900181"

    phecode = "290.11"
    term = "MP:0005179"

    from data_services import StatsRegistry

    phenotype_service = PhenotypeService(Config.PHENOTYPES_HDF)
    genotype_service = GenotypeService(Config.GENOTYPES_HDF)
    nomaly_score_service = NomalyScoreService(Config.NOMALY_SCORES_H5)
    stats_registry = StatsRegistry(Config.STATS_SELECTOR)

    data = get_top_variants(
        phecode,
        term,
        phenotype_service,
        genotype_service,
        nomaly_score_service,
        stats_registry.get(run, ancestry),
        run,
        ancestry,
        StreamLogger(Queue()),
        no_cache=True,
        # protective=True,
    )

    data = convert_numpy_types(data)

    print(json.dumps(data, indent=2, sort_keys=True))

    exit(0)

    #
    #
    # DOING SOMEETHING ELSE BELOW
    #
    #

    phecodes = phenotype_service.phecodes

    # HACK, but fine here
    stats_service = stats_registry.get(run, ancestry)
    terms = stats_service._hdf.terms

    for _ in range(100_000):
        phecode = np.random.choice(phecodes)
        term = np.random.choice(terms)
        run = np.random.choice(["Run-v1", "Run-v2"])
        ancestry = np.random.choice(["EUR", "AFR", "EAS", "SAS"])

        print(f"Running {run} {ancestry} for {phecode} {term}")

        try:
            print(f"Running {run} {ancestry} for {phecode} {term}")
            print(f"Running {run} {ancestry} for {phecode} {term}")
            print(f"Running {run} {ancestry} for {phecode} {term}")
            data = get_top_variants(
                phecode,
                term,
                phenotype_service,
                genotype_service,
                nomaly_score_service,
                stats_registry.get(run, ancestry),
                run,
                ancestry,
                no_cache=True,
                # protective=True,
            )
        except Exception as e:
            print(f"Error in main: {e}")
            logger.error(f"Error in main: {e}")

        print("OK")


if __name__ == "__main__":
    main()
