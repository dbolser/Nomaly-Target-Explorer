"""
Prioritisation of individual variants by Nomaly prediction scores and disease status.
per disease and term.
Genes are prioritised by the sum of Nomaly scores of their variants.
"""

import json
import logging
import pickle
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
    NomalyDataService,
    NomalyScoreService,
    PhenotypeService,
    ServiceRegistry,
    StatsRegistry,
    StatsService,
)
from db import get_term_domains, get_term_genes, get_term_names, get_term_variants

logger = logging.getLogger(__name__)


# # Create a 'dummy' profile decorator if we don't have line_profiler installed
# try:
#     from line_profiler import profile  # type: ignore
# except ImportError:

#     def profile(func):
#         return func


# TODO: Move to a new variant scores (or 'nomaly_data') data service?

# In memory data
# variant_scores_old = pd.read_csv(
#     "/data/clu/ukbb/variantscores.tsv", sep="\t", index_col="variant_id"
# )

variant_scores = pd.read_csv(
    "/data/clu/ukbb/genotype_counts_with_vs.tsv",
    sep="\t",
    index_col="nomaly_variant_id",
)

variant2gene = pd.read_csv(
    "/data/clu/ukbb/variant2gene.tsv", sep="\t", header=None, index_col=0
)

# Global thread pool for variant processing
# Limit concurrent processing to avoid overwhelming the server
MAX_WORKERS = 20  # Adjust based on server capacity
variant_processor = ThreadPoolExecutor(
    max_workers=MAX_WORKERS, thread_name_prefix="variant_processor"
)


def log_and_stream(message, stream_logger=None, level="info"):
    if stream_logger:
        stream_logger.info(message)
    else:
        getattr(logger, level)(message)


# Convert NumPy types to native Python types
def convert_numpy_types(obj):
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
    """Helper class to capture and stream progress messages."""

    def __init__(self, message_queue):
        self.message_queue = message_queue

    def info(self, message):
        """Send a progress message immediately to the queue."""
        msg = {"type": "progress", "data": message}
        self.message_queue.put(msg)


def read_cases_for_disease_code(phecode: str) -> dict:
    """Read the case information for the disease code."""

    # TODO: Switch to the phenotype service here!
    ukbb_pheno_dir = Config.UKBB_PHENO_DIR
    with open(
        f"{ukbb_pheno_dir}/phecode_cases_excludes/phecode_{phecode}.pkl", "rb"
    ) as f:
        cases = pickle.load(f)
    return cases


def read_nomaly_filtered_genotypes_new(
    eids: np.ndarray,
    vids: np.ndarray,
    genotype_service: GenotypeService,
) -> Dict[str, Any]:
    """Read genotypes using the new genotype service method."""
    genotypes = genotype_service.get_genotypes(eids=eids, vids=vids, nomaly_ids=True)

    # Transpose the genotypes matrix
    genotypes = genotypes.T

    return {
        "row_eids": eids,
        "col_variants": vids,
        "data": genotypes,
        "error_variants": [],
    }


# @profile
def individual_variant_prioritisation(row, term_variant_scores) -> pd.DataFrame:
    """Return numpy array of variant scores for the selected variants."""
    indices = term_variant_scores.index
    sel_vs = term_variant_scores.to_numpy()

    geno_matrix = np.zeros((len(row), 3))
    for i, val in enumerate(row):
        if val != -9:
            geno_matrix[i, int(val)] = 1

    scores = np.dot(geno_matrix, sel_vs.T)
    diagonal_scores = np.diagonal(scores)

    sorted_indices = np.argsort(diagonal_scores)[::-1]
    sorted_scores = diagonal_scores[sorted_indices]
    sorted_variants = indices[sorted_indices]

    top_variants = pd.DataFrame({"vs": sorted_scores}, index=sorted_variants)
    top_variants = top_variants.assign(rank=range(1, len(top_variants) + 1))
    top_variants = top_variants[(top_variants["vs"] > 1) | (top_variants["rank"] <= 5)]

    return top_variants


# @profile
def process_individual_variants(sel_genotypes, term_variant_scores):
    """Process variants for all individuals more efficiently"""
    ind_top_variants = defaultdict(int)
    ind_top_variants_scores = defaultdict(float)

    for individual in sel_genotypes["data"]:
        top_variants = individual_variant_prioritisation(
            individual, term_variant_scores
        )
        for variant, vs in zip(top_variants.index, top_variants["vs"]):
            ind_top_variants[variant] += 1
            ind_top_variants_scores[variant] += vs

    return pd.DataFrame(
        {
            "num_individuals": pd.Series(ind_top_variants),
            "vs": pd.Series(ind_top_variants_scores),
        }
    )


# @profile
def term_variant_prioritisation(
    vids: pd.DataFrame,
    eids: np.ndarray,
    ancestry: str,
    genotype_service: GenotypeService,
    stream_logger=None,
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
        eids, vids["variant_id"].to_numpy(), genotype_service
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
    # top_variants["gene"] = top_variants["gene"].map(lambda x: ", ".join(x))

    log_and_stream(f"Term variants: {len(vids)}", stream_logger)
    log_and_stream(f"Top variants: {len(top_variants)}", stream_logger)
    log_and_stream(
        f"Individual top variants: {len(ind_top_variants_df)}", stream_logger
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
    services: ServiceRegistry,
    run_version: str,
    ancestry: str,
    stream_logger=None,
    protective: bool = False,
    no_cache: bool = False,
) -> dict:
    """Get the top variants for the disease and term."""

    # Try to load from cache first (unless no_cache is True)
    if not no_cache:
        try:
            cached_results = load_cached_results(phecode, term, run_version, ancestry)
            if cached_results:
                if stream_logger:
                    stream_logger.info("Loaded results from cache")
                return cached_results
        except Exception as e:
            logger.error(f"Error loading cached results: {e}")

    # Compute results if not cached...
    log_and_stream(
        f"Computing results (not found in cache, flush={no_cache}, protective={protective})",
        stream_logger,
    )

    # Get phenotype data for the given phecode...
    phenotypes = services.phenotype.get_cases_for_phecode(phecode, ancestry)

    log_and_stream(
        f"Got {len(phenotypes)} samples for phecode {phecode} with ancestry {ancestry}",
        stream_logger,
    )

    # Doing this outside the loop makes sense, but it doesn't save much time.
    term_variant_ids = get_term_variants(term)

    log_and_stream(
        f"Got {len(term_variant_ids)} variants for term {term}", stream_logger
    )

    # cases_info = read_cases_for_disease_code(disease_code)
    # cases_eids = list(cases_info["cases"])
    case_eids = phenotypes[phenotypes["phenotype"] == 1]["eid"].to_numpy()
    ctrl_eids = phenotypes[phenotypes["phenotype"] == 0]["eid"].to_numpy()
    excl_eids = phenotypes[phenotypes["phenotype"] == 9]["eid"]

    assert len(phenotypes) == len(case_eids) + len(ctrl_eids) + len(excl_eids)

    log_and_stream(
        f"Got {len(case_eids)} cases for {phecode} with ancestry {ancestry}",
        stream_logger,
    )

    # Assert that the case eids are sorted
    assert np.all(np.diff(case_eids) > 0), "case_eids are not sorted!"

    case_scores = services.nomaly_score.get_scores_by_eids_unsorted(
        case_eids, terms=np.array([term])
    )

    assert len(case_eids) == len(case_scores)

    ctrl_scores = services.nomaly_score.get_scores_by_eids_unsorted(
        ctrl_eids, terms=np.array([term])
    )

    assert len(ctrl_eids) == len(ctrl_scores)

    # TODO: Not sure the stats service can access this funciton, we need to use
    # the registry here to get the appropriate service for hte run and
    # population

    stats_registry: StatsRegistry = services.stats_registry
    stats_service: StatsService = stats_registry.get(run_version, ancestry)

    stats = stats_service.get_term_stats(term=term, phecode=phecode)

    # Quick cleanup of stats...

    # Just ignore these!
    stats[stats.filter(like="_threshold").columns] = stats.filter(
        like="_threshold"
    ).replace(float("inf"), 0)
    stats[stats.filter(like="_or").columns] = stats.filter(like="_or").fillna(1)

    # Convert stats to dict..
    stats = stats.to_dict(orient="records")[0]

    # Ensure NaN values are converted to None
    if pd.isna(stats["metric1_pvalue"]):
        stats["metric1_pvalue"] = 1


    # TODO: Something is wrong somewhere...
    if len(case_eids) != stats["num_rp"] or len(ctrl_eids) != stats["num_rn"]:
        # TODO: DEBUG THIS!
        logger.warning(
            f"Phecode {phecode} has {len(case_eids)} / {len(ctrl_eids)} cases and controls"
            f" from phenotype service but {stats['num_rp']} / {stats['num_rn']} cases and controls"
            f" from stats service! (Excludes = {len(excl_eids)}).\n"
            "I THINK THIS IS TODO WITH EXCLUDED SEX-MATCHED CONTROLS!"
        )

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
        # FUCK... How did this happen?
        if np.any(stats[f"{stat}_threshold"] <= 0) and np.any(stats[f"{stat}_tp"] == 0):
            # I think this should have been set this way?
            stats[f"{stat}_threshold"] = 1e308

        case_eids_above_threshold = case_eids[case_scores >= stats[f"{stat}_threshold"]]

        if stat.endswith("protective"):
            case_eids_above_threshold = ctrl_eids[
                ctrl_scores >= stats[f"{stat}_threshold"]
            ]

        print(f"Selected for statistic {stat}: {len(case_eids_above_threshold)}")

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
            term_variant_ids,
            case_eids_above_threshold,
            ancestry,
            services.genotype,
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
    services = current_app.extensions["nomaly_services"]
    phenotype_service = services.phenotype

    """Show the variant scores page."""
    # Get flush parameter
    flush = bool(request.args.get("flush", False))

    # Get run version and ancestry from the session
    run_version = session.get("run_version", "Run-v1")
    ancestry = session.get("ancestry", "EUR")

    # Get protecive parameter
    protective = bool(request.args.get("protective", False))

    print(f"Flush: {flush}, Protective: {protective}")

    # Get phecode data
    data = get_phecode_data(disease_code, phenotype_service, ancestry)

    # Get term data
    term_names = get_term_names([term])
    term_domains = get_term_domains([term])
    term_genes = get_term_genes([term])

    # Add term details to data
    data["term"] = term
    data["termname"] = term_names.get(term, "")
    data["domainlen"] = len(term_domains.get(term, []))

    genes = term_genes[term_genes["term"] == term]["gene"].tolist()
    data["genelen"] = len(genes)
    data["genes"] = ", ".join(genes) if len(genes) < 5 else f"{len(genes)} genes"

    return render_template(
        "variant_scores.html",
        disease_code=disease_code,
        term=term,
        data=data,
        flush=flush,  # Pass to template
        protective=protective,
        top_variants=pd.DataFrame(),
        top_gene_set=pd.DataFrame(),
    )


# Called by variant_scores.html?
@prioritisation_bp.route("/stream_progress/<disease_code>/<term>")
def stream_progress(
    disease_code: str,
    term: str,
    population: str | None = None,
    biological_sex: str | None = None,
):
    """Stream progress updates and final results."""
    # Get no_cache from the flush parameter
    no_cache = bool(request.args.get("flush", False))

    # Get protective parameter
    protective = bool(request.args.get("protective", False))

    print(f"Flush: {no_cache}, Protective: {protective}")

    def process_variants(
        disease_code: str,
        term: str,
        services: ServiceRegistry,
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
                services,
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
                            services,
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

        try:
            # Submit to thread pool with the service
            print("Submitting to thread pool...")
            _ = variant_processor.submit(
                process_variants,
                disease_code,
                term,
                services,
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
        services = current_app.extensions["nomaly_services"]

        data = get_top_variants(
            phecode,
            term,
            services,
            run_version="Run-v1",
            ancestry="EUR",
            no_cache=True,
            protective=True,
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
