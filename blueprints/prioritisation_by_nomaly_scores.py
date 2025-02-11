"""
Prioritisation of individual variants by Nomaly prediction scores and disease status.
per disease and term.
Genes are prioritised by the sum of Nomaly scores of their variants.
"""

import json
import logging
import pickle
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from queue import Queue

import numpy as np
import pandas as pd
from flask import (
    Blueprint,
    Response,
    current_app,
    render_template,
    request,
    stream_with_context,
)

from blueprints.phecode import get_phecode_data
from config import Config
from db import get_term_domains, get_term_genes, get_term_names, get_term_variants
from data_services.stats import StatsHDF5

logger = logging.getLogger(__name__)


# TODO: Move to a new variant scores (or 'nomaly_data') data service?

# In memory data
variant_scores = pd.read_csv(
    "/data/clu/ukbb/variantscores.tsv", sep="\t", index_col="variant_id"
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


def read_nomaly_filtered_genotypes(
    sorted_eids, short_listed_variants, genotype_service
) -> dict:
    """Read genotypes for the individuals and variants."""

    genotype_eids = genotype_service.individual
    assert np.all(genotype_eids[:-1] <= genotype_eids[1:])

    # sort the genotype eids (actually, they should be sorted already...)
    sorted_indices = np.argsort(genotype_eids)
    sorted_genotype_eids = genotype_eids[sorted_indices]

    # search
    indices = np.searchsorted(sorted_genotype_eids, sorted_eids)
    valid_indices = indices[indices < len(sorted_genotype_eids)]
    matched_eids = sorted_genotype_eids[valid_indices]

    rows = len(matched_eids)
    columns = len(short_listed_variants)
    geno_matrix = np.empty((rows, columns))
    error_variants = []

    for vindex, variant_id in enumerate(short_listed_variants):
        genotype_result = genotype_service.query_variantID_genotypes(variant_id)
        if genotype_result is None:
            logger.warning(f"No genotype data found for variant {variant_id}")
            error_variants.append(variant_id)
            continue

        sorted_genotype_eids, sorted_genotypes = genotype_result
        matched_genotypes = sorted_genotypes[valid_indices]
        assert (sorted_genotype_eids[valid_indices] == matched_eids).all()
        geno_matrix[:, vindex] = matched_genotypes

    return {
        "row_eids": matched_eids,
        "col_variants": short_listed_variants,
        "data": geno_matrix,
        "error_variants": error_variants,
    }


def individual_variant_prioritisation(row, term_variant_scores) -> pd.DataFrame:
    """Return numpy array of variant scores for the selected variants."""
    indices = term_variant_scores.index
    sel_vs = term_variant_scores.to_numpy()

    geno_matrix = np.zeros((len(row), 3))
    for i, val in enumerate(row):
        if val != -1:
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


def term_variant_prioritisation(
    sorted_eids, variant_scores, term, genotype_service, stream_logger=None
):
    """For each term, prioritise the variants for selected individuals."""
    term_variants = get_term_variants(term)

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

    if stream_logger:
        stream_logger.info(f"Reading genotypes for {len(term_variants)} variants")
    else:
        logger.info(f"Reading genotypes for {len(term_variants)} variants")

    sel_genotypes = read_nomaly_filtered_genotypes(
        sorted_eids, term_variants["variant_id"], genotype_service
    )

    if stream_logger:
        stream_logger.info(
            f"Genotypes for {len(sel_genotypes['row_eids'])} individuals and {len(sel_genotypes['col_variants'])} variants are read."
        )
        if len(sel_genotypes["error_variants"]) > 0:
            stream_logger.warning(
                f"Failed to get genotype data for {len(sel_genotypes['error_variants'])} variants for term {term}."
            )
    else:
        logger.info(
            f"Genotypes for {len(sel_genotypes['row_eids'])} individuals and {len(sel_genotypes['col_variants'])} variants are read."
        )
        if len(sel_genotypes["error_variants"]) > 0:
            logger.warning(
                f"Failed to get genotype data for {len(sel_genotypes['error_variants'])} variants for term {term}."
            )

    term_variant_scores = variant_scores.loc[sel_genotypes["col_variants"]][
        ["VS00", "VS01", "VS11"]
    ]

    print(
        f"We read {len(term_variant_scores)} variant scores for our {len(term_variants)} variants"
    )

    ind_top_variants_df = process_individual_variants(
        sel_genotypes, term_variant_scores
    )

    print(
        f"For some reason, we get just {len(ind_top_variants_df)} top variants for {len(sel_genotypes['row_eids'])} individuals"
    )

    top_variants = term_variants.join(ind_top_variants_df, on="variant_id", how="right")
    top_variants = top_variants.join(term_variant_scores, on="variant_id")

    if stream_logger:
        stream_logger.info(f"Term variants: {len(term_variants)}")
        stream_logger.info(f"Top variants: {len(top_variants)}")
        stream_logger.info(f"Individual top variants: {len(ind_top_variants_df)}")
    else:
        logger.info(f"Term variants: {len(term_variants)}")
        logger.info(f"Individual top variants: {len(ind_top_variants_df)}")
        logger.info(f"Top variants: {len(top_variants)}")

    return top_variants


def get_cache_path(disease_code: str, term: str) -> Path:
    """Get the cache file path for this disease/term combination."""
    cache_dir = Path(Config.VARIANT_SCORES_DIR)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"variant_prioritization_{disease_code}_{term}.json"


def load_cached_results(disease_code: str, term: str) -> dict | None:
    """Load cached results if they exist."""
    cache_path = get_cache_path(disease_code, term)
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


def save_results_to_cache(disease_code: str, term: str, data: dict):
    """Save results to cache file."""
    try:
        cache_path = get_cache_path(disease_code, term)

        # Convert data before saving
        converted_data = convert_numpy_types(data)

        # Save as JSON
        with open(cache_path, "w") as f:
            json.dump(converted_data, f, indent=2)  # indent for readability
        # Meh
        # return load_cached_results(disease_code, term)
    except Exception as e:
        logger.error(f"Error saving results to cache: {e}")


def get_top_variants(
    disease_code: str,
    term: str,
    phenotype_service,
    genotype_service,
    nomaly_scores_service,
    stats_service,
    population: str | None = None,
    biological_sex: str | None = None,
    stream_logger=None,
    no_cache: bool = False,
) -> dict:
    """Get the top variants for the disease and term."""

    # Try to load from cache first (unless no_cache is True)
    if not no_cache:
        cached_results = load_cached_results(disease_code, term)
        if cached_results:
            if stream_logger:
                stream_logger.info("Loaded results from cache")
            return cached_results

    # Compute results if not cached...

    if stream_logger:
        stream_logger.info("Computing results (not found in cache)")

    # Get phenotype data for the given phecode...
    eids, phenotypes = phenotype_service.get_cases_for_phecode(
        disease_code, population, biological_sex
    )

    print(
        f"Got {len(eids)} cases for phecode {disease_code} with population {population} and sex {biological_sex}"
    )

    # cases_info = read_cases_for_disease_code(disease_code)
    # cases_eids = list(cases_info["cases"])
    case_eids = eids[phenotypes == 1]
    control_eids = eids[phenotypes == 0]
    exclude_eids = eids[phenotypes == 9]

    assert len(eids) == len(case_eids) + len(control_eids) + len(exclude_eids)

    print(
        f"Got {len(case_eids)} cases for {disease_code} with population {population} and sex {biological_sex}"
    )

    case_eids_sorted = np.sort(case_eids)

    case_scores = nomaly_scores_service.get_scores_by_eids_unsorted(
        case_eids_sorted, terms=[term]
    )

    assert len(case_eids_sorted) == len(case_scores)

    stats = stats_service.get_stats_by_term_phecode(
        term=term,
        phecode=disease_code,
        statstype=[
            "num_rp",
            "num_rn",
            #
            "metric1_pvalue",
            #
            "roc_stats_mcc_pvalue",
            "roc_stats_mcc_or",
            "roc_stats_mcc_threshold",
            "roc_stats_mcc_tp",
            "roc_stats_mcc_fp",
            "roc_stats_mcc_fn",
            "roc_stats_mcc_tn",
            #
            "roc_stats_yjs_pvalue",
            "roc_stats_yjs_or",
            "roc_stats_yjs_threshold",
            "roc_stats_yjs_tp",
            "roc_stats_yjs_fp",
            "roc_stats_yjs_fn",
            "roc_stats_yjs_tn",
            #
            "roc_stats_lrp_pvalue",
            "roc_stats_lrp_or",
            "roc_stats_lrp_threshold",
            "roc_stats_lrp_tp",
            "roc_stats_lrp_fp",
            "roc_stats_lrp_fn",
            "roc_stats_lrp_tn",
            #
            "roc_stats_lrn_protective_pvalue",
            "roc_stats_lrn_protective_or",
            "roc_stats_lrn_protective_threshold",
            "roc_stats_lrn_protective_tp",
            "roc_stats_lrn_protective_fp",
            "roc_stats_lrn_protective_fn",
            "roc_stats_lrn_protective_tn",
        ],
    )

    # Something is wrong...

    if len(case_eids) != stats["num_rp"] or len(control_eids) != stats["num_rn"]:
        # TODO: DEBUG THIS!
        logger.warning(
            f"Phecode {disease_code} has {len(case_eids)} / {len(control_eids)} cases and controls"
            + f" from phenotype service but {stats['num_rp']} / {stats['num_rn']} cases and controls"
            + f" from stats service! (Excludes = {len(exclude_eids)})"
        )

    # Fill in the 'missing' values for metric1

    stats["metric1_threshold"] = 0.022

    stats["metric1_tp"] = np.sum(case_scores >= stats["metric1_threshold"])
    stats["metric1_fp"] = np.sum(case_scores < stats["metric1_threshold"])
    stats["metric1_fn"] = len(case_eids) - stats["metric1_tp"]
    stats["metric1_tn"] = len(control_eids) - stats["metric1_fp"]

    for stat in [
        "metric1",
        "roc_stats_mcc",
        "roc_stats_yjs",
        "roc_stats_lrp",
        # "roc_stats_lrn_protective",
    ]:
        selected_eids = case_eids_sorted[case_scores >= stats[f"{stat}_threshold"]]

        print(f"Selected for statistic {stat}: {len(selected_eids)}")

        if len(selected_eids) != stats[f"{stat}_tp"]:
            logger.warning(
                f"Phecode {disease_code} has {len(selected_eids)} True Positives for {stat} from "
                + f"nomaly scores but {stats[f'{stat}_tp']} True Positives from stats service!"
            )

        top_variants = term_variant_prioritisation(
            selected_eids,
            variant_scores,
            term,
            genotype_service,
            stream_logger,
        )

        # Calculate True Positive Rate (Sensitivity (TP/(TP+FN) or TP/RP))
        # Calculate False Positive Rate (1 - Specificity (FP/(FP+TN) or FP/RN))
        stats[f"{stat}_tpr"] = stats[f"{stat}_tp"] / stats["num_rp"]
        stats[f"{stat}_fpr"] = stats[f"{stat}_fp"] / stats["num_rn"]

        # Calculate Positive Likelihood Ratio (TPR/FPR) (TPR/(1-TNR))
        stats[f"{stat}_lrp"] = stats[f"{stat}_tpr"] / (stats[f"{stat}_fpr"] + 1e-10)

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

    # Ensure NaN values are converted to None
    if pd.isna(stats["metric1_pvalue"]):
        stats["metric1_pvalue"] = 1

    # Save to cache
    save_results_to_cache(disease_code, term, stats)

    return stats


prioritisation_bp = Blueprint("prioritisation", __name__)


@prioritisation_bp.route("/variant_scores/<disease_code>/<term>")
def show_variant_scores(disease_code: str, term: str):
    """Show the variant scores page."""
    # Get flush parameter
    flush = bool(request.args.get("flush", False))

    # Get phecode data
    data = get_phecode_data(disease_code)

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

    # Initialize empty stats values
    data.update(
        {
            "metric1_pvalue": 1,
            "metric1_tpr": 0,
            "metric1_fpr": 0,
            "metric1_lrp": 0,
            "metric1_tp": 0,
        }
    )

    return render_template(
        "variant_scores.html",
        disease_code=disease_code,
        term=term,
        data=data,
        flush=flush,  # Pass to template
        top_variants=pd.DataFrame(),
        top_gene_set=pd.DataFrame(),
    )


def get_stats_data(stats_service: StatsHDF5, disease_code: str, term: str) -> dict:
    """Get the stats data for the disease and term."""
    stats_data = stats_service.get_stats_by_term_phecode(term, disease_code)
    return stats_data


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

    def process_variants(
        disease_code,
        term,
        phenotype_service,
        genotype_service,
        nomaly_scores_service,
        stats_service,
        population,
        biological_sex,
        message_queue,
        no_cache,
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
                population,
                biological_sex,
                stream_logger,
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
                            "pvalue": float(data["metric1_pvalue"]),
                            "tpr": float(data["metric1_tpr"]),
                            "fpr": float(data["metric1_fpr"]),
                            "lrp": float(data["metric1_lrp"]),
                            "tp": int(data["metric1_tp"]),
                            "threshold": float(data["metric1_threshold"]),
                            "meaning": "True Positives, True Negatives, False Positives and False Negatives determined using a specific nomaly score threshold",
                        },
                    },
                    "roc_stats_mcc": {
                        "top_variants": data["roc_stats_mcc_top_variants"],
                        "top_gene_set": data["roc_stats_mcc_top_gene_set"],
                        "stats": {
                            "pvalue": float(data["roc_stats_mcc_pvalue"]),
                            "tpr": float(data["roc_stats_mcc_tpr"]),
                            "fpr": float(data["roc_stats_mcc_fpr"]),
                            "lrp": float(data["roc_stats_mcc_lrp"]),
                            "tp": int(data["roc_stats_mcc_tp"]),
                            "threshold": float(data["roc_stats_mcc_threshold"]),
                            "meaning": "Threshold determined by maximum 'best-separation' p-value (MCC). Note this is a two-sided P value. To see if people above the threshold have a predisposition to be affected or not-affected, check the TPR and FPR below.",
                        },
                    },
                    "roc_stats_yjs": {
                        "top_variants": data["roc_stats_yjs_top_variants"],
                        "top_gene_set": data["roc_stats_yjs_top_gene_set"],
                        "stats": {
                            "pvalue": float(data["roc_stats_yjs_pvalue"]),
                            "tpr": float(data["roc_stats_yjs_tpr"]),
                            "fpr": float(data["roc_stats_yjs_fpr"]),
                            "lrp": float(data["roc_stats_yjs_lrp"]),
                            "tp": int(data["roc_stats_yjs_tp"]),
                            "threshold": float(data["roc_stats_yjs_threshold"]),
                            "meaning": "Threshold by maximum fraction of cases over threshold to that of controls.",
                        },
                    },
                    "roc_stats_lrp": {
                        "top_variants": data["roc_stats_lrp_top_variants"],
                        "top_gene_set": data["roc_stats_lrp_top_gene_set"],
                        "stats": {
                            "pvalue": float(data["roc_stats_lrp_pvalue"]),
                            "tpr": float(data["roc_stats_lrp_tpr"]),
                            "fpr": float(data["roc_stats_lrp_fpr"]),
                            "lrp": float(data["roc_stats_lrp_lrp"]),
                            "tp": int(data["roc_stats_lrp_tp"]),
                            "threshold": float(data["roc_stats_lrp_threshold"]),
                            "meaning": "Threshold by maximum fraction of cases over threshold to that of controls.",
                        },
                    },
                    # "roc_stats_lrn_protective": {
                    #     "top_variants": data["roc_stats_lrn_protective_top_variants"],
                    #     "top_gene_set": data["roc_stats_lrn_protective_top_gene_set"],
                    #     "stats": {
                    #         "pvalue": float(data["roc_stats_lrn_protective_pvalue"]),
                    #         "tpr": float(data["roc_stats_lrn_protective_tpr"]),
                    #         "fpr": float(data["roc_stats_lrn_protective_fpr"]),
                    #         "lrp": float(data["roc_stats_lrn_protective_lrp"]),
                    #         "tp": int(data["roc_stats_lrn_protective_tp"]),
                    #         "threshold": float(
                    #             data["roc_stats_lrn_protective_threshold"]
                    #         ),
                    #         "meaning": "Threshold by maximum fraction of cases over threshold to that of controls.",
                    #     },
                    # },
                },
            }

            print(json.dumps(message_data, indent=2, sort_keys=True))

            # Add debug print to see the stats values
            print("Stats values being sent:")
            print(f"metric1_pvalue: {data['metric1_pvalue']}")
            print(f"metric1_tpr: {data['metric1_tpr']}")
            print(f"metric1_fpr: {data['metric1_fpr']}")
            print(f"metric1_lrp: {data['metric1_lrp']}")
            print(f"metric1_tp: {data['metric1_tp']}")

            message_queue.put(message_data)
            print("Results message sent")

        except Exception as e:
            print(f"Error in process_variants: {str(e)}")
            logger.error(f"Error processing variants: {e}")
            message_queue.put({"type": "error", "data": str(e)})
        finally:
            print("Sending done message...")
            message_queue.put({"type": "done"})
            print("Done message sent")

    def generate():
        message_queue = Queue()
        services = current_app.extensions["nomaly_services"]

        try:
            # Submit to thread pool with the service
            print("Submitting to thread pool...")
            _ = variant_processor.submit(
                process_variants,
                disease_code,
                term,
                services.phenotype._hdf,
                services.genotype._hdf,
                services.nomaly_score._hdf,
                services.stats._hdf,
                population,
                biological_sex,
                message_queue,
                no_cache,
            )
            print("Submitted to thread pool")

            # Stream messages as they arrive
            while True:
                try:
                    print("Waiting for message...")
                    msg = message_queue.get(timeout=60)  # 1 minute timeout
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

    disease_code = "290.11"
    term = "CD:MESH:D009139"

    # disease_code = "334.2"
    # term = "GO:0044559"

    # disease_code = "334.2"
    # term = "GO:0044559"

    # disease_code = "324.1"
    # term = "GO:0009225"

    disease_code = "332"
    term = "GO:0030800"

    from app import create_app

    app = create_app("development")
    with app.app_context():
        services = current_app.extensions["nomaly_services"]

        data = get_top_variants(
            disease_code,
            term,
            services.phenotype._hdf,
            services.genotype._hdf,
            services.nomaly_score._hdf,
            services.stats._hdf,
            no_cache=True,
        )

        # Convert NumPy types to native Python types
        data = convert_numpy_types(data)

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

        data = convert_numpy_types(data)

        for key, value in data.items():
            if "top_variants" in key or "top_gene_set" in key:
                data[key] = len(data[key])

        print(json.dumps(data, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
