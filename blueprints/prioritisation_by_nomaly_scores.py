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

# Create a 'dummy' profile decorator if we don't have line_profiler installed
try:
    from line_profiler import profile
except ImportError:

    def profile(func):
        return func


from blueprints.phecode import get_phecode_data
from config import Config
from db import get_term_domains, get_term_genes, get_term_names, get_term_variants

logger = logging.getLogger(__name__)


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


def read_nomaly_filtered_genotypes_new(eids, vids, genotype_service) -> dict:
    genotypes = genotype_service.get_genotypes(eids=eids, vids=vids, nomaly_ids=True)

    # FUCK ME!
    genotypes = genotypes.T

    return {
        "row_eids": eids,
        "col_variants": vids,
        "data": genotypes,
        "error_variants": [],
    }


def read_nomaly_filtered_genotypes(eids, vids, genotype_service) -> dict:
    """Read genotypes for the individuals and variants.

    TODO:
        1) All EIDs are now already sorted, so we can remove the 'alignment'
           step.
        2) Instead of the messy id mapping, we can use the new
           nomaly_variant_id index.
        3) Essentially this entire function should be replaced with a
           single call to the genotype service.
    """

    genotype_eids = genotype_service.individual

    # Because we need to 'align' the eids with genotype eids, we need to jump
    # through a few hoops here ...

    # 1) Get the sort order of the genotype eids
    sorted_indices = np.argsort(genotype_eids)

    # 2) Remove eids that are not found in the genotype data
    eids = eids[np.isin(eids, genotype_eids)]

    # 3) Use searchsorted (fast) to find the indices of the eids in the
    #    sorted genotype eids
    eidx = np.searchsorted(genotype_eids[sorted_indices], eids)

    geno_matrix = np.empty((len(eids), len(vids)))
    failed_variants = []

    for v, nomaly_variant_id in enumerate(vids):
        genotype_result = genotype_service.query_variantID_genotypes(nomaly_variant_id)
        if genotype_result is None:
            logger.warning(f"No genotype data found for variant {nomaly_variant_id}")
            failed_variants.append(nomaly_variant_id)
            continue

        sorted_genotype_eids, sorted_genotypes = genotype_result

        # HOrrible code!
        assert np.all(sorted_genotype_eids[eidx] == eids)

        geno_matrix[:, v] = sorted_genotypes[eidx]

    return {
        "row_eids": genotype_eids[eidx],
        "col_variants": vids,
        "data": geno_matrix,
        "error_variants": failed_variants,
    }


@profile
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


@profile
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


@profile
def term_variant_prioritisation(vids, eids, genotype_service, stream_logger=None):
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

    if stream_logger:
        stream_logger.info(f"Reading genotypes for {len(vids)} variants")
    else:
        logger.info(f"Reading genotypes for {len(vids)} variants")

    # sel_genotypes2 = read_nomaly_filtered_genotypes_new(
    #     eids, vids["variant_id"], genotype_service
    # )

    # The CRITICAL difference betwween the new version is that I dont' flip
    # alleles. This means that the genotypes and the genotype frequencies used
    # to calculate vs are consistent. Everything is consistent with the genotype
    # ref/alt definition EXCEPT nomaly variant IDs... (conveniently, the latter
    # are used in the displya...)
    sel_genotypes = read_nomaly_filtered_genotypes_new(
        eids, vids["variant_id"], genotype_service
    )

    if stream_logger:
        stream_logger.info(
            f"Genotypes for {len(sel_genotypes['row_eids'])} individuals and {len(sel_genotypes['col_variants'])} variants are read."
        )
        if len(sel_genotypes["error_variants"]) > 0:
            stream_logger.warning(
                f"Failed to get genotype data for {len(sel_genotypes['error_variants'])} variants."
            )
    else:
        logger.info(
            f"Genotypes for {len(sel_genotypes['row_eids'])} individuals and {len(sel_genotypes['col_variants'])} variants are read."
        )
        if len(sel_genotypes["error_variants"]) > 0:
            logger.warning(
                f"Failed to get genotype data for {len(sel_genotypes['error_variants'])} variants."
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

    if stream_logger:
        stream_logger.info(f"Term variants: {len(vids)}")
        stream_logger.info(f"Top variants: {len(top_variants)}")
        stream_logger.info(f"Individual top variants: {len(ind_top_variants_df)}")
    else:
        logger.info(f"Term variants: {len(vids)}")
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
        return load_cached_results(disease_code, term)
    except Exception as e:
        logger.error(f"Error saving results to cache: {e}")


@profile
def get_top_variants(
    phecode: str,
    term: str,
    phenotype_service,
    genotype_service,
    nomaly_scores_service,
    stats_service,
    population: str | None = None,
    biological_sex: str | None = None,
    stream_logger=None,
    protective: bool = False,
    no_cache: bool = False,
) -> dict:
    """Get the top variants for the disease and term."""

    # Try to load from cache first (unless no_cache is True)
    if not no_cache:
        cached_results = load_cached_results(phecode, term)
        if cached_results:
            if stream_logger:
                stream_logger.info("Loaded results from cache")
            return cached_results

    # Compute results if not cached...

    if stream_logger:
        stream_logger.info(
            f"Computing results (not found in cache, flush={no_cache}, protective={protective})"
        )

    # Get phenotype data for the given phecode...
    eids, phenotypes = phenotype_service.get_cases_for_phecode(
        phecode, population, biological_sex
    )

    if stream_logger:
        stream_logger.info(
            f"Got {len(eids)} samples for phecode {phecode} with population {population} and sex {biological_sex}"
        )
    else:
        logger.info(
            f"Got {len(eids)} samples for phecode {phecode} with population {population} and sex {biological_sex}"
        )

    # Doing this outside the loop makes sense, but it doesn't save much time.
    term_variant_ids = get_term_variants(term)

    if stream_logger:
        stream_logger.info(f"Got {len(term_variant_ids)} variants for term {term}")
    else:
        logger.info(f"Got {len(term_variant_ids)} variants for term {term}")

    # cases_info = read_cases_for_disease_code(disease_code)
    # cases_eids = list(cases_info["cases"])
    case_eids = eids[phenotypes == 1]
    control_eids = eids[phenotypes == 0]
    exclude_eids = eids[phenotypes == 9]

    assert len(eids) == len(case_eids) + len(control_eids) + len(exclude_eids)

    if stream_logger:
        stream_logger.info(
            f"Got {len(case_eids)} cases for {phecode} with population {population} and sex {biological_sex}"
        )
    else:
        logger.info(
            f"Got {len(case_eids)} cases for {phecode} with population {population} and sex {biological_sex}"
        )

    # Assert that the case eids are sorted
    assert np.all(np.diff(case_eids) > 0), "case_eids are not sorted!"

    case_scores = nomaly_scores_service.get_scores_by_eids_unsorted(
        case_eids, terms=[term]
    )

    assert len(case_eids) == len(case_scores)

    control_scores = nomaly_scores_service.get_scores_by_eids_unsorted(
        control_eids, terms=[term]
    )

    assert len(control_eids) == len(control_scores)

    stats = stats_service.get_stats_by_term_phecode(
        term=term,
        phecode=phecode,
        statstype=[
            "num_rp",
            "num_rn",
            #
            "metric1_pvalue",
            #
            "roc_stats_mcc_pvalue",
            "roc_stats_mcc_or",
            "roc_stats_mcc_threshold",
            "roc_stats_mcc_all_index",
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

    # Quick cleanup of stats...

    # Ensure NaN values are converted to None
    if pd.isna(stats["metric1_pvalue"]):
        stats["metric1_pvalue"] = 1

    for key in stats.keys():
        if key.endswith("_threshold"):
            if stats[key] == float("inf"):
                # Just ignore these!
                stats[key] = 0
        if key.endswith("_or"):
            if np.isnan(stats[key]):
                # This is probably the wrong thing to do...
                stats[key] = 1

    # TODO: Something is wrong somewhere...
    if len(case_eids) != stats["num_rp"] or len(control_eids) != stats["num_rn"]:
        # TODO: DEBUG THIS!
        logger.warning(
            f"Phecode {phecode} has {len(case_eids)} / {len(control_eids)} cases and controls"
            + f" from phenotype service but {stats['num_rp']} / {stats['num_rn']} cases and controls"
            + f" from stats service! (Excludes = {len(exclude_eids)})."
        )

    # Fill in the 'missing' values for metric1

    stats["metric1_threshold"] = 0.022

    stats["metric1_tp"] = np.sum(case_scores >= stats["metric1_threshold"])
    stats["metric1_fp"] = np.sum(control_scores < stats["metric1_threshold"])
    stats["metric1_fn"] = len(case_eids) - stats["metric1_tp"]
    stats["metric1_tn"] = len(control_eids) - stats["metric1_fp"]

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
        if stats[f"{stat}_threshold"] <= 0 and stats[f"{stat}_tp"] == 0:
            # I think this should have been set this way?
            stats[f"{stat}_threshold"] = 1e308

        case_eids_above_threshold = case_eids[case_scores >= stats[f"{stat}_threshold"]]

        if stat.endswith("protective"):
            case_eids_above_threshold = control_eids[
                control_scores >= stats[f"{stat}_threshold"]
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
    save_results_to_cache(phecode, term, stats)

    return stats


prioritisation_bp = Blueprint("prioritisation", __name__)


@prioritisation_bp.route("/variant_scores/<disease_code>/<term>")
def show_variant_scores(disease_code: str, term: str):
    """Show the variant scores page."""
    # Get flush parameter
    flush = bool(request.args.get("flush", False))

    # Get protecive parameter
    protective = bool(request.args.get("protective", False))

    print(f"Flush: {flush}, Protective: {protective}")

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
        disease_code,
        term,
        phenotype_service,
        genotype_service,
        nomaly_scores_service,
        stats_service,
        population,
        biological_sex,
        message_queue,
        protective,
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
                            population,
                            biological_sex,
                            message_queue,
                            protective,
                            no_cache=True,
                        )

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

    app = create_app("development")
    with app.app_context():
        services = current_app.extensions["nomaly_services"]

        data = get_top_variants(
            phecode,
            term,
            services.phenotype._hdf,
            services.genotype._hdf,
            services.nomaly_score._hdf,
            services.stats._hdf,
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
