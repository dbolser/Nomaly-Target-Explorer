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

    # sort the genotype eids

    # sort the genotype eids (actually, they should be sorted already...)
    genotype_eids = genotype_service._hdf.individual
    assert np.all(genotype_eids[:-1] <= genotype_eids[1:])

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
        genotype_result = genotype_service._hdf.query_variantID_genotypes(variant_id)
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

    # Filter but get the largest hmm_score of duplicated variants.. OR NOT
    # TODO: DECIDE!
    # term_variants = (
    #     term_variants.sort_values(by="hmm_score", ascending=False)
    #     .drop_duplicates(subset="variant_id", keep="first")
    #     .reset_index(drop=True)
    # )

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

    ind_top_variants_df = process_individual_variants(
        sel_genotypes, term_variant_scores
    )

    # NOTE, if we keep 'duplicate' variants above, I think we lose them here!
    top_variants = term_variants.join(ind_top_variants_df, on="variant_id", how="right")
    top_variants = top_variants.join(term_variant_scores, on="variant_id")
    top_variants = top_variants.sort_values(by="hmm_score", ascending=False)

    if stream_logger:
        stream_logger.info(f"Term variants: {len(term_variants)}")
        stream_logger.info(f"Top variants: {len(top_variants)}")
        stream_logger.info(f"Individual top variants: {len(ind_top_variants_df)}")
    else:
        logger.info(f"Term variants: {len(term_variants)}")
        logger.info(f"Top variants: {len(top_variants)}")
        logger.info(f"Individual top variants: {len(ind_top_variants_df)}")

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

        # Save as JSON
        with open(cache_path, "w") as f:
            json.dump(data, f, indent=2)  # indent for readability
        # Meh
        return load_cached_results(disease_code, term)
    except Exception as e:
        logger.error(f"Error saving results to cache: {e}")


def get_top_variants(
    disease_code: str,
    term: str,
    genotype_service,
    nomaly_scores_service,
    stats_service,
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

    if stream_logger:
        stream_logger.info("Computing results (not found in cache)")

    # Compute results if not cached
    cases_info = read_cases_for_disease_code(disease_code)
    cases_eids = list(cases_info["cases"])

    sorted_cases_eids = np.sort(cases_eids)

    # TODO: FILTER BY SCORE HERE!
    # WE WANT HIGH SCORING AFFECTED INDIVIDUALS

    cases_scores_for_term = nomaly_scores_service._hdf.get_scores_by_eids_unsorted(
        sorted_cases_eids, terms=[term]
    )

    stats = stats_service._hdf.get_stats_by_term_phecode(
        term=term,
        phecode=disease_code,
        statstype=[
            "num_rp",
            "num_rn",
            "metric1_pvalue",
            "roc_stats_mcc_pvalue",
            "roc_stats_mcc_threshold",
            "roc_stats_yjs_pvalue",
            "roc_stats_yjs_threshold",
        ],
    )

    stats["total_cases"] = len(cases_scores_for_term)

    stats["metric1_threshold"] = 0.022

    # Lets recalculate! Yay!

    stats["metric1_tpr"] = (
        np.sum(cases_scores_for_term > stats["metric1_threshold"]) / stats["num_rp"]
    )
    stats["metric1_fpr"] = (
        np.sum(cases_scores_for_term <= stats["metric1_threshold"]) / stats["num_rn"]
    )

    stats["roc_stats_mcc_tpr"] = (
        np.sum(cases_scores_for_term > stats["roc_stats_mcc_threshold"])
        / stats["num_rp"]
    )
    stats["roc_stats_mcc_fpr"] = (
        np.sum(cases_scores_for_term <= stats["roc_stats_mcc_threshold"])
        / stats["num_rn"]
    )

    stats["roc_stats_yjs_tpr"] = (
        np.sum(cases_scores_for_term > stats["roc_stats_yjs_threshold"])
        / stats["num_rp"]
    )
    stats["roc_stats_yjs_fpr"] = (
        np.sum(cases_scores_for_term <= stats["roc_stats_yjs_threshold"])
        / stats["num_rn"]
    )

    # Get the specific EIDs above the threshold for Variant Prioritization

    stats["metric1_eids"] = sorted_cases_eids[
        cases_scores_for_term > stats["metric1_threshold"]
    ]
    stats["roc_stats_mcc_eids"] = sorted_cases_eids[
        cases_scores_for_term > stats["roc_stats_mcc_threshold"]
    ]
    stats["roc_stats_yjs_eids"] = sorted_cases_eids[
        cases_scores_for_term > stats["roc_stats_yjs_threshold"]
    ]

    # FUCK meeeeeeeeeeeeeeeeeeeeee
    for stat in ["metric1", "roc_stats_mcc", "roc_stats_yjs"]:
        eids = stats[f"{stat}_eids"]
        top_variants = term_variant_prioritisation(
            eids,
            variant_scores,
            term,
            genotype_service,
            stream_logger,
        )
        stats[f"{stat}_top_variants"] = top_variants.drop(
            columns=["term", "aa"]
        ).to_dict(orient="records")

        # Create gene set with properly formatted variant lists
        gene_variants = top_variants.groupby("gene")["variant_id"].agg(list)
        gene_hmm_score = top_variants.groupby("gene")["hmm_score"].sum()
        gene_total_vs = top_variants.groupby("gene")["vs"].sum()
        gene_num_individual = top_variants.groupby("gene")["num_individuals"].sum()
        gene_num_variants = top_variants.groupby("gene").size()

        top_gene_set = pd.DataFrame(
            {
                "variant_id": gene_variants.map(lambda x: ", ".join(x)),
                "hmm_score": gene_hmm_score,
                "total_vs": gene_total_vs,
                "variant_num": gene_num_variants,
                "num_individuals": gene_num_individual,
            }
        )
        top_gene_set = top_gene_set.sort_values("hmm_score", ascending=False)

        stats[f"{stat}_top_gene_set"] = (
            # Reset index and rename the index column to 'gene'
            top_gene_set.reset_index()
            .rename(columns={"index": "gene"})
            .to_dict(orient="records"),
        )

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

    stats_service = current_app.extensions["nomaly_services"].stats._hdf
    stats_data = get_stats_data(stats_service, disease_code, term)
    # Merge teh two dictionaries:
    data = {**data, **stats_data}

    data["stats"] = stats_data
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
def stream_progress(disease_code: str, term: str):
    """Stream progress updates and final results."""
    # Get no_cache from the flush parameter
    no_cache = bool(request.args.get("flush", False))

    def process_variants(
        disease_code,
        term,
        genotype_service,
        nomaly_scores_service,
        stats_service,
        message_queue,
        no_cache,
    ):
        """Process variants using the thread pool."""
        try:
            stream_logger = StreamLogger(message_queue)
            top_variants, top_gene_set = get_top_variants(
                disease_code,
                term,
                genotype_service,
                nomaly_scores_service,
                stats_service,
                stream_logger,
                no_cache,
            )

            # Format the gene set variant lists with proper comma separation
            top_gene_set_dict = top_gene_set.reset_index().to_dict(orient="records")
            for record in top_gene_set_dict:
                if "variant_id" in record and record["variant_id"]:
                    variants = str(record["variant_id"]).split(",")
                    record["variant_id"] = ", ".join(v.strip() for v in variants)

            # Send results
            message_queue.put(
                {
                    "type": "results",
                    "data": {
                        "top_variants": top_variants.to_dict(orient="records"),
                        "top_gene_set": top_gene_set_dict,
                    },
                }
            )

        except Exception as e:
            logger.error(f"Error processing variants: {e}")
            message_queue.put({"type": "error", "data": str(e)})
        finally:
            message_queue.put({"type": "done"})

    def generate():
        message_queue = Queue()
        services = current_app.extensions["nomaly_services"]

        try:
            # Submit to thread pool with the service
            _ = variant_processor.submit(
                process_variants,
                disease_code,
                term,
                services.genotype,
                services.nomaly_score,
                services.stats,
                message_queue,
                no_cache,
            )

            # Stream messages as they arrive
            while True:
                try:
                    msg = message_queue.get(timeout=60)  # 1 minute timeout
                    yield f"data: {json.dumps(msg)}\n\n"
                    if msg["type"] == "done":
                        break
                except Exception as e:
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

    disease_code = "332"
    term = "GO:0030800"

    # disease_code = "256"
    # term = "GO:0036265"  # <- Dis term booog!

    disease_code = "290.11"
    term = "GO:0016861"

    from app import create_app

    app = create_app("development")
    with app.app_context():
        services = current_app.extensions["nomaly_services"]

        data = get_top_variants(
            disease_code,
            term,
            services.genotype,
            services.nomaly_score,
            services.stats,
            no_cache=True,
        )

        # FUCK ME!
        print(data)


if __name__ == "__main__":
    main()
