"""
Prioritisation of individual variants by Nomaly prediction scores and disease status.
per disease and term.
Genes are prioritised by the sum of Nomaly scores of their variants.
"""

import logging
import pickle
import os
import numpy as np
import pandas as pd
from flask import (
    Blueprint,
    Response,
    render_template,
    stream_with_context,
)
import json
from flask_login import login_required
from queue import Queue
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from pathlib import Path

from blueprints.nomaly import nomaly_genotype
from config import Config
from db import get_term_variants

logger = logging.getLogger(__name__)

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
    ukbb_pheno_dir = Config.UKBB_PHENO_DIR
    with open(
        f"{ukbb_pheno_dir}/phecode_cases_excludes/phecode_{phecode}.pkl", "rb"
    ) as f:
        cases = pickle.load(f)
    return cases


def read_nomaly_filtered_genotypes(sorted_eids, short_listed_variants) -> dict:
    """Read genotypes for the individuals and variants."""
    # sort the genotype eids
    sorted_indices = np.argsort(nomaly_genotype.individual)
    sorted_genotype_eids = nomaly_genotype.individual[sorted_indices]

    # search
    indices = np.searchsorted(sorted_genotype_eids, sorted_eids)
    valid_indices = indices[indices < len(sorted_genotype_eids)]
    matched_eids = sorted_genotype_eids[valid_indices]

    rows = len(matched_eids)
    columns = len(short_listed_variants)
    geno_matrix = np.empty((rows, columns))
    error_variants = []

    for vindex, variant_id in enumerate(short_listed_variants):
        genotype_result = nomaly_genotype.query_variantID_genotypes(variant_id)
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


def individual_variant_prioritisation(row, term_variant_scores):
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

    return top_variants.index


def term_variant_prioritisation(sorted_eids, variant_scores, term, stream_logger=None):
    """For each term, prioritise the variants for selected individuals."""
    term_variants = get_term_variants(term).drop(columns=["term"])

    sel_genotypes = read_nomaly_filtered_genotypes(
        sorted_eids, term_variants["variant_id"]
    )

    if stream_logger:
        stream_logger.info(
            f"Genotypes for {len(sel_genotypes['row_eids'])} individuals and {len(sel_genotypes['col_variants'])} variants are read."
        )

    if len(sel_genotypes["error_variants"]) > 0 and stream_logger:
        stream_logger.info(
            f"Failed to get genotype data for {len(sel_genotypes['error_variants'])} variants for term {term}."
        )

    term_variant_scores = variant_scores.loc[sel_genotypes["col_variants"]][
        ["VS00", "VS01", "VS11"]
    ]

    ind_top_variants = set()
    for row in sel_genotypes["data"]:
        ind_top_variants = ind_top_variants.union(
            individual_variant_prioritisation(row, term_variant_scores)
        )

    top_variants = term_variants[
        term_variants.variant_id.isin(ind_top_variants)
    ].sort_values(by="hmm_score", ascending=False)

    if stream_logger:
        stream_logger.info(f"Term variants: {len(term_variants)}")
        stream_logger.info(f"Top variants: {len(top_variants)}")
        stream_logger.info(f"Individual top variants: {len(ind_top_variants)}")

    return top_variants


def get_cache_path(disease_code: str, term: str) -> Path:
    """Get the cache file path for this disease/term combination."""
    cache_dir = Path(Config.VARIANT_SCORES_DIR)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"variant_prioritization_{disease_code}_{term}.json"


def load_cached_results(
    disease_code: str, term: str
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Load cached results if they exist."""
    cache_path = get_cache_path(disease_code, term)
    if not cache_path.exists():
        return None

    try:
        # Load the cached JSON data
        with open(cache_path) as f:
            data = json.load(f)

        # Convert back to DataFrames
        top_variants = pd.DataFrame.from_records(data["top_variants"])
        top_gene_set = pd.DataFrame.from_records(data["top_gene_set"]).set_index("gene")

        return top_variants, top_gene_set
    except Exception as e:
        logger.error(f"Error loading cached results: {e}")
        return None


def save_results_to_cache(
    disease_code: str, term: str, top_variants: pd.DataFrame, top_gene_set: pd.DataFrame
):
    """Save results to cache file."""
    try:
        cache_path = get_cache_path(disease_code, term)

        # Convert to the same format we send to the client
        data = {
            "top_variants": top_variants.to_dict(orient="records"),
            "top_gene_set": top_gene_set.reset_index().to_dict(orient="records"),
        }

        # Save as JSON
        with open(cache_path, "w") as f:
            json.dump(data, f, indent=2)  # indent for readability
    except Exception as e:
        logger.error(f"Error saving results to cache: {e}")


def get_top_variants(
    disease_code: str, term: str, stream_logger=None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Get the top variants for the disease and term."""
    # Try to load from cache first
    cached_results = load_cached_results(disease_code, term)
    if cached_results is not None:
        if stream_logger:
            stream_logger.info("Loaded results from cache")
        return cached_results

    if stream_logger:
        stream_logger.info("Computing results (not found in cache)")

    # Compute results if not cached
    cases_info = read_cases_for_disease_code(disease_code)
    cases_eids = list(cases_info["cases"])

    if stream_logger:
        stream_logger.info(f"Processing {len(cases_eids)} cases")

    sorted_eids = np.sort(cases_eids)
    top_variants = term_variant_prioritisation(
        sorted_eids, variant_scores, term, stream_logger
    )

    top_gene_set = (
        top_variants[["gene", "variant_id", "hmm_score"]]
        .drop_duplicates()
        .groupby("gene")
        .sum()
        .sort_values(by="hmm_score", ascending=False)
    )
    top_gene_set = top_gene_set.assign(variant_num=top_variants.groupby("gene").size())

    # Save to cache
    save_results_to_cache(disease_code, term, top_variants, top_gene_set)

    return top_variants, top_gene_set


prioritisation_bp = Blueprint("prioritisation", __name__)


@prioritisation_bp.route("/variant_scores/<disease_code>/<term>")
@login_required
def show_variant_scores(disease_code: str, term: str):
    """Show the variant scores page."""
    return render_template(
        "variant_scores.html",
        disease_code=disease_code,
        term=term,
        top_variants=pd.DataFrame(),
        top_gene_set=pd.DataFrame(),
    )


@prioritisation_bp.route("/stream_progress/<disease_code>/<term>")
@login_required
def stream_progress(disease_code: str, term: str):
    """Stream progress updates and final results."""

    def process_variants(disease_code, term, message_queue):
        """Process variants using the thread pool."""
        try:
            stream_logger = StreamLogger(message_queue)
            top_variants, top_gene_set = get_top_variants(
                disease_code, term, stream_logger
            )

            # Send results
            message_queue.put(
                {
                    "type": "results",
                    "data": {
                        "top_variants": top_variants.to_dict(orient="records"),
                        "top_gene_set": top_gene_set.to_dict(orient="records"),
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

        try:
            # Submit to thread pool
            future = variant_processor.submit(
                process_variants, disease_code, term, message_queue
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


if __name__ == "__main__":
    disease_code = "332"
    term = "GO:0030800"
    top_variants, top_gene_set = get_top_variants(disease_code, term)
    print(top_variants)
    print(top_gene_set)
