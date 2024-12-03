from flask import Blueprint, render_template, jsonify
import threading
import pandas as pd
import os
import traceback

from blueprints.phewas import phecode_level_assoc, PHEWAS_PHENO_DIR
from db import get_term_genes, get_phecode_info

import re

# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #
# GWAS_PHENO_DIR = '/data/clu/ukbb/by_pheno/'

# ----------------------------------------------------- #
# Variant Blueprint
# ----------------------------------------------------- #

# Create the blueprint
variant_bp = Blueprint("variant", __name__, template_folder="../templates")


# Route to render the Variant page
@variant_bp.route("/variant/<string:variant>", methods=["POST", "GET"])
def show_variant(variant):
    try:
        # Parse variant information regardless of format
        if ":" in variant:
            # Format: "5:33951588_C/G"
            chrom = variant.split(":")[0].replace("chr", "")
            rest = variant.split(":")[1]
            pos = rest.split("_")[0]
            alleles = rest.split("_")[1].replace("/", "_")
        else:
            # Format: "5_33951588_C_G"
            parts = variant.split("_")
            chrom = parts[0].replace("chr", "")
            pos = parts[1]
            alleles = f"{parts[2]}_{parts[3]}"  # Reconstruct alleles

        # Normalize the variant format
        normalized_variant = f"{chrom}_{pos}_{alleles}"

        # Validate the normalized format
        format_str = r"(chr)?\d+_\d+_[ACGT]+_[ACGT]+"
        if not re.match(format_str, normalized_variant):
            return render_template(
                "error.html", error="Invalid variant ID format: " + variant
            )

        # Prepare data for the template
        variant_data = {
            "variant_id": f"{chrom}:{pos}_{alleles.replace('_', '/')}",  # Display format
            "chromosome": chrom,
            "position": pos,
            "alleles": alleles.replace("_", "/"),  # Convert back to C/G for display
        }

        return render_template("variant.html", data=variant_data)

    except Exception as e:
        error_msg = f"Error processing variant {variant}: {str(e)}"
        print(error_msg)
        print(traceback.format_exc())
        return render_template("error.html", error=error_msg)

    # ----------------------------------------------------- #
    # Run pheWAS on the variant
    # ----------------------------------------------------- #
    # GenotypeHDF5, ICD10HDF5 and phecodeHDF5 are needed


def add_gene_info_to_DataTable(plot_df, variant):
    # get term to gene mapping
    print("getting term to gene mapping", flush=True)
    term_gene_df = get_term_genes(plot_df["term"].tolist())
    print("got term to gene mapping", flush=True)

    # filter gene by assoc var
    var_assoc_sig = read_gwas(phecode)
    genefilter = set([x["Gene"] for x in var_assoc_sig])
    term_gene_df_sig = term_gene_df[term_gene_df["gene"].isin(genefilter)].rename(
        columns={"gene": "sig gene"}
    )
    # term_gene_df_other = term_gene_df[~term_gene_df['gene'].isin(genefilter)]

    # group by term, no significance filter (to
    term_gene_df = (
        term_gene_df.groupby("term")["gene"]
        .apply(lambda x: ", ".join(x) if len(x) < 5 else f"{len(x)} genes")
        .reset_index()
    )

    # group by term, use sig filter, uncomment above)
    term_gene_df_sig = (
        term_gene_df_sig.groupby("term")["sig gene"]
        .apply(lambda x: ", ".join(x) if len(x) < 50 else f"{len(x)} genes")
        .reset_index()
    )

    # # fill term_gene_df_sig NA with term_gene_df_other
    # term_gene_df_other = term_gene_df.groupby('term')['gene'].apply(
    #     lambda x: ', '.join(x)
    #     ).reset_index().set_index('term')

    # term_gene_df_sig['gene'] = term_gene_df_sig['gene'].fillna(
    #     term_gene_df_sig['term'].map(lambda x: f"None ({term_gene_df_other.loc[x, 'gene']})")
    # )

    plot_df = plot_df.merge(term_gene_df, on="term", how="left")
    plot_df = plot_df.merge(term_gene_df_sig, on="term", how="left")

    # print(plot_df, flush=True)

    return plot_df


# ----------------------------------------------------- #
# Background task to run PheWAS
# ----------------------------------------------------- #

# Dictionary to store background task results
phewas_results = {}


# Background task function
def background_task(variant):
    # ----------------------------------------------------- #
    # Task is to run PheWAS and read the results
    # ----------------------------------------------------- #
    try:
        run_phewas_if_not_done(variant)
    except Exception:
        phewas_results["result"] = (
            f"Failed to get phecode-level stats for Variant {variant}, exception was <br> {traceback.format_exc()}"
        )


# Background task function for PheWAS
def run_phewas_if_not_done(variant):

    # ----------------------------------------------------- #
    # PheWAS result file path
    # ----------------------------------------------------- #
    output_prefix = f"variant_{variant}"
    phewas_path = f"{PHEWAS_PHENO_DIR}{output_prefix}.assoc_nomaly.tsv"

    # ----------------------------------------------------- #
    # Check if PheWAS has been run for this variant
    # ----------------------------------------------------- #
    if not os.path.exists(phewas_path):
        phewas_results[variant] = "No PheWAS data found for this variant. Processing..."
        
        print(variant)
        variant = variant.replace("_", ":")
        print(variant)

        assoc = phecode_level_assoc(variant)
    else:
        assoc = pd.read_csv(phewas_path, sep="\t")

    # ----------------------------------------------------- #
    # If assoc is empty, return error message
    # ----------------------------------------------------- #
    if assoc is None:

        assoc_sig = assoc[assoc["P"] < 0.05]
        result = (
            f"PheWAS identified {assoc_sig.shape[0]} phecodes has association p<0.05."
        )

    # Store the result in the task_results dictionary
    phewas_results[variant] = result


# Endpoint to trigger the PheWAS task
@variant_bp.route("/run-phewas/<string:variant>", methods=["POST"])
def run_phewas(variant):
    # Start the background task using threading
    task_thread = threading.Thread(target=run_phewas_if_not_done, args=(variant,))
    task_thread.start()
    return jsonify({"status": "Task started"}), 202


# Endpoint to get the PheWAS results
@variant_bp.route("/phewas-result/<string:variant>", methods=["GET"])
def get_phewas_result(variant):
    result = phewas_results.get(
        variant, {"result": "Processing...", "associations": []}
    )
    return jsonify(result)
