"""
Prioritisation of individual variants by Nomaly prediction scores and disease status.
per disease and term.
Genes are prioritised by the sum of Nomaly scores of their variants.
"""

from flask import Blueprint, render_template

prioritisation_bp = Blueprint("prioritisation", __name__)


@prioritisation_bp.route("/variant_scores/<disease_code>/<term>")
def show_variant_scores(disease_code: str, term: str):
    top_variants, top_gene_set = get_top_variants(disease_code, term)
    return render_template(
        "variant_scores.html",
        disease_code=disease_code,
        term=term,
        top_variants=top_variants,
        top_gene_set=top_gene_set,
    )


# Turn line profiling on or off
if False:
    from line_profiler import profile
else:
    profile = lambda x: x  # noqa: E731


import numpy as np
import pandas as pd
import pickle

from blueprints.nomaly import nomaly_genotype

from db import get_term_variants

# TODO (DB): can db functions be separated from flask functions? I had
# "RuntimeError: Working outside of application context." when importing db
# functions in the flask app. In order to use the db functions in the flask app,
# I manually set the config in the db.py

from config import Config

import logging

logger = logging.getLogger(__name__)

# In memory data
variant_scores = pd.read_csv(
    "/data/clu/ukbb/variantscores.tsv", sep="\t", index_col="variant_id"
)
variant2gene = pd.read_csv(
    "/data/clu/ukbb/variant2gene.tsv", sep="\t", header=None, index_col=0
)


@profile
def read_cases_for_disease_code(phecode: str) -> dict:
    """
    Read the case information for the disease code.
    dict_keys(['phecode', 'Sex', 'phecode_exclude_range', 'cases', 'exclude', 'controls_num'])
    """
    ukbb_pheno_dir = Config.UKBB_PHENO_DIR
    with open(
        f"{ukbb_pheno_dir}/phecode_cases_excludes/phecode_{phecode}.pkl", "rb"
    ) as f:
        cases = pickle.load(f)
    return cases


# get genotypes for selected individuals and variants
# TODO (DB): add this function to nomaly.py? I think it can be used by others.
@profile
def read_nomaly_filtered_genotypes(sorted_eids, short_listed_variants) -> dict:
    """
    Read genotypes for the individuals and the variants.
    """
    # read genotypes for variants

    # -------------------- Find the eid indices ------------------ #
    # To quickly find the genotype data for the individuals, sort both eids lists

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

    # record the variants that failed to get genotype data
    error_variants = []

    for vindex, variant_id in enumerate(short_listed_variants):
        # Get genotype data for variants one by one.
        genotype_result = nomaly_genotype.query_variantID_genotypes(variant_id)
        if genotype_result is None:
            print(f"ERROR No genotype data found for variant {variant_id}")
            error_variants.append(variant_id)
            continue

        sorted_genotype_eids, sorted_genotypes = genotype_result

        #
        matched_genotypes = sorted_genotypes[valid_indices]
        assert (sorted_genotype_eids[valid_indices] == matched_eids).all()

        geno_matrix[:, vindex] = matched_genotypes

    variants_genotypes = {
        "row_eids": matched_eids,
        "col_variants": short_listed_variants,
        "data": geno_matrix,
        "error_variants": error_variants,
    }
    return variants_genotypes


# read variant_scores from db
@profile
def individual_variant_prioritisation(row, term_variant_scores):
    """
    return numpy array of variant scores for the selected variants by the sequence of the variants.
    """
    indices = term_variant_scores.index

    sel_vs = term_variant_scores.to_numpy()

    # change row into a 2d array where 0 to 1 0 0, 1 to 0 1 0, 2 to 0 0 1, -1 to 0 0 0
    geno_matrix = np.zeros((len(row), 3))
    for i, val in enumerate(row):
        if val != -1:
            geno_matrix[i, int(val)] = 1

    # Calculate the scores by multiplying the genotype matrix with the variant scores
    scores = np.dot(geno_matrix, sel_vs.T)
    diagonal_scores = np.diagonal(scores)

    # Sort the diagonal scores and reorder indices by the sorted order
    sorted_indices = np.argsort(diagonal_scores)[::-1]
    sorted_scores = diagonal_scores[sorted_indices]
    sorted_variants = indices[sorted_indices]

    # define conditions for top contributing variants
    # top 5 variants or score >1
    top_variants = pd.DataFrame({"vs": sorted_scores}, index=sorted_variants)
    top_variants = top_variants.assign(rank=range(1, len(top_variants) + 1))
    top_variants = top_variants[(top_variants["vs"] > 1) | (top_variants["rank"] <= 5)]

    return top_variants.index


@profile
def term_variant_prioritisation(sorted_eids, variant_scores, term):
    """
    For each term, prioritise the variants for selected individuals.
    """
    # individuals list
    # individual_IDs = get_individual_IDs(term_scores, threshold, disease_phenotypes)

    # Variants related to the term (from the database),
    term_variants = get_term_variants(term).drop(
        columns=["term"]
    )  # dataframes with columns: 'term', 'variant_id', 'gene', 'aa', 'hmm_score'

    # DEBUGGING
    # sorted_eids = sorted_eids[:10]
    # term_variants = term_variants.head(10)

    # From the genotypes file, get the genotypes for the individuals and the variants
    # Prune the variants to exclude where genotypes are the same for every individual
    sel_genotypes = read_nomaly_filtered_genotypes(
        sorted_eids, term_variants["variant_id"]
    )

    assert (sel_genotypes["row_eids"] == sorted_eids).all()
    assert (sel_genotypes["col_variants"] == term_variants["variant_id"]).all()

    logger.info(
        f"Genotypes for {len(sel_genotypes['row_eids'])} individuals and {len(sel_genotypes['col_variants'])} variants are read. "
    )

    if len(sel_genotypes["error_variants"]) > 0:
        logger.warning(
            f"Failed to get genotype data for {len(sel_genotypes['error_variants'])} variants for term {term}."
        )

    term_variant_scores = variant_scores.loc[sel_genotypes["col_variants"]][
        ["VS00", "VS01", "VS11"]
    ]
    assert len(term_variant_scores) == len(sel_genotypes["col_variants"])

    # For each individual, get the variant scores by his genotype, and rank variants

    ind_top_variants = set()
    for row in sel_genotypes["data"]:
        ind_top_variants = ind_top_variants.union(
            individual_variant_prioritisation(row, term_variant_scores)
        )

    top_variants = term_variants[
        term_variants.variant_id.isin(ind_top_variants)
    ].sort_values(by="hmm_score", ascending=False)

    logger.info(f"Term variants: {len(term_variants)}")
    logger.info(f"Top variants: {len(top_variants)}")
    logger.info(f"Ind top variants: {len(ind_top_variants)}")

    return top_variants


# ----------------- Main ----------------- #


@profile
def get_top_variants(disease_code: str, term: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Get the top variants for the disease and term.
    """
    ## individuals that we are interested in
    # in disease
    cases_info = read_cases_for_disease_code(disease_code)
    cases_eids = list(cases_info["cases"])

    logger.info(f"Cases eids: {len(cases_eids)}")

    # and above threshold

    # then as numpy array and sorted
    sorted_eids = np.sort(cases_eids)

    top_variants = term_variant_prioritisation(sorted_eids, variant_scores, term)

    top_gene_set = (
        top_variants[["gene", "variant_id", "hmm_score"]]
        .drop_duplicates()
        .groupby("gene")
        .sum()
        .sort_values(by="hmm_score", ascending=False)
    )
    top_gene_set = top_gene_set.assign(variant_num=top_variants.groupby("gene").size())

    return top_variants, top_gene_set


@profile
def main():
    """ahh..."""

    # disease_code = "571.5"
    # term = "UP:UPA00240"

    disease_code = "332"
    term = "CC:MESH:D015766"

    top_variants, top_gene_set = get_top_variants(disease_code, term)

    print(top_variants)
    print(top_gene_set)

    # # top_variants
    #            variant_id           gene           aa  hmm_score
    # 40   10_102830743_G/A        CYP17A1        R496C    8.15141
    # 43   10_102830905_A/G        CYP17A1        C442R    7.48391
    # 170   15_42081094_G/A        PLA2G4D        R333W    6.51697
    # 145   15_41844510_C/T  JMJD7-PLA2G4B        R538W    6.51697
    # 146   15_41844510_C/T        PLA2G4B        R307W    6.51697
    # ..                ...            ...          ...        ...
    # 354    2_38074958_T/C         CYP1B1        Q144R    0.01456
    # 474      7_985219_G/A         CYP2W1  A125T,A181T    0.01109
    # 275   19_40880242_C/T         CYP2A7        A166T    0.01040
    # 204   15_74752237_A/G         CYP1A2        I386V    0.00139
    # 358   20_49524089_C/T          PTGIS        R275Q    0.00000

    # # top_gene_set
    #                                                 variant_id   hmm_score  variant_num
    # gene
    # CYP17A1  10_102830743_G/A10_102830905_A/G10_102830742_C...  106.030020           39
    # CYP24A1  20_54158136_G/A20_54158136_G/A20_54158136_G/A2...   36.453299           17
    # CYP4B1   1_46817100_C/T1_46817100_C/T1_46815212_C/T1_46...   35.867411           12
    # CYP27A1  2_218814186_C/T2_218814408_C/T2_218814701_C/T2...   29.151700           11
    # CYP2B6   19_41012339_C/T19_41012339_C/T19_41012316_T/C1...   24.603921           18
    # ...                                                    ...         ...          ...
    # CYP27C1                     2_127199410_T/C2_127203450_A/G    0.595420            2
    # CYP51A1                       7_92117086_G/A7_92117086_G/A    0.511540            2
    # PNPLA2                          11_823586_C/G11_823729_C/T    0.390940            2
    # PNPLA1                                      6_36295394_G/A    0.350730            1
    # MCAT                                       22_43133308_G/C    0.164270            1

    # [72 rows x 3 columns]


if __name__ == "__main__":
    main()
