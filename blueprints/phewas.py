"""

Calculate PheWAS results for a variant using about the most cumbersome and
stupid method possible.

Two matrix multiplications is probably all we need (see
'try_simplified_verssion' for wip code).

"""

import logging
from typing import List

import numpy as np
import pandas as pd
from fisher import pvalue_npy  # type: ignore

from config import Config
from data_services import GenotypeService, PhenotypeService
from db import get_all_phecodes

# # Create a 'dummy' profile decorator if we don't have line_profiler installed
# try:
#     from line_profiler import profile  # type: ignore
# except ImportError:

#     def profile(func):
#         return func


logger = logging.getLogger(__name__)


# @profile
def format_phewas_row(row: pd.Series) -> dict:
    """
    Format a PheWAS result row for display.
    Returns a dictionary with HTML-formatted values.

    # TODO: Shouldn't this be done in the frontend?
    """
    return {
        "Phecode": row["phecode"],
        "Description": row.get("description", "Unknown"),
        "Sex": row["sex"],
        "Group": row["phecode_group"],
        "P": f"{row['p_value']:.2e}<br/>",
        "OR": f"{row['odds_ratio']:.2f}<br/>",
        "Counts": f"{row['case_num']:,}<br/>{row['ctrl_num']:,}",
        "RefAF": f"{row['case_ref_af']:.5f}<br/>{row['ctrl_ref_af']:.5f}",
        "AltAF": f"{row['case_alt_af']:.5f}<br/>{row['ctrl_alt_af']:.5f}",
        "Ref_HMOZ": f"{row['case_hom_ref_num']:,}<br/>{row['ctrl_hom_ref_num']:,}",
        "Alt_HMOZ": f"{row['case_hom_alt_num']:,}<br/>{row['ctrl_hom_alt_num']:,}",
        "HTRZ": f"{row['case_het_num']:,}<br/>{row['ctrl_het_num']:,}",
    }


def format_phewas_results(phewas_df: pd.DataFrame) -> List[dict[str, str | float]]:
    """
    Format the PheWAS results for display.
    """

    # Get the phecode data (from the database)
    phecode_data = get_all_phecodes()

    # Merge the phewas results with the phecode data
    merged_results = phewas_df.merge(phecode_data, on="phecode", how="left")

    # Sanity check
    assert len(phewas_df) == len(merged_results)

    # Handle missing phecode descriptions
    merged_results["description"] = merged_results["description"].fillna("Unknown")
    merged_results["sex"] = merged_results["sex"].fillna("Unknown")
    merged_results["phecode_group"] = merged_results["phecode_group"].fillna("Unknown")

    # Format the results, row by row, and return a list of dics, one per row
    return [format_phewas_row(row) for _, row in merged_results.iterrows()]


def run_phewas(
    variant: str,
    genotype_service: GenotypeService,
    phenotype_service: PhenotypeService,
    ancestry: str = "EUR",
) -> pd.DataFrame:
    phecodes = phenotype_service.phecodes
    phenotype_eids = phenotype_service.get_eids(ancestry=ancestry)
    phenotype_matrix = phenotype_service.get_phenotypes(ancestry=ancestry)

    # Important to pass phenotype_eids here
    genotype_matrix = genotype_service.get_genotypes(
        eids=phenotype_eids, vids=np.array([variant])
    )

    assert phenotype_matrix.shape == (phenotype_eids.shape[0], phecodes.shape[0])
    assert genotype_matrix.shape == (1, phenotype_eids.shape[0])

    # TODO: Somewhere deep down, it's implicit that the underlying EIDs are all
    #       and always the same...

    # Moving on...

    # OK, this gets a bit slower, but, compared to the original crap, it's lightning fast...

    # Get the individual phenotype masks
    ctrl_phenotype_mask = phenotype_matrix == 0
    case_phenotype_mask = phenotype_matrix == 1

    # Get the individual genotype masks
    het_genotype_mask = genotype_matrix == 1
    ref_genotype_mask = genotype_matrix == 0
    alt_genotype_mask = genotype_matrix == 2

    # Now we get all the counts in two steps...

    # 1) Get the composite mask for all combinations
    ctrl_het_mask = ctrl_phenotype_mask * het_genotype_mask.T
    ctrl_ref_mask = ctrl_phenotype_mask * ref_genotype_mask.T
    ctrl_alt_mask = ctrl_phenotype_mask * alt_genotype_mask.T
    case_het_mask = case_phenotype_mask * het_genotype_mask.T
    case_ref_mask = case_phenotype_mask * ref_genotype_mask.T
    case_alt_mask = case_phenotype_mask * alt_genotype_mask.T

    # 2) Get the counts for each combination
    ctrl_het_num = ctrl_het_mask.sum(axis=0)
    ctrl_hom_ref_num = ctrl_ref_mask.sum(axis=0)
    ctrl_hom_alt_num = ctrl_alt_mask.sum(axis=0)

    case_het_num = case_het_mask.sum(axis=0)
    case_hom_ref_num = case_ref_mask.sum(axis=0)
    case_hom_alt_num = case_alt_mask.sum(axis=0)

    # Sum alleles here for convenience
    ctrl_ref_num = ctrl_hom_ref_num * 2 + ctrl_het_num
    ctrl_alt_num = ctrl_hom_alt_num * 2 + ctrl_het_num
    case_ref_num = case_hom_ref_num * 2 + case_het_num
    case_alt_num = case_hom_alt_num * 2 + case_het_num

    # Moving on...

    # Add pseudocounts prior to calculating the odds ratios
    ctrl_ref_psu = ctrl_ref_num + 1
    ctrl_alt_psu = ctrl_alt_num + 1
    case_ref_psu = case_ref_num + 1
    case_alt_psu = case_alt_num + 1

    # Calculate allele frequencies here for convenience
    ctrl_ref_af = ctrl_ref_psu / (ctrl_ref_psu + ctrl_alt_psu)
    ctrl_alt_af = ctrl_alt_psu / (ctrl_ref_psu + ctrl_alt_psu)
    case_ref_af = case_ref_psu / (case_ref_psu + case_alt_psu)
    case_alt_af = case_alt_psu / (case_ref_psu + case_alt_psu)

    # Sanity check...
    assert np.all(ctrl_ref_af + ctrl_alt_af == 1)
    assert np.all(case_ref_af + case_alt_af == 1)

    # Calculate odds ratios

    # Sanity check...
    # odds_ratio_long = (case_alt_psu / case_ref_psu) / (ctrl_alt_psu / ctrl_ref_psu)
    # np.allclose(odds_ratio, odds_ratio_long)

    odds_ratio = (case_alt_psu * ctrl_ref_psu) / (case_ref_psu * ctrl_alt_psu)

    # Now use vecttorized fishers exact...
    # https://stackoverflow.com/questions/34947578/how-to-vectorize-fishers-exact-test

    _, _, twosided = pvalue_npy(
        case_alt_num.astype(np.uint32),
        case_ref_num.astype(np.uint32),
        ctrl_alt_num.astype(np.uint32),
        ctrl_ref_num.astype(np.uint32),
    )

    results_df = pd.DataFrame(
        {
            "phecode": phecodes,
            "case_num": case_phenotype_mask.sum(axis=0),
            "ctrl_num": ctrl_phenotype_mask.sum(axis=0),
            "odds_ratio": odds_ratio,
            "p_value": twosided,
            # Add counts
            "case_het_num": case_het_num,
            "case_hom_ref_num": case_hom_ref_num,
            "case_hom_alt_num": case_hom_alt_num,
            "ctrl_het_num": ctrl_het_num,
            "ctrl_hom_ref_num": ctrl_hom_ref_num,
            "ctrl_hom_alt_num": ctrl_hom_alt_num,
            # Add some convenience counts
            "case_ref_num": case_ref_num,
            "case_alt_num": case_alt_num,
            "ctrl_ref_num": ctrl_ref_num,
            "ctrl_alt_num": ctrl_alt_num,
            # Go on then...
            "case_ref_af": case_ref_af,
            "case_alt_af": case_alt_af,
            "ctrl_ref_af": ctrl_ref_af,
            "ctrl_alt_af": ctrl_alt_af,
        }
    )

    return results_df


def run_phewas_or_load_from_cache(
    variant_id: str,
    genotype_service: GenotypeService,
    phenotype_service: PhenotypeService,
    ancestry: str = "EUR",
    no_cache: bool = False,
) -> pd.DataFrame:
    """
    Run PheWAS analysis on a variant and save the results to a cache file.
    If the cache file exists, load the results from the cache file.
    """
    # Parse variant information from 'URL format' (e.g. "5_33951588_C_G")
    parts = variant_id.split("_")
    chrom = parts[0]
    pos = parts[1]
    allele1 = parts[2]
    allele2 = parts[3]

    # De-Nom'alize the variant format
    plnkified_variant_id = f"{chrom}:{pos}_{allele1}/{allele2}"

    cache_file = Config.PHEWAS_PHENO_DIR / f"phewas_results_{variant_id}_{ancestry}.tsv"
    if cache_file.exists() and not no_cache:
        return pd.read_csv(cache_file, sep="\t", dtype={"phecode": str})
    else:
        results = run_phewas(
            plnkified_variant_id, genotype_service, phenotype_service, ancestry
        )
        results.to_csv(cache_file, sep="\t")
        return results


def run_full_analysis(
    genotype_service: GenotypeService,
    phenotype_service: PhenotypeService,
    ancestry: str = "EUR",
):
    """Run PheWAS analysis on all variants."""
    from multiprocessing import Pool

    # NOTE: N
    variants = genotype_service.plink_variant_ids

    def run_one_variant_id(variant_id):
        return run_phewas(variant_id, genotype_service, phenotype_service, ancestry)

    with Pool(96) as p:
        p.map(run_one_variant_id, variants)


# @profile
def main():
    """Debug entry point for blueprint development."""

    test_variant1 = "19:15373898:C:T"
    test_variant2 = "5:33951588_C/G"
    test_variant3 = "9_32448929_C_T"

    from data_services import ServiceRegistry

    registry = ServiceRegistry.from_config(Config)

    genotype_service = registry.genotype
    phenotype_service = registry.phenotype

    print(f"Performing PheWAS on variant {test_variant3}")
    results = run_phewas_or_load_from_cache(
        test_variant3,
        genotype_service,
        phenotype_service,
        ancestry="EUR",
        # no_cache=True,
    )

    print(results.iloc[0].T)

    formatted_results = format_phewas_results(results)

    print(formatted_results[0])

    exit(0)


if __name__ == "__main__":
    main()
    print("Done!")
