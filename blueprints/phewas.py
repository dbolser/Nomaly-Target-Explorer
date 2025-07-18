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
        # NOTE: The Sex, description and phecode_group columns all come from the
        # database, not the PheWAS results.
        "Sex": row.get("sex", "Unknown"),
        "Description": row.get("description", "Unknown"),
        "Group": row.get("phecode_group", "Unknown"),
        # NOTE: The rest of the columns come from the PheWAS results.
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
    """Format the results, row by row for display.

    Return a list of dics, one per row.
    """
    # Sort the results by p-value for display...
    phewas_df = phewas_df.sort_values(by="p_value", ascending=True)
    return [format_phewas_row(row) for _, row in phewas_df.iterrows()]


def run_phewas(
    variant: str,
    genotype_service: GenotypeService,
    phenotype_service: PhenotypeService,
    ancestry: str = "EUR",
    sanity: bool = True,
) -> pd.DataFrame:
    phecodes = phenotype_service.phecodes
    phenotype_eids = phenotype_service.get_eids(ancestry=ancestry)
    phenotype_matrix = phenotype_service.get_phenotypes(ancestry=ancestry)

    # Important to pass phenotype_eids here
    genotype_matrix = genotype_service.get_genotypes(
        eids=phenotype_eids, vids=np.array([variant])
    )

    if not sanity:
        assert phenotype_matrix.shape == (phenotype_eids.shape[0], phecodes.shape[0])
        assert genotype_matrix.shape == (1, phenotype_eids.shape[0])

    # TODO: Somewhere deep down, it's implicit that the underlying EIDs are all
    #       and always the same...

    # Moving on...

    # Get the individual phenotype masks
    case_phenotype_mask = phenotype_matrix == 1
    ctrl_phenotype_mask = phenotype_matrix == 0
    # skip_phenotype_mask = phenotype_matrix == 9  # For book keeping...

    # Get the individual genotype masks
    het_genotype_mask = genotype_matrix == 1
    ref_genotype_mask = genotype_matrix == 0
    alt_genotype_mask = genotype_matrix == 2

    # mis_genotype_mask = genotype_matrix == -9  # For book keeping...

    # Now we get all the counts in two steps...

    # 1) Get the composite mask for all combinations
    case_het_mask = case_phenotype_mask * het_genotype_mask.T
    case_ref_mask = case_phenotype_mask * ref_genotype_mask.T
    case_alt_mask = case_phenotype_mask * alt_genotype_mask.T
    ctrl_het_mask = ctrl_phenotype_mask * het_genotype_mask.T
    ctrl_ref_mask = ctrl_phenotype_mask * ref_genotype_mask.T
    ctrl_alt_mask = ctrl_phenotype_mask * alt_genotype_mask.T

    # For book keeping...
    # case_mis_mask = case_phenotype_mask * mis_genotype_mask.T
    # ctrl_mis_mask = ctrl_phenotype_mask * mis_genotype_mask.T

    # 2) Get the counts for each combination
    case_het_num = case_het_mask.sum(axis=0)
    case_hom_ref_num = case_ref_mask.sum(axis=0)
    case_hom_alt_num = case_alt_mask.sum(axis=0)

    ctrl_het_num = ctrl_het_mask.sum(axis=0)
    ctrl_hom_ref_num = ctrl_ref_mask.sum(axis=0)
    ctrl_hom_alt_num = ctrl_alt_mask.sum(axis=0)

    # For book keeping...
    # case_mis_num = case_mis_mask.sum(axis=0)
    # ctrl_mis_num = ctrl_mis_mask.sum(axis=0)

    # Sum alleles here for convenience
    case_ref_num = case_hom_ref_num * 2 + case_het_num
    case_alt_num = case_hom_alt_num * 2 + case_het_num
    ctrl_ref_num = ctrl_hom_ref_num * 2 + ctrl_het_num
    ctrl_alt_num = ctrl_hom_alt_num * 2 + ctrl_het_num

    # Calculate allele frequencies here for convenience
    with np.errstate(divide="ignore", invalid="ignore"):
        # TODO: We get zeros in the denominator due to missing genotypes...
        #       At least that's what I think... (See book keeping above).
        case_ref_af = np.divide(case_ref_num, case_ref_num + case_alt_num)
        case_alt_af = np.divide(case_alt_num, case_ref_num + case_alt_num)
        ctrl_ref_af = np.divide(ctrl_ref_num, ctrl_ref_num + ctrl_alt_num)
        ctrl_alt_af = np.divide(ctrl_alt_num, ctrl_ref_num + ctrl_alt_num)

    if not sanity:
        # Create masks for non-NaN values
        case_valid = ~np.isnan(case_ref_af) & ~np.isnan(case_alt_af)
        ctrl_valid = ~np.isnan(ctrl_ref_af) & ~np.isnan(ctrl_alt_af)

        # Sanity check only on positions with valid values
        assert np.allclose(case_ref_af[case_valid] + case_alt_af[case_valid], 1)
        assert np.allclose(ctrl_ref_af[ctrl_valid] + ctrl_alt_af[ctrl_valid], 1)

    # Calculate odds ratios
    with np.errstate(divide="ignore", invalid="ignore"):
        # TODO: We get zeros in the denominator due to missing genotypes...
        #       At least that's what I think... (See book keeping above).
        odds_ratio = (case_alt_num * ctrl_ref_num) / (case_ref_num * ctrl_alt_num)

    # Sanity check...
    if not sanity:
        odds_ratio_long = (case_alt_num / case_ref_num) / (ctrl_alt_num / ctrl_ref_num)

        # Create masks for non-NaN values
        assert np.all(np.isnan(odds_ratio) == np.isnan(odds_ratio_long))
        odds_ratio_valid = ~np.isnan(odds_ratio)

        # Sanity check only on positions with valid values
        assert np.allclose(
            odds_ratio[odds_ratio_valid], odds_ratio_long[odds_ratio_valid]
        )

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
    sanity: bool = True,
) -> pd.DataFrame:
    """
    Run PheWAS analysis on a variant and save the results to a cache file.
    If the cache file exists, load the results from the cache file.
    """

    plnkified_variant_id = variant_id
    flat_variant_id = variant_id

    if ":" not in variant_id:
        # Parse variant information from 'URL format' (e.g. "5_33951588_C_G")
        parts = variant_id.split("_")
        chrom = parts[0]
        pos = parts[1]
        allele1 = parts[2]
        allele2 = parts[3]

        plnkified_variant_id = f"{chrom}:{pos}_{allele1}/{allele2}"
    else:
        flat_variant_id = variant_id.replace(":", "_").replace("/", "_")

    cache_file = (
        Config.PHEWAS_PHENO_DIR / f"phewas_results_{flat_variant_id}_{ancestry}.tsv"
    )
    if cache_file.exists() and not no_cache:
        return pd.read_csv(cache_file, sep="\t", dtype={"phecode": str})
    else:
        results = run_phewas(
            plnkified_variant_id,
            genotype_service,
            phenotype_service,
            ancestry,
            sanity=sanity,
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

    variants = genotype_service.plink_variant_ids

    def run_one_variant_id(variant_id):
        return run_phewas(variant_id, genotype_service, phenotype_service, ancestry)

    with Pool(96) as p:
        p.map(run_one_variant_id, variants)


# @profile
def main():
    """Debug entry point for blueprint development."""

    from data_services import ServiceRegistry

    services = ServiceRegistry.from_config(Config)
    genotype_service = services.genotype
    phenotype_service = services.phenotype

    test_variant = "19:15373898:C:T"
    test_variant = "5:33951588_C/G"
    test_variant = "9_32448929_C_T"
    test_variant = "5_151467708_T_C"
    test_variant = "11_116836316_A_G"

    ancestry = "AFR"

    print(f"Performing PheWAS on variant {test_variant}")
    results = run_phewas_or_load_from_cache(
        test_variant,
        genotype_service,
        phenotype_service,
        ancestry=ancestry,
        no_cache=True,
        sanity=False,
    )

    # DEBUGGING
    results = results[results["phecode"].isin(["008", "627.4", "627.5", "627.6"])]

    print(results.iloc[0].T)

    formatted_results = format_phewas_results(results)

    print(formatted_results[0])

    exit(0)

    # Run some random PheWAS for fun...
    variants = genotype_service.plink_variant_ids
    ancestries = np.array(["AFR", "EAS", "EUR", "SAS"])

    done = set()
    for _ in range(1_000_000):
        variant = np.random.choice(variants)
        ancestry = np.random.choice(ancestries)
        if (variant, ancestry) in done:
            continue
        done.add((variant, ancestry))
        print(f"Running PheWAS for {variant} ({ancestry})")
        try:
            run_phewas_or_load_from_cache(
                variant,
                genotype_service,
                phenotype_service,
                ancestry=ancestry,
                no_cache=True,
            )
            print(f"Successfully ran PheWAS for {variant} ({ancestry})")
        except Exception as e:
            print(f"Error running PheWAS for {variant} ({ancestry}): {e}")
            continue


if __name__ == "__main__":
    main()
    print("Done!")
