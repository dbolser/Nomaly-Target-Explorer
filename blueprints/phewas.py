"""

Calculate PheWAS results for a variant using about the most cumbersome and
stupid method possible.

Two matrix multiplications is probably all we need (see
'try_simplified_verssion' for wip code).

"""

import logging
import os
from collections import Counter
from typing import Any, List, NamedTuple

import numpy as np
import pandas as pd
from fisher import pvalue_npy
from flask import current_app
from scipy.stats import fisher_exact
from tqdm import tqdm

from config import Config
from data_services.phenotype import PhenotypeService
from db import get_all_phecodes, get_all_variants

# Create a 'dummy' profile decorator if we don't have line_profiler installed
try:
    from line_profiler import profile  # type: ignore
except ImportError:

    def profile(func):
        return func


logger = logging.getLogger(__name__)


PHEWAS_PHENO_DIR = Config.PHEWAS_PHENO_DIR


# New style of declaring a named tuple
class PheWASItem(NamedTuple):
    case: int
    genotype: int


class PhecodeCounts:
    def __init__(self, variant: str = "AB") -> None:
        self.counts = Counter()
        self.variant = variant

    @property
    def total(self) -> int:
        return self.counts.total()

    @property
    def case_total(self) -> int:
        return sum(count for item, count in self.counts.items() if item.case == 1)

    @property
    def control_total(self) -> int:
        return sum(count for item, count in self.counts.items() if item.case == 0)

    def add(self, case: int, genotype: int, count: int = 1) -> None:
        if genotype > 2:
            print(f"Warning: genotype {genotype} is greater than 2")
            return
        if genotype >= 0:  # Ignore -1 (missing)
            self.counts[PheWASItem(case, genotype)] += count

    # @profile
    def get_stats(self) -> dict:
        """Calculate the odds ratio and p-value for the PheWAS results

        The code is a bit hairy, but:

             | Counts    |
             | alt | ref |
        Ctrl |   a |   b |
        Case |   c |   d |

        Odds ratio = (a/b) / (c/d)
        p-value = fisher_exact(table)

        and additional statistics for allele frequencies and genotype counts.

        """
        table = np.zeros((2, 2)).astype(int)
        ref_allele_count_cases = 0
        alt_allele_count_cases = 0
        ref_allele_count_controls = 0
        alt_allele_count_controls = 0
        homozygous_ref_cases = 0
        homozygous_alt_cases = 0
        heterozygous_cases = 0
        homozygous_ref_controls = 0
        homozygous_alt_controls = 0
        heterozygous_controls = 0

        for item, count in self.counts.items():
            ref_alleles = 2 - item.genotype
            alt_alleles = item.genotype
            table[item.case, 0] += count * alt_alleles
            table[item.case, 1] += count * ref_alleles

            if item.case == 1:  # Case
                ref_allele_count_cases += count * ref_alleles
                alt_allele_count_cases += count * alt_alleles
                if item.genotype == 0:
                    homozygous_ref_cases += count
                elif item.genotype == 2:
                    homozygous_alt_cases += count
                elif item.genotype == 1:
                    heterozygous_cases += count
            else:  # Control
                ref_allele_count_controls += count * ref_alleles
                alt_allele_count_controls += count * alt_alleles
                if item.genotype == 0:
                    homozygous_ref_controls += count
                elif item.genotype == 2:
                    homozygous_alt_controls += count
                elif item.genotype == 1:
                    heterozygous_controls += count

        # Calculate allele frequencies
        total_alleles_cases = ref_allele_count_cases + alt_allele_count_cases
        total_alleles_controls = ref_allele_count_controls + alt_allele_count_controls
        ref_allele_freq_cases = (
            ref_allele_count_cases / total_alleles_cases
            if total_alleles_cases > 0
            else 0
        )
        alt_allele_freq_cases = (
            alt_allele_count_cases / total_alleles_cases
            if total_alleles_cases > 0
            else 0
        )
        ref_allele_freq_controls = (
            ref_allele_count_controls / total_alleles_controls
            if total_alleles_controls > 0
            else 0
        )
        alt_allele_freq_controls = (
            alt_allele_count_controls / total_alleles_controls
            if total_alleles_controls > 0
            else 0
        )

        # Return the table and additional statistics
        return {
            "table": table,
            "fisher_stats": fisher_exact(table),
            "ref_allele_freq_cases": ref_allele_freq_cases,
            "alt_allele_freq_cases": alt_allele_freq_cases,
            "ref_allele_freq_controls": ref_allele_freq_controls,
            "alt_allele_freq_controls": alt_allele_freq_controls,
            "homozygous_ref_cases": homozygous_ref_cases,
            "homozygous_alt_cases": homozygous_alt_cases,
            "heterozygous_cases": heterozygous_cases,
            "homozygous_ref_controls": homozygous_ref_controls,
            "homozygous_alt_controls": homozygous_alt_controls,
            "heterozygous_controls": heterozygous_controls,
        }

    def __str__(self):
        string = ""
        for item, count in sorted(self.counts.items()):
            string += f"Case={item.case}, Genotype={item.genotype}: {count}\n"
        return string


# @profile
def get_genotype_data(variant: str, services) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Get genotype data for a variant.

    Args:
        variant (str): Variant ID in format "CHR:POS:REF:ALT"
        services: Service registry containing required services

    Returns:
        tuple: (sorted_eids, sorted_genotypes) or None if error
    """
    genotype_service = services.genotype
    if genotype_service is None:
        raise ValueError("Genotype service is not initialized!")

    try:
        # Get genotype data
        genotype_eids = genotype_service._hdf.individual
        genotypes = genotype_service._hdf.query_variants(variant)

        # NOTE THE FORMER WILL FLIP! THe latter wont!

        # SOME PAGES NEED TO FILP FOR SOME REASON!

        # TODO: MAKE SUER WE NEVER NEED TO FILP!

        # genotypes2 = genotype_service._hdf.get_genotypes(vids=[variant])

        # Well... this is one way of testing...
        # assert np.all(genotypes == genotypes2)

        if genotypes is None or len(genotypes) == 0:
            print(f"No genotype data found for variant {variant}")
            return None

        # Get first row of genotypes (for single variant)
        genotypes = genotypes[0]
        if len(genotypes) != len(genotype_eids):
            print(
                f"Mismatch between genotypes ({len(genotypes)}) and IDs ({len(genotype_eids)})"
            )
            return None

        # I WOULD REFACTOR THIS, but see try_simplified_verssion

        # Sort the data
        sorted_indices = np.argsort(genotype_eids)
        sorted_genotype_eids = genotype_eids[sorted_indices]
        sorted_genotypes = genotypes[sorted_indices]

        return sorted_genotype_eids, sorted_genotypes

    except Exception as e:
        print(
            f"Error in get_genotype_data for variant (get_genotype_data) {variant}: {str(e)}"
        )
        return None


@profile
def process_phecode(
    phecode: str,
    sorted_genotype_eids: np.ndarray,
    sorted_genotypes: np.ndarray,
    phenotype_service: PhenotypeService,  # Type hint the interface
    all_phecodes: pd.DataFrame,
) -> dict | None:
    """Process a single phecode."""
    try:
        eids, cases = phenotype_service._hdf.get_cases_for_phecode(phecode)
    except ValueError:
        print(f"Phecode {phecode} not found in the data matrix")
        return None

    indices = np.searchsorted(sorted_genotype_eids, eids)
    valid_indices = indices[indices < len(sorted_genotype_eids)]
    matched_genotypes = sorted_genotypes[valid_indices]
    matched_cases = cases

    phecode_counts = PhecodeCounts()

    for genotype in [0, 1, 2]:
        for case in [0, 1]:
            count = np.sum((matched_genotypes == genotype) & (matched_cases == case))
            if count > 0:
                phecode_counts.add(case, genotype, count)

    stats = phecode_counts.get_stats()
    odds_ratio, p_value = stats["fisher_stats"]

    return {
        "phecode": phecode,
        "description": all_phecodes[all_phecodes.phecode == phecode].description.values[
            0
        ],
        "phecode_group": all_phecodes[
            all_phecodes.phecode == phecode
        ].phecode_group.values[0],
        "sex": all_phecodes[all_phecodes.phecode == phecode].sex.values[0],
        "n_cases": phecode_counts.case_total,
        "n_controls": phecode_counts.control_total,
        "n_cases_alt": stats["table"][1, 0],
        "n_cases_ref": stats["table"][1, 1],
        "n_controls_alt": stats["table"][0, 0],
        "n_controls_ref": stats["table"][0, 1],
        "odds_ratio": odds_ratio,
        "p_value": p_value,
        "ref_allele_freq_cases": stats["ref_allele_freq_cases"],
        "alt_allele_freq_cases": stats["alt_allele_freq_cases"],
        "ref_allele_freq_controls": stats["ref_allele_freq_controls"],
        "alt_allele_freq_controls": stats["alt_allele_freq_controls"],
        "homozygous_ref_cases": stats["homozygous_ref_cases"],
        "homozygous_alt_cases": stats["homozygous_alt_cases"],
        "heterozygous_cases": stats["heterozygous_cases"],
        "homozygous_ref_controls": stats["homozygous_ref_controls"],
        "homozygous_alt_controls": stats["homozygous_alt_controls"],
        "heterozygous_controls": stats["heterozygous_controls"],
    }


@profile
def phecode_level_assoc(variant: str, services=None) -> pd.DataFrame:
    """
    Run PheWAS analysis for a variant and save results to file.
    Returns the DataFrame of results.
    """
    # Get genotype data and handle failure case
    genotype_result = get_genotype_data(variant, services)
    if genotype_result is None:
        print(f"No genotype data found for variant (phecode_level_assoc) {variant}")
        return pd.DataFrame()

    sorted_genotype_eids, sorted_genotypes = genotype_result

    all_phecodes = get_all_phecodes()

    services = services if services else current_app.extensions["nomaly_services"]
    phenotype_data = services.phenotype
    assert phenotype_data is not None

    data_to_return = []

    for phecode in tqdm(
        all_phecodes.phecode, desc=f"Counting cases for each PheCode for {variant}"
    ):
        result = process_phecode(
            phecode,
            sorted_genotype_eids,
            sorted_genotypes,
            phenotype_data,
            all_phecodes,
        )
        if result:
            data_to_return.append(result)

    # Create DataFrame
    results_df = pd.DataFrame(data_to_return)

    # Save results to file if we have data
    if not results_df.empty:
        # Convert variant format from CHR:POS:REF:ALT to CHR_POS_REF_ALT
        variant_underscore = variant.replace(":", "_")
        output_prefix = f"variant_{variant_underscore}"
        output_path = PHEWAS_PHENO_DIR / f"{output_prefix}.assoc_nomaly.tsv"

        # Save to file
        results_df.to_csv(output_path, sep="\t", index=False)
        print(f"Saved PheWAS results to {output_path}")

    return results_df


@profile
def get_phewas_results(
    variant: str,
    services,
    phecode: str | None = None,
    no_cache: bool = False,
) -> pd.DataFrame:
    """
    Get PheWAS results for a variant, either from cache or by running analysis.
    If phecode is provided and no cached results exist, only analyze that specific phecode.

    Args:
        variant (str): Variant ID in format "CHR:POS:REF:ALT"
        phecode (str | None): Optional phecode to analyze
        no_cache (bool): Whether to ignore cached results
        services: Service registry containing required services

    Returns:
        pd.DataFrame | None: DataFrame with PheWAS results or None if analysis fails
    """
    # Convert variant format for filename
    variant_underscore = variant.replace(":", "_")
    phewas_file = PHEWAS_PHENO_DIR / f"variant_{variant_underscore}.assoc_nomaly.tsv"

    # Check if full results already exist
    if os.path.exists(phewas_file) and not no_cache:
        try:
            df = pd.read_csv(phewas_file, sep="\t", dtype={"phecode": str})
            if phecode is not None:
                return df[df["phecode"] == phecode]
            return df
        except Exception as e:
            logger.error(
                f"Error reading existing PheWAS file for variant {variant}: {e}"
            )
            return pd.DataFrame()

    # If no cached results and specific phecode requested, analyze just that phecode
    if phecode is not None:
        try:
            print(
                f"Running single-phecode analysis for variant {variant} and phecode {phecode}"
            )
            # Get genotype data
            genotype_result = get_genotype_data(variant, services)
            if genotype_result is None:
                return pd.DataFrame()

            sorted_genotype_eids, sorted_genotypes = genotype_result
            all_phecodes = get_all_phecodes()

            # Process just the requested phecode
            result = process_phecode(
                phecode,
                sorted_genotype_eids,
                sorted_genotypes,
                services.phenotype,
                all_phecodes,
            )

            if result:
                result_df = pd.DataFrame([result])
                result_df["ref_allele"] = variant.split(":")[2]
                result_df["alt_allele"] = variant.split(":")[3]
                return result_df
            return pd.DataFrame()

        except Exception as e:
            logger.error(
                f"Error in single-phecode analysis for variant '{variant}': {e}"
            )
            return pd.DataFrame()

    # No cached results and no specific phecode - run full analysis
    try:
        results_df = phecode_level_assoc(variant, services)

        if results_df is None or results_df.empty:
            print(f"PheWAS analysis returned no results for variant {variant}")
            return pd.DataFrame()

        results_df["ref_allele"] = variant_underscore.split("_")[2]
        results_df["alt_allele"] = variant_underscore.split("_")[3]

        return results_df
    except Exception as e:
        print(f"Error running PheWAS for variant {variant}: {e}")
        return pd.DataFrame()


# @profile
def format_phewas_row_for_display(row: pd.Series) -> dict:
    """
    Format a PheWAS result row for display.
    Returns a dictionary with HTML-formatted values.
    """
    return {
        "Phecode": row["phecode"],
        "Description": row["description"],
        "Sex": row["sex"],
        "Group": row["phecode_group"],
        "P": f"{row['p_value']:.2e}<br/>",
        "OR": f"{row['odds_ratio']:.2f}<br/>",
        "Counts": f"{row['n_cases']}<br/>{row['n_controls']}",
        "RefAF": f"{row['ref_allele_freq_cases']:.5f}<br/>{row['ref_allele_freq_controls']:.5f}",
        "AltAF": f"{row['alt_allele_freq_cases']:.5f}<br/>{row['alt_allele_freq_controls']:.5f}",
        "Ref_HMOZ": f"{row.get('homozygous_ref_cases', 0)}<br/>{row.get('homozygous_ref_controls', 0)}",
        "Alt_HMOZ": f"{row.get('homozygous_alt_cases', 0)}<br/>{row.get('homozygous_alt_controls', 0)}",
        "HTRZ": f"{row.get('heterozygous_cases', 0)}<br/>{row.get('heterozygous_controls', 0)}",
        # Keep raw values for filtering
        "p_value": row["p_value"],
    }


def get_formatted_phewas_data(
    variant_id: str, services, phecode: str | None = None, no_cache: bool = False
) -> List[dict[str, Any]]:
    """
    Get formatted PheWAS data for a variant, optionally filtered by phecode.

    Args:
        variant_id (str): Variant ID
        phecode (str | None): Optional phecode to filter results
        services: Service registry containing required services

    Returns:
        List[dict]: List of formatted PheWAS results (one dict per row)
    """
    phewas_df = get_phewas_results(variant_id, services, phecode, no_cache)

    if phewas_df is None or phewas_df.empty:
        print(
            f"No PheWAS results available for variant (get_formatted_phewas_data) {variant_id}"
        )
        return []

    if phecode is not None:
        # Filter for specific phecode
        phewas_df = phewas_df[phewas_df["phecode"] == phecode]
        if phewas_df.empty:
            print(
                f"PheWAS row not found for variant {variant_id} and phecode {phecode}"
            )
            return []

    # Format all rows (whether filtered or not)
    return [format_phewas_row_for_display(row) for _, row in phewas_df.iterrows()]


def process_variant(variant: str):
    return get_phewas_results(variant, None)


def run_full_analysis():
    """Run PheWAS analysis on all variants."""
    from multiprocessing import Pool

    variants_df = get_all_variants()
    variants_df["variant_id_standard"] = variants_df.variant_id.apply(
        lambda x: x.replace("/", "_")
    )

    with Pool(96) as p:
        p.map(process_variant, variants_df.variant_id_standard)


def try_simplified_version(variant: str) -> pd.DataFrame:
    # 1) Data loading:

    from data_services import ServiceRegistry

    registry = ServiceRegistry.from_config(Config)

    genotype_service = registry.genotype
    phenotype_service = registry.phenotype

    if not genotype_service or not phenotype_service:
        raise ValueError("Genotype or phenotype service not initialized")

    phenotype_eids = phenotype_service.eids
    phenotype_phec = phenotype_service.phecodes
    phenotype_matrix = phenotype_service.get_phenotypes()  # Get them all

    # Important to pass phenotype_eids here
    genotype_matrix = genotype_service.get_genotypes(
        eids=phenotype_eids, vids=np.array([variant])
    )

    assert phenotype_matrix.shape == (486145, 1855)
    assert genotype_matrix.shape == (1, phenotype_eids.shape[0])

    # TODO: Somewhere deep down, it's implicit that the underlying EIDs are the
    # same...

    # Moving on...

    # OK, this gets a bit slower, but, compared to the original crap, it's lightning fast...
    ctrl_phenotype_mask = phenotype_matrix == 0
    case_phenotype_mask = phenotype_matrix == 1

    het_genotype_mask = genotype_matrix == 1
    ref_genotype_mask = genotype_matrix == 0
    alt_genotype_mask = genotype_matrix == 2

    # Get all the counts in two stesp
    ctrl_het_mask = ctrl_phenotype_mask * het_genotype_mask.T
    ctrl_ref_mask = ctrl_phenotype_mask * ref_genotype_mask.T
    ctrl_alt_mask = ctrl_phenotype_mask * alt_genotype_mask.T
    case_het_mask = case_phenotype_mask * het_genotype_mask.T
    case_ref_mask = case_phenotype_mask * ref_genotype_mask.T
    case_alt_mask = case_phenotype_mask * alt_genotype_mask.T

    ctrl_het_num = ctrl_het_mask.sum(axis=0)
    ctrl_ref_num = ctrl_ref_mask.sum(axis=0) * 2 + ctrl_het_num
    ctrl_alt_num = ctrl_alt_mask.sum(axis=0) * 2 + ctrl_het_num
    case_het_num = case_het_mask.sum(axis=0)
    case_ref_num = case_ref_mask.sum(axis=0) * 2 + case_het_num
    case_alt_num = case_alt_mask.sum(axis=0) * 2 + case_het_num

    case_ref_psu = case_ref_num + 1
    case_alt_psu = case_alt_num + 1
    ctrl_ref_psu = ctrl_ref_num + 1
    ctrl_alt_psu = ctrl_alt_num + 1

    # case_ref_af = case_ref_psu / (case_ref_psu + case_alt_psu)
    # case_alt_af = case_alt_psu / (case_ref_psu + case_alt_psu)

    # assert np.all(case_ref_af + case_alt_af == 1)

    # Calculate odds ratios
    odds_ratio = (case_alt_psu * ctrl_ref_psu) / (case_ref_psu * ctrl_alt_psu)

    # odds_ratio_long = (case_alt_psu / case_ref_psu) / (ctrl_alt_psu / ctrl_ref_psu)
    # np.allclose(odds_ratio, odds_ratio_long)

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
            "phecodes": phenotype_phec,
            "odds_ratios": odds_ratio,
            "p_values": twosided,
            "ref_case": case_ref_num,
            "alt_case": case_alt_num,
            "ref_ctrl": ctrl_ref_num,
            "alt_ctrl": ctrl_alt_num,
        }
    )

    return results_df


@profile
def main():
    """Debug entry point for blueprint development."""

    test_variant1 = "19:15373898:C:T"
    test_variant2 = "5:33951588:C:G"

    _ = try_simplified_version(test_variant2)

    from app import create_app

    app = create_app("development")

    services = app.extensions["nomaly_services"]

    with app.app_context():
        # Quick test of single variant/phecode
        results = get_phewas_results(
            test_variant1, services, phecode="642.1", no_cache=True
        )
        print(results.T)

        results = get_phewas_results(test_variant1, services, no_cache=True)
        print(results)

        # Optionally run full analysis
        DO_FULL_RUN = False
        if DO_FULL_RUN:
            run_full_analysis()


if __name__ == "__main__":
    main()
    print("Done!")
