import logging
import os
from collections import Counter
from typing import Any, List, NamedTuple

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from tqdm import tqdm

from blueprints.nomaly import PhenotypesHDF5
from blueprints.nomaly_services import services
from config import Config
from db import get_all_phecodes

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
            ref_alleles = item.genotype
            alt_alleles = 2 - item.genotype
            table[item.case, 0] += count * alt_alleles
            table[item.case, 1] += count * ref_alleles

            if item.case == 1:  # Case
                ref_allele_count_cases += count * ref_alleles
                alt_allele_count_cases += count * alt_alleles
                if item.genotype == 2:
                    homozygous_ref_cases += count
                elif item.genotype == 0:
                    homozygous_alt_cases += count
                elif item.genotype == 1:
                    heterozygous_cases += count
            else:  # Control
                ref_allele_count_controls += count * ref_alleles
                alt_allele_count_controls += count * alt_alleles
                if item.genotype == 2:
                    homozygous_ref_controls += count
                elif item.genotype == 0:
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


def get_genotype_data(variant: str) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Get genotype data for a variant.

    Args:
        variant (str): Variant ID in format "CHR:POS:REF:ALT"

    Returns:
        tuple: (sorted_eids, sorted_genotypes) or None if error
    """
    genotype_service = services.genotype
    if genotype_service is None:
        raise ValueError("Genotype service is not initialized!")

    try:
        # Get genotype data
        genotype_eids = genotype_service.individual
        genotypes = genotype_service.query_variants(variant)
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

        # Sort the data
        sorted_indices = np.argsort(genotype_eids)
        sorted_genotype_eids = genotype_eids[sorted_indices]
        sorted_genotypes = genotypes[sorted_indices]

        return sorted_genotype_eids, sorted_genotypes

    except Exception as e:
        print(f"Error in get_genotype_data for variant {variant}: {str(e)}")
        return None


def process_phecode(
    phecode, sorted_genotype_eids, sorted_genotypes, phenotype_data, all_phecodes
):
    try:
        eids, cases = phenotype_data.get_cases_for_phecode(phecode)
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


def phecode_level_assoc(variant: str) -> pd.DataFrame:
    """
    Run PheWAS analysis for a variant and save results to file.
    Returns the DataFrame of results.
    """
    # Get genotype data and handle failure case
    genotype_result = get_genotype_data(variant)
    if genotype_result is None:
        print(f"No genotype data found for variant {variant}")
        return pd.DataFrame()

    sorted_genotype_eids, sorted_genotypes = genotype_result

    all_phecodes = get_all_phecodes()
    phenotype_data = PhenotypesHDF5()
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


def get_phewas_results(variant: str, phecode: str | None = None) -> pd.DataFrame | None:
    """
    Get PheWAS results for a variant, either from cache or by running analysis.
    If phecode is provided and no cached results exist, only analyze that specific phecode.

    Args:
        variant (str): Variant ID in format "CHR:POS:REF:ALT"
        phecode (str | None): Optional phecode to analyze

    Returns:
        pd.DataFrame | None: DataFrame with PheWAS results or None if analysis fails
    """
    # Convert variant format for filename
    variant_underscore = variant.replace(":", "_")
    phewas_file = PHEWAS_PHENO_DIR / f"variant_{variant_underscore}.assoc_nomaly.tsv"

    # Check if full results already exist
    if os.path.exists(phewas_file):
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
            genotype_result = get_genotype_data(variant)
            if genotype_result is None:
                return None

            sorted_genotype_eids, sorted_genotypes = genotype_result
            all_phecodes = get_all_phecodes()
            phenotype_data = PhenotypesHDF5()

            # Process just the requested phecode
            result = process_phecode(
                phecode,
                sorted_genotype_eids,
                sorted_genotypes,
                phenotype_data,
                all_phecodes,
            )

            if result:
                return pd.DataFrame([result])
            return None

        except Exception as e:
            logger.error(
                f"Error in single-phecode analysis for variant '{variant}': {e}"
            )
            return None

    # No cached results and no specific phecode - run full analysis
    try:
        results_df = phecode_level_assoc(variant)
        if results_df is None or results_df.empty:
            print(f"PheWAS analysis returned no results for variant {variant}")
            return None
        return results_df
    except Exception as e:
        print(f"Error running PheWAS for variant {variant}: {e}")
        return None


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
    variant_id: str, phecode: str | None = None
) -> List[dict[str, Any]]:
    """
    Get formatted PheWAS data for a variant, optionally filtered by phecode.

    Args:
        variant_id (str): Variant ID
        phecode (str | None): Optional phecode to filter results

    Returns:
        List[dict]: List of formatted PheWAS results (one dict per row)
    """
    phewas_df = get_phewas_results(variant_id, phecode)

    if phewas_df is None or phewas_df.empty:
        print(f"No PheWAS results available for variant {variant_id}")
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


def main():
    test_variant = "11:69083946:T:C"
    test_variant = "9:35066710:A:G"
    test_variant = "19:44908684:C:T"
    test_variant = "8:7055492:C:T"  # Just kill me
    test_variant = "8:6870776:T:C"
    test_variant = "19:15373898:C:T"

    # phecode_level_assoc(test_variant)

    # Some setup to do outside of the app context
    from nomaly_services import services
    from nomaly import GenotypeHDF5

    global services
    services.genotype = GenotypeHDF5(Config.GENOTYPES_H5)

    get_phewas_results(test_variant, "642.1")

    DO_FULL_RUN = False
    if DO_FULL_RUN:
        from db import get_all_variants

        variants_df = get_all_variants()
        variants_df["variant_id_standard"] = variants_df.variant_id.apply(
            lambda x: x.replace("/", "_")
        )

        from multiprocessing import Pool

        # Create a pool of 10 workers
        with Pool(96) as p:
            # Map get_phewas_results to all variants
            p.map(process_variant, variants_df.variant_id_standard)
            print("Done!")


if __name__ == "__main__":
    main()
    print("Done!")
