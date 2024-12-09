from blueprints.nomaly import nomaly_genotype, PhenotypesHDF5
from db import get_all_phecodes

from collections import Counter
from typing import NamedTuple

from tqdm import tqdm

import pandas as pd
import numpy as np

from scipy.stats import fisher_exact

PHEWAS_PHENO_DIR = '/data/clu/ukbb/by_variant/'
UKBB_PHENO_DIR = '/data/general/UKBB/Phenotypes/'


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

    def get_stats(self) -> float:
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


def get_genotype_data(variant: str):
    genotype_eids = nomaly_genotype.individual.astype(int)
    genotypes = nomaly_genotype.query_variants(variant)[0]
    sorted_indices = np.argsort(genotype_eids)
    sorted_genotype_eids = genotype_eids[sorted_indices]
    sorted_genotypes = genotypes[sorted_indices]


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
                phecode_counts.add(case, int(genotype), count)

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
    sorted_genotype_eids, sorted_genotypes = get_genotype_data(variant)
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

    return pd.DataFrame(data_to_return)


def main():
    test_variant = "11:69083946:T:C"
    test_variant = "9:35066710:A:G"
    test_variant = "19:44908684:C:T"

    phecode_level_assoc(test_variant)


if __name__ == "__main__":
    main()
    print("Done!")
