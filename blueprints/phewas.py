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
        if genotype >= 0: # Ignore -1 (missing)
            self.counts[PheWASItem(case, genotype)] += count

    def get_stats(self) -> float:
        """Calculate the odds ratio and p-value for the PheWAS results

            The code is a bit hairy, but:

                Counts
                 | alt | ref |
            Ctrl |   a |   b |
            Case |   c |   d |

            Odds ratio = (a/b) / (c/d)
            p-value = fisher_exact(table)
        """
        table = np.zeros((2, 2))
        for item, count in self.counts.items():
            ref_alleles = item.genotype
            alt_alleles = 2 - item.genotype
            table[item.case, 0] += count * alt_alleles
            table[item.case, 1] += count * ref_alleles

        return fisher_exact(table)

    def __str__(self):
        string = ""
        for item, count in sorted(self.counts.items()):
            string += f"Case={item.case}, Genotype={item.genotype}: {count}\n"
        return string


def phecode_level_assoc(variant: str) -> pd.DataFrame:
    output_prefix = f'phecode_{variant}'

    output_path = f'{PHEWAS_PHENO_DIR}{output_prefix}'

    # Get genotype data
    genotype_eids = nomaly_genotype.individual

    genotype_matrix = nomaly_genotype.query_variants(variant)

    # Convert eids to integers and create a mapping array
    genotype_eids_array = genotype_eids.astype(int)
    genotype_array = genotype_matrix[0]  # First (and only) variant

    all_phecodes = get_all_phecodes()
    phenotype_data = PhenotypesHDF5()

    data_to_return = []

    for phecode in tqdm(all_phecodes.phecode, desc=f"Counting cases for each PheCode for {variant}"):
        # Get case data for current phecode
        try:    
            eids, cases = phenotype_data.get_cases_for_phecode(phecode)
        except ValueError:
            print(f"Phecode {phecode} not found in the data matrix")
            continue

        # Find genotypes for eids
        mask = np.isin(genotype_eids_array, eids)

        # Get corresponding genotypes and cases
        matched_genotypes = genotype_array[mask]

        phecode_counts = PhecodeCounts()

        # Count combinations using numpy
        for genotype in np.unique(matched_genotypes):
            for case in [0, 1]:  # Using binary case/control
                count = np.sum((matched_genotypes == genotype) & (cases == case))
                if count > 0:
                    phecode_counts.add(case, int(genotype), count)
        odds_ratio, p_value = phecode_counts.get_stats()

        data_to_return.append(
            {
                "phecode": phecode,
                "description": all_phecodes[
                    all_phecodes.phecode == phecode
                ].description.values[0],
                "phecode_group": all_phecodes[
                    all_phecodes.phecode == phecode
                ].phecode_group.values[0],
                "sex": all_phecodes[all_phecodes.phecode == phecode].sex.values[0],
                "n_cases": phecode_counts.case_total,
                "n_controls": phecode_counts.control_total,
                # "n_cases_ref": phecode_counts.case_ref_total,
                # "n_cases_alt": phecode_counts.case_alt_total,
                # "n_controls_ref": phecode_counts.control_ref_total,
                # "n_controls_alt": phecode_counts.control_alt_total,
                # "cases_alt_frequency": phecode_counts.case_alt_total
                # / phecode_counts.case_total,
                # "controls_alt_frequency": phecode_counts.control_alt_total
                # / phecode_counts.control_total,
                "odds_ratio": odds_ratio,
                "p_value": p_value,
            }
        )

    return pd.DataFrame(data_to_return)


def main():
    test_variant = "11:69083946:T:C"
    phecode_level_assoc(test_variant)

if __name__ == '__main__':
    main()
    print("Done!")
