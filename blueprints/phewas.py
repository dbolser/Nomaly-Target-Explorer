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

        # Format the table as regular integers
        print(table.astype(int))

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
    genotype_matrix = nomaly_genotype.query_variants(test_variant)
    
    # Convert eids to integers and create a mapping array
    eid_array = genotype_eids.astype(int)
    genotype_array = genotype_matrix[0]  # First (and only) variant
    
    all_phecodes = get_all_phecodes()
    phenotype_data = PhenotypesHDF5()

    for phecode in tqdm(all_phecodes.phecode, desc=f"Counting cases for each PheCode for {test_variant}"):
        # Get case data for current phecode
        try:    
            eids_eur, cases_eur = phenotype_data.get_cases_for_phecode(phecode, "EUR")
        except ValueError:
            print(f"Phecode {phecode} not found in the data matrix")
            continue
               
        # Find indices where eids_eur exists in eid_array
        mask = np.isin(eid_array, eids_eur)
        indices = np.where(mask)[0]
        
        # Get corresponding genotypes and cases
        matched_genotypes = genotype_array[indices]

        phecode_counts = PhecodeCounts()

        # Count combinations using numpy
        for genotype in np.unique(matched_genotypes):
            for case in [0, 1]:  # Assuming binary case/control
                count = np.sum((matched_genotypes == genotype) & (cases_eur == case))
                if count > 0:
                    phecode_counts.add(case, int(genotype), count)
        print(f"Phecode: {phecode}")
        print(phecode_counts)
        odds_ratio, p_value = phecode_counts.get_stats()
        print(f"Odds ratio: {odds_ratio}, p-value: {p_value}")

    print("Done!")


def main():
    test_variant = "11:69083946:T:C"

    phecode_level_assoc(test_variant)

if __name__ == '__main__':
    main()
    print("Done!")