import os
import pickle
import subprocess
import pandas as pd

# ----------------------------------------------------- #
# GLOBAL VARIABLES
# ----------------------------------------------------- #
GWAS_PHENO_DIR = '/data/clu/ukbb/by_pheno/'
UKBB_PHENO_DIR = '/data/general/UKBB/Phenotypes/'

# ----------------------------------------------------- #
# LOCAL
# ----------------------------------------------------- #
# source_plink_genome = '/data/general/UKBB/Run-v1/DatabaseInputs/genotypes' # the original plink files
source_plink_genome = '/data/clu/ukbb/genotypes_nomaly'
plink_binary = '/data/clu/ukbb/plink' # PLINK v1.9.0-b.7.6 64-bit (13 Oct 2024)

def variant_level_assoc(pheno_type, code) -> pd.DataFrame:

    # ------------------ Load the cases ------------------ #
    if pheno_type == 'PheCode':
        phecode = code
        with open(UKBB_PHENO_DIR + f'phecode_cases_excludes/phecode_{phecode}.pkl', 'rb') as f:
            cases = pickle.load(f)
        output_prefix = f'phecode_{phecode}'
    else:
        raise ValueError('type not recognized')
    

    # -------------------- run GWAS with plink ------------------ #
    output_path = f'{GWAS_PHENO_DIR}{output_prefix}'

    # make new fam file if not exists
    import pandas as pd
    if not os.path.exists(f'{output_path}.fam'):
        fam = pd.read_csv(f'{source_plink_genome}.fam', header=None, sep=r'\s+')
        fam.columns = ['FID', 'IID', 'Father', 'Mother','sex', 'phenotype'] 
        # Sex code ('1' = male, '2' = female, '0' = unknown)
        # Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
        fam['phenotype'] = 1
        fam.loc[fam['IID'].isin(cases['cases']), 'phenotype'] = 2
        if cases['Sex'] == 'Female':
            fam.loc[fam['sex'] == 1, 'phenotype'] = -9 # if cases only applies to Female, set all Male to missing
        elif cases['Sex'] == 'Male':
            fam.loc[fam['sex'] == 2, 'phenotype'] = -9 # if cases only applies to Male, set all Female to missing
        if cases['exclude']:
            fam.loc[fam['IID'].isin(cases['exclude']), 'phenotype'] = -9
        fam.to_csv(f'{output_path}.fam', sep=' ', header=False, index=False)
        print('new fam file created')
    else:
        print('fam file exists')

    if not os.path.exists(f'{output_path}.assoc'):
        _ = subprocess.run(f'{plink_binary} --bed {source_plink_genome}.bed --bim {source_plink_genome}.bim --fam {output_path}.fam --assoc --out {output_path} --silent', shell=True)
        assoc_path = f'{output_path}.assoc'
        if not os.path.exists(assoc_path):
            print('GWAS failed') # ERROR
            return None
        else:
            print('GWAS done')
    else:
        print('GWAS results exist')

    # -------------------- Load information ------------------ #
    nomaly_variants = pd.read_csv('/data/clu/ukbb/nomaly_variants.tsv', sep='\t')

    # -------------------- parse GWAS results ------------------ #

    assoc = pd.read_csv(assoc_path, sep=r'\s+', dtype={'CHR':str})
    # reconcile the Chromosome: 23:X, 24:Y, 25:XY, 26:MT
    assoc['CHR'] = assoc['CHR'].replace({'23':'X', '24':'Y', '25':'XY', '26':'MT'})
    # print('chromosomes are', assoc['CHR'].unique(), 'total unique variants', assoc.shape[0])

    assoc['CHR_BP_A1_A2'] = [f'{x}:{y}_{A1}/{A2}' for x, y, A1, A2 in zip(assoc.CHR, assoc.BP, assoc.A1, assoc.A2)]
    assoc = assoc.merge(nomaly_variants[['CHR_BP_A1_A2', 'gene_id', 'nomaly_variant', 'RSID']], on='CHR_BP_A1_A2', how='left')
    # print('unique variants processed by Nomaly', assoc.nomaly_variant.nunique())

    assoc = assoc[['nomaly_variant', 'gene_id', 'RSID', 'CHR_BP_A1_A2', 'F_A', 'F_U', 'OR', 'P']].sort_values('P')
    assoc.to_csv(f'{assoc_path}_nomaly.tsv', sep='\t', index=False)
    return assoc

