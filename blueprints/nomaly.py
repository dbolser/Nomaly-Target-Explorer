import os
import h5py
import numpy as np
import pandas as pd

from functools import cached_property

import plotly.express as px

UKBB_PHENO_DIR = '/data/general/UKBB/Phenotypes/'
icd10_cases_h5 = UKBB_PHENO_DIR + 'ukbbrun_icd10_2024-09-07_any.h5'
phenotypes_h5 = UKBB_PHENO_DIR + 'phecode_cases_excludes.h5'

RESOURCE_DATA_DIR = '/data/general/Data/'
pharos_path = RESOURCE_DATA_DIR + 'pharos_api_query.out'
pp_path = RESOURCE_DATA_DIR + 'pp_by_gene.tsv'
pharos = pd.read_csv(pharos_path, sep="\t", encoding='ISO-8859-1').rename(columns={'#symbol': 'gene'})
pp = pd.read_csv(pp_path, sep="\t", index_col=0).rename(columns={'Description': 'drug_program_indication'})
pp['gene'] = pp.index

# ------------------------------------------------------------------------------#
# Nomaly stats pre-calculated
# ------------------------------------------------------------------------------#

class StatsHDF5:
    def __init__(self, hdf5_file):
        try:
            self.hdf5_file = hdf5_file
            if not os.path.exists(hdf5_file):
                raise FileNotFoundError(f"HDF5 file not found: {hdf5_file}")
            self.f = h5py.File(hdf5_file, 'r')

            # Load the data matrix and index information
            self.data_matrix = self.f['data']
            self.rows = self.f['rows_term'][...].astype(str)  # Load term data
            self.columns = self.f['columns_phecode'][...].astype(str)  # Load phecode data
            self.third_dim = self.f['3rddim_statstype'][...].astype(str)  # Load statstype data

            # rename third_dim that ends with _p_value to _pvalue
            self.third_dim = np.array([x.replace('_p_value', '_pvalue') for x in self.third_dim])
        except Exception as e:
            current_app.logger.error(f'StatsHDF5 initialization error: {e}\n{traceback.format_exc()}')
            raise

    def get_stats_by_disease(self, disease, statstype=None):
        try:
            if disease not in self.columns:
                disease = disease + '.0'
                if disease not in self.columns:
                    raise ValueError(f"Disease {disease} not found in dataset")

            mask_column = self.columns == disease
            mask_column_indices = np.where(mask_column)[0][0]
            if statstype is None:
                # if statstype is not specified, return all statstypes for the disease
                disease_data = self.data_matrix[:, mask_column_indices, :]
                disease_DF = pd.DataFrame(
                    disease_data, index=self.rows.astype(str), columns=self.third_dim
                )
                return disease_DF
            elif isinstance(statstype, str) and statstype in self.third_dim:
                # if statstype is specified, return the data for the statstype
                mask_3rd_dim = self.third_dim == statstype
                mask_3rd_dim_indices = np.where(mask_3rd_dim)[0][0]
                return self.data_matrix[:,mask_column_indices, mask_3rd_dim_indices]
            elif isinstance(statstype, list):
                # if statstype is a list of strings, return the data for the statstypes
                mask_3rd_dim = np.isin(self.third_dim, statstype)
                mask_3rd_dim_indices = np.where(mask_3rd_dim)[0]
                return self.data_matrix[:,mask_column_indices, mask_3rd_dim_indices]
            else:
                raise ValueError(f"statstype must be a string or a list of strings. statstype: {statstype}")
        except Exception as e:
            current_app.logger.error(f'Error getting stats for disease {disease}: {e}')
            raise
    
    def get_stats_by_term_disease(self, term, disease):
        if disease not in self.columns:
            disease = disease + '.0'
            if disease not in self.columns:
                raise ValueError(f"{disease} not found in the columns of the data matrix.")
        
        mask_column = self.columns == disease
        mask_column_indices = np.where(mask_column)[0][0]

        if term not in self.rows:
            raise ValueError(f"{term} not found in the rows of the data matrix.")
        
        mask_row = self.rows == term
        mask_row_indices = np.where(mask_row)[0][0]

        # return a dictionary with the statstype as key and the value as the data
        statsdict = {
            self.third_dim[i]: self.data_matrix[mask_row_indices, mask_column_indices, i] for i in range(len(self.third_dim))
        }
        return statsdict


# stats = StatsHDF5(nomaly_run_stats_h5)
# data = stats.get_stats_by_disease(tocheck['code'], 'mwu_pvalue')
# # save data and row names
# data = pd.DataFrame(data, index=stats.rows.astype(str))
# data.columns = ['mwu_pvalue']
# data = data.sort_values(by='mwu_pvalue')

# ------------------------------------------------------------------------------#
# Genotypes
# ------------------------------------------------------------------------------#

class GenotypeHDF5:
    def __init__(self, hdf5_file):
        try:
            self.hdf5_file = hdf5_file
            if not os.path.exists(hdf5_file):
                raise FileNotFoundError(f"Genotype HDF5 file not found: {hdf5_file}")

            self.f = h5py.File(hdf5_file, "r")

            # Verify required datasets exist
            required_datasets = ["genotype_matrix", "fam", "bim"]
            for dataset in required_datasets:
                if dataset not in self.f:
                    raise KeyError(
                        f"Required dataset '{dataset}' not found in genotype HDF5 file"
                    )

            # Load the genotype matrix and index information
            self.genotype_matrix = self.f["genotype_matrix"]
            self.individual: np.array = self.f["fam"][...].astype(int)
            self.variants = self.f["bim"][...].astype(str)

        except Exception as e:
            print(f"Error initializing GenotypeHDF5: {str(e)}")
            raise

    def __del__(self):
        """Cleanup method to ensure file is closed."""
        try:
            if hasattr(self, "f"):
                self.f.close()
        except Exception as e:
            print(f"Error closing genotype file: {str(e)}")

    def _flip_alleles(self, variant: str) -> str:
        """
        Flip the ref and alt alleles in a variant ID.

        Args:
            variant (str): Variant ID in format "CHR:POS:REF:ALT" or "CHR_POS_REF/ALT"

        Returns:
            str: Flipped variant ID in same format as input
        """
        try:
            if "_" in variant:  # Format: CHR_POS_REF/ALT
                chrom, pos, alleles = variant.split("_")
                if "/" not in alleles:
                    return None
                ref, alt = alleles.split("/")
                return f"{chrom}_{pos}_{alt}/{ref}"
            else:  # Format: CHR:POS:REF:ALT
                chrom, pos, ref, alt = variant.split(":")
                return f"{chrom}:{pos}:{alt}:{ref}"
        except Exception as e:
            print(f"Error flipping alleles for variant {variant}: {str(e)}")
            return None

    def _standardize_variant_format(self, variant: str) -> str | None:
        """
        Convert any variant format to CHR:POS:REF:ALT format.

        Handles formats:
        - CHR_POS_REF/ALT (e.g., "8_6870776_C/T")
        - CHR:POS:REF:ALT (e.g., "8:6870776:C:T")
        - CHR-POS-REF-ALT (e.g., "8-6870776-C-T")
        - CHR.POS.REF.ALT (e.g., "8.6870776.C.T")

        Args:
            variant (str): Variant in any supported format

        Returns:
            str|None: Standardized variant string or None if invalid format
        """
        try:
            if not isinstance(variant, str):
                print(f"Invalid variant type: {type(variant)}, expected string")
                return None

            # First, extract the components regardless of format
            if "_" in variant and "/" in variant:  # CHR_POS_REF/ALT
                chrom, pos, alleles = variant.split("_")
                ref, alt = alleles.split("/")
            elif ":" in variant:  # CHR:POS:REF:ALT
                chrom, pos, ref, alt = variant.split(":")
            elif "-" in variant:  # CHR-POS-REF-ALT
                chrom, pos, ref, alt = variant.split("-")
            elif "." in variant:  # CHR.POS.REF.ALT
                chrom, pos, ref, alt = variant.split(".")
            else:
                print(f"Unrecognized variant format: {variant}")
                return None

            # Validate components
            if not chrom or not pos or not ref or not alt:
                print(f"Missing component in variant: {variant}")
                return None

            if not pos.isdigit():
                print(f"Position must be numeric: {pos}")
                return None

            # Clean chromosome format (e.g., 'chr8' -> '8')
            chrom = chrom.lower().replace("chr", "")

            # Standardize to colon format
            return f"{chrom}:{pos}:{ref}:{alt}"

        except Exception as e:
            print(f"Error standardizing variant format {variant}: {str(e)}")
            return None

    def query_variants(self, variant: str) -> np.ndarray | None:
        """
        Query genotypes for a variant, trying flipped alleles if original not found.

        Args:
            variant (str): Variant ID in any supported format

        Returns:
            np.ndarray|None: Array of genotypes or None if error/not found
        """
        try:
            std_variant = self._standardize_variant_format(variant)
            if std_variant is None:
                print(f"Variant {variant} could not be standardized!")
                return None

            # Try original variant
            result = self._query_variants_internal(std_variant)
            if result is not None:
                return result

            # Try flipped alleles
            flipped = self._flip_alleles(std_variant)
            if flipped:
                genotypes = self._query_variants_internal(flipped)
                # Because we flipped the alleles, we need to flip the genotypes...
                genotypes[genotypes != -1] = 2 - genotypes[genotypes != -1]
                return genotypes

            print(f"Variant not found in either orientation: {variant}")
            return None

        except Exception as e:
            print(f"Error in query_variants: {str(e)}")
            return None

    def _query_variants_internal(self, variant: str) -> np.ndarray | None:
        """Internal method to query a single variant."""
        try:
            variant_mask = self._single_variant_mask(variant)
            if variant_mask is None:
                return None

            # Find matching variant
            selected_variant_indices = np.where(variant_mask)[0]
            if len(selected_variant_indices) == 0:
                return None

            # Query the submatrix
            try:
                submatrix = self.genotype_matrix[selected_variant_indices, :]
                if submatrix.size == 0:
                    return None
                return submatrix
            except Exception as e:
                print(f"Error accessing genotype matrix: {str(e)}")
                return None

        except Exception as e:
            print(f"Error in _query_variants_internal: {str(e)}")
            return None

    def _single_variant_mask(self, varsel: str) -> np.ndarray:
        """Create mask for single variant."""
        try:
            mask = self.variants == varsel
            if mask.sum() == 0:
                return None
            return mask
        except Exception as e:
            print(f"Error in _single_variant_mask: {str(e)}")
            raise

# ------------------------------------------------------------------------------#
# icd10 = 'E831'
# ------------------------------------------------------------------------------#


class ICD10HDF5:
    def __init__(self, hdf5_file):
        self.hdf5_file = hdf5_file
        self.f = h5py.File(hdf5_file, "r")

        # Load the data matrix and index information
        self.phenotype_data = self.f["phenotype_data"]
        self.row_indices = self.f["row_indices"][...]  # Load row data
        self.col_indices = self.f["col_indices"][...]  # Load column data

    def get_cases_for_icd10(self, icd10):
        # get data for col indice E831
        indices = np.where(self.col_indices.astype(str) == icd10)[0]

        # get people phenotype for col indice E831
        indices_data = self.phenotype_data[:, indices]

        # get the row indices for the data = 1
        row_indices = np.where(indices_data == 1)[
            0
        ]  # for row indices, take the second index
        cases = self.row_indices[row_indices]

        cases = cases.astype(int)
        return cases


# ------------------------------------------------------------------------------#
# Nomaly ScoreHDF5
# ------------------------------------------------------------------------------#


class ScoreHDF5:
    def __init__(self, hdf5_file):
        self.hdf5_file = hdf5_file
        self.f = h5py.File(hdf5_file, "r")

        # Load the data matrix and index information
        self.data_matrix = self.f["scores"]
        self.rows = self.f["eid"][...]  # Load row data
        self.columns = self.f["term"][...]  # Load column data

    def get_score_by_eid(self, eid):
        mask_row = self.rows.astype(str) == eid
        mask_row_indices = np.where(mask_row)[0][0]
        return self.data_matrix[mask_row_indices, :]

    def get_score_by_term(self, term):
        mask_column = self.columns.astype(str) == term
        mask_column_indices = np.where(mask_column)[0][0]
        return self.data_matrix[:, mask_column_indices]


class PhenotypesHDF5:
    def __init__(self, hdf5_file=phenotypes_h5):
        self.hdf5_file = hdf5_file
        self.f = h5py.File(hdf5_file, "r")

        self._population_mask_cache = {}  # Add cache dictionary
        self._population_eids_cache = {}  # Cache for filtered eids

        # Load the data matrix and index information as numpy arrays
        self.phenotype_data = self.f["phenotype_data"][...][...]
        self.eids = self.f['eids'][...]  # Load row data
        self.phecodes = self.f['phecodes'][...]  # Load column data
        self.populations = self.f['populations'][...]  # Load population data

        # Convert to strings
        self.phecodes: np.ndarray = self.phecodes.astype(str)
        self.populations: np.ndarray = self.populations.astype(str)

    def __str__(self):
        return f"PhenotypesHDF5: {self.hdf5_file}"

    @cached_property
    def phecode_to_index(self):
        phecodes = self.phecodes.astype(str)
        return {phecode: index for index, phecode in enumerate(phecodes)}

    def get_population_mask(self, population):
        if population not in self._population_mask_cache:
            # Note: We use 'equal' here for speed.
            # Be aware that there are separate EUR and EUR_S populations!
            self._population_mask_cache[population] = np.char.equal(
                self.populations, population
            )
        return self._population_mask_cache[population]

    def get_population_eids(self, population):
        if population not in self._population_eids_cache:
            population_mask = self.get_population_mask(population)
            self._population_eids_cache[population] = self.eids[population_mask]
        return self._population_eids_cache[population]

    def get_cases_for_phecode(self, phecode, population=None):
        try:
            phecode_index = self.phecode_to_index[phecode]
        except KeyError:
            raise ValueError(f"Phecode '{phecode}' not found in the data matrix!")

        if population is None:
            return self.eids, self.phenotype_data[:, phecode_index]
        else:
            population_mask = self.get_population_mask(population)
            if population_mask.sum() == 0:
                raise ValueError(
                    f"Population '{population}' not found in the data matrix!"
                )
            eids = self.get_population_eids(population)
            cases_data = self.phenotype_data[:, phecode_index][population_mask]
            return eids, cases_data


# ------------------------------------------------------------------------------#
# Initiate classes
# ------------------------------------------------------------------------------#
NOMALY_RESULTS_DIR = '/data/general/UKBB/Run-v1/DatabaseInputs/'
nomaly_stats_h5 = NOMALY_RESULTS_DIR + 'stats.h5'
nomaly_genotype_h5 = NOMALY_RESULTS_DIR + 'genotypes.h5'
nomaly_scores_h5 = NOMALY_RESULTS_DIR + 'float16_scores.h5'

nomaly_stats = StatsHDF5(nomaly_stats_h5)
nomaly_genotype = GenotypeHDF5(nomaly_genotype_h5)
nomaly_scores = ScoreHDF5(nomaly_scores_h5)

NOMALY_RESULTS_DIR_V2 = '/data/general/UKBB/Run-v2/DatabaseInputs/'
nomaly_stats_h5_v2 = NOMALY_RESULTS_DIR_V2 + 'stats.h5'
nomaly_genotype_h5_v2 = NOMALY_RESULTS_DIR_V2 + 'genotypes.h5'
nomaly_scores_h5_v2 = NOMALY_RESULTS_DIR_V2 + 'float16_scores.h5'

nomaly_stats_v2 = StatsHDF5(nomaly_stats_h5_v2)
nomaly_genotype_v2 = GenotypeHDF5(nomaly_genotype_h5_v2)
nomaly_scores_v2 = ScoreHDF5(nomaly_scores_h5_v2)

# Don't think we need this?
# nomaly_phenotypes = PhenotypesHDF5(phenotypes_h5)

icd10_cases = ICD10HDF5(icd10_cases_h5)
# ------------------------------------------------------------------------------#
# plotting functions
# ------------------------------------------------------------------------------#

def qqstats(dfstats):
    tags = dfstats.columns[dfstats.columns.str.endswith('_pvalue')]

    melt_stats = pd.melt(dfstats,
                         id_vars=['term'], value_vars=tags, var_name='tag', value_name='P_obs')
    melt_stats.dropna(inplace=True)
    # add metric column as tag
    melt_stats['test'] = melt_stats['tag'].apply(lambda x: x.split('_pvalue')[0])
    # # drop the tag column
    # melt_stats.drop(columns='tag', inplace=True)
    # add -log10(P_obs) column
    try:
        melt_stats['-log10(observed)'] = -np.log10(melt_stats['P_obs'])
    except Exception as e:
        print(f'Exception "{e}" encourterd for {melt_stats["term"][0]}')
    # sort the table by P_obs
    melt_stats.sort_values('P_obs', inplace=True)
    # add -log10(expected) column
    for tag in tags:
        len_tag = len(melt_stats[melt_stats['tag']==tag])
        melt_stats.loc[melt_stats['tag']==tag, '-log10(expected)'] = -np.log10(np.linspace(0+1/len_tag, 1-1/len_tag, len_tag))

    return melt_stats

def make_qqplot(plot_df):
    try:
        melt_stats = qqstats(plot_df)
    except Exception as e:
        print(f'Exception "{e}" encourterd')
        # save plot_df to a temp file
        plot_df.to_csv('temp.csv', sep='\t', index=None)
        return None
    # Add a scatter plot with plotly
    xlabel = '-log10(expected)'
    ylabel = '-log10(observed)'
    fig = px.scatter(
        melt_stats, x=xlabel, y=ylabel, color='test', 
        hover_name='term', 
        # title=f'{disease_select} QQ plot'
    )
    # add the diagonal line
    lims =[0, melt_stats['-log10(expected)'].max()]
    fig.add_scatter(x=lims, y=lims, mode='lines', 
                    name ='Expected', line=dict(color='gray', dash='dash'))
    
    # figure size
    fig.update_layout(
        width=600, 
        height = 400,
    )

    return fig

# ------------------------------------------------------------------------------#
# Add termnames, term2genes
# ------------------------------------------------------------------------------#

# from db import get_db_connection
