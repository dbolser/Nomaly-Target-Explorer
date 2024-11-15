import os
import h5py
import numpy as np
import pandas as pd

import plotly.express as px

UKBB_PHENO_DIR = '/data/general/UKBB/Phenotypes/'
icd10_cases_h5 = UKBB_PHENO_DIR + 'ukbbrun_icd10_2024-09-07_any.h5'

NOMALY_RESULTS_DIR = '/data/general/UKBB/Run-v1/DatabaseInputs/'
nomaly_stats_h5 = NOMALY_RESULTS_DIR + 'stats.h5'
nomaly_genotype_h5 = NOMALY_RESULTS_DIR + 'genotypes.h5'
nomaly_scores_h5 = NOMALY_RESULTS_DIR + 'float16_scores.h5'
# ------------------------------------------------------------------------------#
# Nomaly stats pre-calculated
# ------------------------------------------------------------------------------#

class StatsHDF5:
    def __init__(self, hdf5_file):
        self.hdf5_file = hdf5_file
        self.f = h5py.File(hdf5_file, 'r')

        # Load the data matrix and index information
        self.data_matrix = self.f['data']
        self.rows = self.f['rows_term'][...].astype(str)  # Load term data
        self.columns = self.f['columns_phecode'][...].astype(str)  # Load phecode data
        self.third_dim = self.f['3rddim_statstype'][...].astype(str)  # Load statstype data

    def get_stats_by_disease(self, disease, statstype=None):
        if disease not in self.columns:
            disease = disease + '.0'
            if disease not in self.columns:
                raise ValueError(f"{disease} not found in the columns of the data matrix.")

        mask_column = self.columns == disease
        mask_column_indices = np.where(mask_column)[0][0]
        if statstype is None:
            # return a dataframe
            disease_data = self.data_matrix[:, mask_column_indices, :]
            disease_DF = pd.DataFrame(
                disease_data, index=self.rows.astype(str), columns=self.third_dim
            )
            return disease_DF
        elif isinstance(statstype, str) and statstype in self.third_dim:
            mask_3rd_dim = self.third_dim == statstype
            mask_3rd_dim_indices = np.where(mask_3rd_dim)[0][0]
            return self.data_matrix[:,mask_column_indices, mask_3rd_dim_indices]
        elif isinstance(statstype, list):
            mask_3rd_dim = np.isin(self.third_dim, statstype)
            mask_3rd_dim_indices = np.where(mask_3rd_dim)[0]
            return self.data_matrix[:,mask_column_indices, mask_3rd_dim_indices]
        else:
            raise ValueError(f"statstype must be a string or a list of strings. statstype: {statstype}")


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
        self.hdf5_file = hdf5_file
        self.f = h5py.File(hdf5_file, 'r')

        # Load the genotype matrix and index information
        self.genotype_matrix = self.f['genotype_matrix']
        self.individual = self.f['fam'][...]  # Load FAM data (individual information)
        self.variants = self.f['bim'][...]  # Load BIM data (variant information)

        # 6:42963149:A:C 0 is AA, 1 is AC, 2 is CC, -1 is missing
    
    def get_genotype(self, i, j):
        return self.genotype_matrix[i, j]
    
    def query_variants(self, variantsel) -> np.ndarray:
        if isinstance(variantsel, str):
            variant_mask = self._single_variant_mask(variantsel)
        elif isinstance(variantsel, list):
            variant_mask = self._multiple_variant_mask(variantsel)
        else:
            raise ValueError("variantsel must be a string or a list of strings.")

        # Apply mask to Find indices of selected variants and individuals
        selected_variant_indices = np.where(variant_mask)[0]

        if len(selected_variant_indices) == 0:
            return None

        # Query the submatrix using selected indices
        submatrix = self.genotype_matrix[selected_variant_indices, :]

        return submatrix
        
    def _single_variant_mask(self, varsel = '6:42963149:A:C'):
        return self.variants.astype(str) == varsel

    def _multiple_variant_mask(self, varsel = ['6:42963149:A:C', '1:930245:A:G']):
        return np.isin(self.variants.astype(str), varsel)
    
# ------------------------------------------------------------------------------#
# icd10 = 'E831'
# ------------------------------------------------------------------------------#

class ICD10HDF5:
    def __init__(self, hdf5_file):
        self.hdf5_file = hdf5_file
        self.f = h5py.File(hdf5_file, 'r')

        # Load the data matrix and index information
        self.phenotype_data = self.f['phenotype_data']
        self.row_indices = self.f['row_indices'][...]  # Load row data
        self.col_indices = self.f['col_indices'][...]  # Load column data

    def get_cases_for_icd10(self, icd10):
        # get data for col indice E831
        indices = np.where(self.col_indices.astype(str) == icd10)[0]

        # get people phenotype for col indice E831
        indices_data = self.phenotype_data[:, indices]

        # get the row indices for the data = 1
        row_indices = np.where(indices_data == 1)[0] # for row indices, take the second index
        cases = self.row_indices[row_indices]

        cases = cases.astype(int)
        return cases

# ------------------------------------------------------------------------------#
# Nomaly ScoreHDF5
# ------------------------------------------------------------------------------#

class ScoreHDF5:
    def __init__(self, hdf5_file):
        self.hdf5_file = hdf5_file
        self.f = h5py.File(hdf5_file, 'r')

        # Load the data matrix and index information
        self.data_matrix = self.f['scores']
        self.rows = self.f['eid'][...]  # Load row data
        self.columns = self.f['term'][...]  # Load column data

    def get_score_by_eid(self, eid):
        mask_row = self.rows.astype(str) == eid
        mask_row_indices = np.where(mask_row)[0][0]
        return self.data_matrix[mask_row_indices, :]

    def get_score_by_term(self, term):
        mask_column = self.columns.astype(str) == term
        mask_column_indices = np.where(mask_column)[0][0]
        return self.data_matrix[:, mask_column_indices]

# ------------------------------------------------------------------------------#
# Initiate classes
# ------------------------------------------------------------------------------#
nomaly_stats = StatsHDF5(nomaly_stats_h5)
nomaly_genotype = GenotypeHDF5(nomaly_genotype_h5)
nomaly_scores = ScoreHDF5(nomaly_scores_h5)

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
