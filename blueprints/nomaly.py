import os
import re
import traceback

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import plotly.express as px
from flask import current_app

from config import Config


# UKBB_PHENO_DIR = "/data/general/UKBB/Phenotypes/"
# icd10_cases_h5 = UKBB_PHENO_DIR + "ukbbrun_icd10_2024-09-07_any.h5"
# phenotypes_h5 = Config.PHENOTYPES_H5

RESOURCE_DATA_DIR = Config.RESOURCE_DATA_DIR
pharos_path = RESOURCE_DATA_DIR / "pharos_api_query.out"
pp_path = RESOURCE_DATA_DIR / "pp_by_gene.tsv"

pharos = pd.read_csv(pharos_path, sep="\t", encoding="ISO-8859-1").rename(
    columns={"#symbol": "gene"}
)

pp = pd.read_csv(pp_path, sep="\t", index_col=0).rename(
    columns={"Description": "drug_program_indication"}
)
pp["gene"] = pp.index


class StatsHDF5:
    def __init__(self, hdf5_file):
        try:
            self.hdf5_file = hdf5_file
            if not os.path.exists(hdf5_file):
                raise FileNotFoundError(f"HDF5 file not found: {hdf5_file}")
            self.f = h5py.File(hdf5_file, "r")

            # Load the data matrix and index information
            self.data_matrix: h5py.Dataset = self.f["data"]
            self.rows: np.ndarray = self.f["rows_term"][...].astype(
                str
            )  # Load term data
            self.columns: np.ndarray = self.f["columns_phecode"][...].astype(
                str
            )  # Load phecode data
            self.third_dim: np.ndarray = self.f["3rddim_statstype"][...].astype(
                str
            )  # Load statstype data

            # rename third_dim that ends with _p_value to _pvalue
            self.third_dim = np.array(
                [x.replace("_p_value", "_pvalue") for x in self.third_dim]
            )
        except Exception as e:
            current_app.logger.error(
                f"StatsHDF5 initialization error: {e}\n{traceback.format_exc()}"
            )
            raise

    def get_stats_by_disease(self, disease, statstype=None):
        try:
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
                return self.data_matrix[:, mask_column_indices, mask_3rd_dim_indices]
            elif isinstance(statstype, list):
                # if statstype is a list of strings, return the data for the statstypes
                mask_3rd_dim = np.isin(self.third_dim, statstype)
                mask_3rd_dim_indices = np.where(mask_3rd_dim)[0]
                return self.data_matrix[:, mask_column_indices, mask_3rd_dim_indices]
            else:
                raise ValueError(
                    f"statstype must be a string or a list of strings. statstype: {statstype}"
                )
        except Exception as e:
            current_app.logger.error(f"Error getting stats for disease {disease}: {e}")
            raise

    def get_stats_by_term_disease(self, term, disease):
        if disease not in self.columns:
            disease = disease + ".0"
            if disease not in self.columns:
                raise ValueError(
                    f"{disease} not found in the columns of the data matrix."
                )

        mask_column = self.columns == disease
        mask_column_indices = np.where(mask_column)[0][0]

        if term not in self.rows:
            raise ValueError(f"{term} not found in the rows of the data matrix.")

        mask_row = self.rows == term
        mask_row_indices = np.where(mask_row)[0][0]

        # return a dictionary with the statstype as key and the value as the data
        statsdict = {
            self.third_dim[i]: self.data_matrix[
                mask_row_indices, mask_column_indices, i
            ]
            for i in range(len(self.third_dim))
        }
        return statsdict


class ICD10HDF5:
    def __init__(self, hdf5_file):
        self.hdf5_file = hdf5_file
        self.f = h5py.File(hdf5_file, "r")

        # Load the data matrix and index information
        self.phenotype_data: np.ndarray = self.f["phenotype_data"]
        self.row_indices: np.ndarray = self.f["row_indices"][...]  # Load row data
        self.col_indices: np.ndarray = self.f["col_indices"][...]  # Load column data

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


class ScoreHDF5:
    def __init__(self, hdf5_file):
        self.hdf5_file = hdf5_file
        self.f = h5py.File(hdf5_file, "r")

        # Load the data matrix and index information
        self.data_matrix: np.ndarray = self.f["scores"]
        self.rows: np.ndarray = self.f["eid"][...]  # Load row data
        self.columns: np.ndarray = self.f["term"][...]  # Load column data

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

# Trying to not use globals...

# nomaly_genotype = GenotypesHDF5(Config.GENOTYPES_H5)

# nomaly_stats = StatsHDF5(Config.STATS_H5)
# nomaly_stats_v2 = StatsHDF5(Config.STATS_H5_V2)

# NOT USED!
# nomaly_scores = ScoreHDF5(Config.SCORES_H5)
# nomaly_scores_v2 = ScoreHDF5(Config.SCORES_H5_V2)

# icd10_cases = ICD10HDF5(icd10_cases_h5)


# ------------------------------------------------------------------------------#
# plotting functions
# ------------------------------------------------------------------------------#


def qqstats(dfstats):
    tags = dfstats.columns[dfstats.columns.str.endswith("_pvalue")]

    melt_stats = pd.melt(
        dfstats, id_vars=["term"], value_vars=tags, var_name="tag", value_name="P_obs"
    )
    melt_stats.dropna(inplace=True)
    # add metric column as tag
    melt_stats["test"] = melt_stats["tag"].apply(lambda x: x.split("_pvalue")[0])
    # # drop the tag column
    # melt_stats.drop(columns='tag', inplace=True)
    # add -log10(P_obs) column
    try:
        melt_stats["-log10(observed)"] = -np.log10(melt_stats["P_obs"])
    except Exception as e:
        print(f'Exception "{e}" encourterd for {melt_stats["term"][0]}')
    # sort the table by P_obs
    melt_stats.sort_values("P_obs", inplace=True)
    # add -log10(expected) column
    for tag in tags:
        len_tag = len(melt_stats[melt_stats["tag"] == tag])
        melt_stats.loc[melt_stats["tag"] == tag, "-log10(expected)"] = -np.log10(
            np.linspace(0 + 1 / len_tag, 1 - 1 / len_tag, len_tag)
        )

    return melt_stats


def make_qqplot(plot_df):
    try:
        melt_stats = qqstats(plot_df)
    except Exception as e:
        print(f'Exception "{e}" encourterd')
        # save plot_df to a temp file
        plot_df.to_csv("temp.csv", sep="\t", index=None)
        return None
    # Add a scatter plot with plotly
    xlabel = "-log10(expected)"
    ylabel = "-log10(observed)"
    fig = px.scatter(
        melt_stats,
        x=xlabel,
        y=ylabel,
        color="test",
        hover_name="term",
        # title=f'{disease_select} QQ plot'
    )
    # add the diagonal line
    lims = [0, melt_stats["-log10(expected)"].max()]
    fig.add_scatter(
        x=lims,
        y=lims,
        mode="lines",
        name="Expected",
        line=dict(color="gray", dash="dash"),
    )

    # figure size
    fig.update_layout(
        width=600,
        height=400,
    )

    return fig
