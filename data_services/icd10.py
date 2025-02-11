# type: ignore

import h5py
import numpy as np


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
