from pathlib import Path

import h5py
import numpy as np
import pandas as pd


class StatsService:
    def __init__(self, hdf5_file: Path | str):
        self._hdf = StatsHDF5(hdf5_file)  # Existing implementation


class StatsHDF5:
    """Handles access to stats data stored in HDF5 format."""

    def __init__(self, hdf5_file):
        self.f = h5py.File(hdf5_file, "r")

        # Jumping through hoops to fix the type checker...
        # TODO: Fix the naming of the fields
        data = self.f["data"]
        terms = self.f["rows_term"]
        phecodes = self.f["columns_phecode"]
        statistics = self.f["3rddim_statstype"]

        # Sanity checks
        assert isinstance(data, h5py.Dataset)
        assert isinstance(terms, h5py.Dataset)
        assert isinstance(phecodes, h5py.Dataset)
        assert isinstance(statistics, h5py.Dataset)

        # Sanity checks
        # TODO: Move these to integration tests?
        try:
            assert data.shape[0] == terms.shape[0]
            assert data.shape[1] == phecodes.shape[0]
            assert data.shape[2] == statistics.shape[0]
        except Exception as e:
            print(f"Error in sanity checks: {str(e)}")
            raise

        # Load the data matrix and index information as numpy arrays
        self.data: h5py.Dataset = data
        self.terms: np.ndarray = terms[...].astype(str)
        self.phecodes: np.ndarray = phecodes[...].astype(str)
        self.statistics: np.ndarray = statistics[...].astype(str)

    def get_stats_by_phecode(
        self, phecode: str, statstype: None | str | list = None
    ) -> pd.DataFrame | np.ndarray:
        if phecode not in self.phecodes:
            raise ValueError(f"Disease {phecode} not found in dataset")

        mask_column = self.phecodes == phecode
        mask_column_indices = np.where(mask_column)[0][0]

        if statstype is None:
            # Return all statstypes for the disease
            disease_data = self.data[:, mask_column_indices, :]
            return pd.DataFrame(
                disease_data, index=self.terms.astype(str), columns=self.statistics
            )

        elif isinstance(statstype, str):
            # Return the data for the given statstype
            mask_3rd_dim = self.statistics == statstype
            mask_3rd_dim_indices = np.where(mask_3rd_dim)[0][0]
            return self.data[:, mask_column_indices, mask_3rd_dim_indices]

        elif isinstance(statstype, list):
            # Return the data for the given statstypes
            mask_3rd_dim = np.isin(self.statistics, statstype)
            mask_3rd_dim_indices = np.where(mask_3rd_dim)[0]
            return self.data[:, mask_column_indices, mask_3rd_dim_indices]

    def get_stats_by_term_phecode(
        self, term: str, phecode: str
    ) -> dict[str, np.ndarray]:
        mask_column = self.phecodes == phecode
        mask_column_indices = np.where(mask_column)[0][0]

        mask_row = self.terms == term
        mask_row_indices = np.where(mask_row)[0][0]

        # return a dictionary with the statstype as key and the value as the data
        statsdict = {
            self.statistics[i]: self.data[mask_row_indices, mask_column_indices, i]
            for i in range(len(self.statistics))
        }

        return statsdict
