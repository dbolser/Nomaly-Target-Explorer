from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import pandas as pd


class StatsService:
    def __init__(self, hdf5_file: Path | str):
        self._hdf = StatsHDF5(hdf5_file)  # Existing implementation

    def get_stats_by_term_phecode(
        self, term: str, phecode: str, statstype: Optional[str] = None
    ) -> np.ndarray:
        return self._hdf.get_stats_by_term_phecode(term, phecode, statstype)


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

    # TODO: ALL THESE FUNCTIONS SHOULD BE ONE GENERAL FUNTION WITH 3 OPTIONAL ARGUMENTS
    # TODO: 1. phecode (str|list), 2. term (str|list), 3. statstype (str|list)
    # TODO: Return type should be consistent.

    def get_stats_wip(
        self,
        phecode: None | str | list = None,
        term: None | str | list = None,
        statstype: None | str | list = None,
    ) -> pd.DataFrame:
        if isinstance(phecode, str):
            phecode = [phecode]
        if isinstance(term, str):
            term = [term]
        if isinstance(statstype, str):
            statstype = [statstype]

        if not phecode and not term:
            raise ValueError(
                "What the fuck are you doing? At least one of phecode or term must be provided!"
            )
        elif phecode and term:
            if not len(phecode) == 1 or not len(term) == 1:
                raise ValueError(
                    "What the fuck are you doing? Either phecode or term must be a single value!"
                )

        if phecode is None:
            phecode_mask = np.ones(len(self.phecodes), dtype=bool)
        else:
            phecode_mask = np.isin(self.phecodes, phecode)

        if term is None:
            term_mask = np.ones(len(self.terms), dtype=bool)
        else:
            term_mask = np.isin(self.terms, term)

        if statstype is None:
            statstype_mask = np.ones(len(self.statistics), dtype=bool)
        else:
            statstype_mask = np.isin(self.statistics, statstype)

        # NOTE: All data is stored as 'float' in the HDF, however, some columns
        # are actually integers. TODO: This is fine, but we could convert to int
        # here.
        return pd.DataFrame(
            self.data[term_mask][:, phecode_mask][:, :, statstype_mask],
            index=self.terms[term_mask],
            columns=self.phecodes[phecode_mask],
            dtype=float,
        )

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

        # NOTE: The return types of these alternative calls are different!
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

    # TODO: Perhaps allow phecode to also be a list like above?
    # TODO: Think about making return type consistent.
    def get_stats_by_term_phecode(
        self, term: str, phecode: str, statstype: None | str | np.ndarray = None
    ) -> dict[str, np.ndarray]:
        if statstype is None:
            statstype = self.statistics
        elif isinstance(statstype, str):
            statstype = np.array([statstype])
        elif isinstance(statstype, np.ndarray):
            pass

        term_mask = self.terms == term
        mask_row_indices = np.where(term_mask)[0][0]

        phecode_mask = self.phecodes == phecode
        mask_column_indices = np.where(phecode_mask)[0][0]

        statstype_mask = np.isin(self.statistics, statstype)
        statstype_mask_indices = np.where(statstype_mask)[0]

        # return a dictionary with the statstype as key and the value as the data
        statsdict = {
            self.statistics[i]: self.data[mask_row_indices, mask_column_indices, i]
            for i in statstype_mask_indices
        }

        return statsdict
