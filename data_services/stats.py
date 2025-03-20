from pathlib import Path
from typing import Optional, Dict, Tuple

import h5py
import numpy as np
import pandas as pd

class StatsRegistry:
    """Registry for managing multiple StatsService instances."""

    def __init__(self, stats_selector: Dict):
        self.stats_selector = stats_selector
        self._services: Dict[Tuple[str, str], StatsService] = {}

    def get(self, run_version: str, ancestry: str) -> "StatsService":
        """Get a StatsService instance for the specified run version and ancestry.

        Args:
            run_version: Version of the run (e.g., "Run-v1")
            ancestry: Ancestry code (e.g., "AFR", "EUR")

        Returns:
            StatsService for the requested data

        Raises:
            ValueError: If the requested combination doesn't exist
        """
        key = (run_version, ancestry)

        if key not in self._services:
            try:
                file_path = self.stats_selector[run_version][ancestry]
                self._services[key] = StatsService(file_path)
            except KeyError:
                raise ValueError(f"No stats file found for {run_version}/{ancestry}")

        return self._services[key]


class StatsService:
    def __init__(self, hdf5_file: Path | str):
        self._hdf = StatsHDF5(hdf5_file)

    # Expose the StatsHDF5 methods here
    def get_stats_by_term_phecode(
        self, term: str, phecode: str, statstype: Optional[str] = None
    ) -> np.ndarray:
        return self._hdf.get_stats_by_term_phecode(term, phecode, statstype)


class StatsHDF5:
    """Handles access to stats data stored in HDF5 format."""

    def __init__(self, hdf5_file):
        self.f = h5py.File(hdf5_file, "r")

        # TODO: Fix the naming of the fields
        data = self.f["data"]
        term = self.f["term"]
        phecode = self.f["phecode"]
        stats_type = self.f["stats_type"]

        # Jumping through hoops to fix the type checker...
        assert isinstance(data, h5py.Dataset)
        assert isinstance(term, h5py.Dataset)
        assert isinstance(phecode, h5py.Dataset)
        assert isinstance(stats_type, h5py.Dataset)

        # Sanity checks
        # TODO: Move these to integration tests?
        try:
            print(f"RUNNING SANITY CHECK ON {hdf5_file}! ONLY SEE THIS ONCE!")
            assert data.shape[0] == term.shape[0]
            assert data.shape[1] == phecode.shape[0]
            assert data.shape[2] == stats_type.shape[0]
        except Exception as e:
            print(f"Error in sanity checks: {str(e)}")
            raise

        # Load the data matrix and index information as numpy arrays
        self.data: h5py.Dataset = data
        self.term: np.ndarray = term[...].astype(str)
        self.phecode: np.ndarray = phecode[...].astype(str)
        self.stats_type: np.ndarray = stats_type[...].astype(str)

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
            phecode_mask = np.ones(len(self.phecode), dtype=bool)
        else:
            phecode_mask = np.isin(self.phecode, phecode)

        if term is None:
            term_mask = np.ones(len(self.term), dtype=bool)
        else:
            term_mask = np.isin(self.term, term)

        if statstype is None:
            statstype_mask = np.ones(len(self.stats_type), dtype=bool)
        else:
            statstype_mask = np.isin(self.stats_type, statstype)

        # NOTE: All data is stored as 'float' in the HDF, however, some columns
        # are actually integers. TODO: This is fine, but we could convert to int
        # here.
        return pd.DataFrame(
            self.data[term_mask][:, phecode_mask][:, :, statstype_mask],
            index=self.term[term_mask],
            columns=self.phecode[phecode_mask],
            dtype=float,
        )

    def get_stats_by_phecode(
        self, phecode: str, statstype: None | str | list = None
    ) -> pd.DataFrame | np.ndarray:
        if phecode not in self.phecode:
            raise ValueError(f"Disease {phecode} not found in dataset")

        mask_column = self.phecode == phecode
        mask_column_indices = np.where(mask_column)[0][0]

        if statstype is None:
            # Return all statstypes for the disease
            disease_data = self.data[:, mask_column_indices, :]
            return pd.DataFrame(
                disease_data, index=self.term.astype(str), columns=self.stats_type
            )

        # NOTE: The return types of these alternative calls are different!
        elif isinstance(statstype, str):
            # Return the data for the given statstype
            mask_3rd_dim = self.stats_type == statstype
            mask_3rd_dim_indices = np.where(mask_3rd_dim)[0][0]
            return self.data[:, mask_column_indices, mask_3rd_dim_indices]

        elif isinstance(statstype, list):
            # Return the data for the given statstypes
            mask_3rd_dim = np.isin(self.stats_type, statstype)
            mask_3rd_dim_indices = np.where(mask_3rd_dim)[0]
            return self.data[:, mask_column_indices, mask_3rd_dim_indices]

    # TODO: Perhaps allow phecode to also be a list like above?
    # TODO: Think about making return type consistent.
    def get_stats_by_term_phecode(
        self, term: str, phecode: str, statstype: None | str | np.ndarray = None
    ) -> dict[str, np.ndarray]:
        if statstype is None:
            statstype = self.stats_type
        elif isinstance(statstype, str):
            statstype = np.array([statstype])
        elif isinstance(statstype, np.ndarray):
            pass

        term_mask = self.term == term
        mask_row_indices = np.where(term_mask)[0][0]

        phecode_mask = self.phecode == phecode
        mask_column_indices = np.where(phecode_mask)[0][0]

        statstype_mask = np.isin(self.stats_type, statstype)
        statstype_mask_indices = np.where(statstype_mask)[0]

        # return a dictionary with the statstype as key and the value as the data
        statsdict = {
            self.stats_type[i]: self.data[mask_row_indices, mask_column_indices, i]
            for i in statstype_mask_indices
        }

        return statsdict
