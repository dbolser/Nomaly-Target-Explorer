from pathlib import Path
from typing import Dict, List, Optional, Tuple, Literal

import h5py
import numpy as np
import pandas as pd


class StatsRegistry:
    """Registry for managing multiple StatsService instances."""

    def __init__(self, stats_selector: Dict | None = None):
        self.stats_selector = stats_selector
        self.initialized = stats_selector is not None
        self._services: Dict[Tuple[str, str], StatsService] = {}

    def _check_initialized(self):
        if not self.initialized:
            raise ValueError("Service not properly initialized: missing filename")

    def get(self, run_version: str = "Run-v1", ancestry: str = "EUR") -> "StatsService":
        """Get a StatsService instance for the specified run version and ancestry.

        Args:
            run_version: Version of the run (e.g., "Run-v1")
            ancestry: Ancestry code (e.g., "AFR", "EUR")

        Returns:
            StatsService for the requested data

        Raises:
            ValueError: If the requested combination doesn't exist
        """
        self._check_initialized()

        key = (run_version, ancestry)

        if key not in self._services:
            try:
                assert self.stats_selector is not None
                file_path = self.stats_selector[run_version][ancestry]
                self._services[key] = StatsService(file_path)
            except KeyError:
                raise ValueError(f"No stats file found for {run_version}/{ancestry}")

        return self._services[key]


class StatsService:
    def __init__(self, hdf5_file: Path | str):
        self._hdf = StatsHDF5(hdf5_file)

    def get_term_stats(
        self,
        term: str,
        phecode: str | None = None,
        phecodes: List[str] | None = None,
        stats_types: List[str] | None = None,
    ) -> pd.DataFrame:
        """Get data for all terms for a specific phecode.

        Args:
            term: The term to get data for
            phecode: The phecode to get data for OR
            phecodes: One or more phecodes to include, or None for all
            stats_types: One or more stats types to include, or None for all

        Returns:
            DataFrame for the term with phecodes as index, stats_types as columns
        """
        if phecode and phecodes:
            raise ValueError("Don't be silly.")
        phecodes = [phecode] if phecode is not None else phecodes
        return self._hdf.get_data_slice(
            t=[term],
            p=phecodes,
            s=stats_types,
            invoker="get_term_stats",
        )

    def get_phecode_stats(
        self,
        phecode: str,
        term: str | None = None,
        terms: List[str] | None = None,
        stats_types: List[str] | None = None,
    ) -> pd.DataFrame:
        """Get data for all phecodes for a specific term.

        Args:
            phecode: The phecode to get data for
            term: The term to get data for
            terms: One or more terms to include, or None for all
            stats_types: One or more stats types to include, or None for all

        Returns:
            DataFrame with phecodes as index, stats_types as columns
        """
        if term and terms:
            raise ValueError("Don't be silly.")
        terms = [term] if term is not None else terms
        return self._hdf.get_data_slice(
            t=terms,
            p=[phecode],
            s=stats_types,
            invoker="get_phecode_stats",
        )


class StatsHDF5:
    """Handles access to stats data stored in HDF5 format."""

    def __init__(self, hdf5_file):
        self.f = h5py.File(hdf5_file, "r")

        data = self.f["data"]
        terms = self.f["term"]  # TODO: Fix the naming of the fields
        phecodes = self.f["phecode"]  # TODO: Fix the naming of the fields
        stats_types = self.f["stats_type"]  # TODO: Fix the naming of the fields

        # Jumping through hoops to fix the type checker...
        assert isinstance(data, h5py.Dataset)
        assert isinstance(terms, h5py.Dataset)
        assert isinstance(phecodes, h5py.Dataset)
        assert isinstance(stats_types, h5py.Dataset)

        # Sanity checks
        # TODO: Move these to integration tests?
        try:
            print(f"RUNNING SANITY CHECK ON {self.f.filename}! ONLY SEE THIS ONCE!")
            assert data.shape[0] == terms.shape[0]
            assert data.shape[1] == phecodes.shape[0]
            assert data.shape[2] == stats_types.shape[0]
        except Exception as e:
            print(f"Error in sanity checks: {str(e)}")
            raise

        # Load the data matrix and index information as numpy arrays
        self.data: h5py.Dataset = data
        self.terms: np.ndarray = terms[...].astype(str)
        self.phecodes: np.ndarray = phecodes[...].astype(str)
        self.stats_types: np.ndarray = stats_types[...].astype(str)

    def get_data_slice(
        self,
        t: Optional[List[str]] = None,
        p: Optional[List[str]] = None,
        s: Optional[List[str]] = None,
        invoker: Literal["get_term_stats", "get_phecode_stats"] | None = None,
    ) -> pd.DataFrame:
        """Get a slice of the data based on the provided filters.

        Args:
            terms: One or more terms to select, or None for all
            phecodes: One or more phecodes to select, or None for all
            stats_types: One or more stats types to select, or None for all

        Returns:
            DataFrame with the requested data slice

        Raises:
            ValueError: If neither terms nor phecodes are provided
        """

        # At least one of terms or phecodes must be provided
        if t is None and p is None:
            raise ValueError("At least one of terms or phecodes must be provided")

        # Create masks for each dimension
        t_msk = (
            np.ones(len(self.terms), dtype=bool)
            if t is None
            else np.isin(self.terms, t)
        )
        p_msk = (
            np.ones(len(self.phecodes), dtype=bool)
            if p is None
            else np.isin(self.phecodes, p)
        )
        s_msk = (
            np.ones(len(self.stats_types), dtype=bool)
            if s is None
            else np.isin(self.stats_types, s)
        )

        # Get the filtered indices
        t_idx = np.where(t_msk)[0]
        p_idx = np.where(p_msk)[0]
        s_idx = np.where(s_msk)[0]

        if len(t_idx) == 0:
            raise ValueError(
                f"No data found for the provided terms in {self.f.filename}\n"
                f"Terms provided: {t}\n"
                f"Terms in dataset: {self.terms}\n"
            )
        if len(p_idx) == 0:
            raise ValueError(
                f"No data found for the provided phecodes in {self.f.filename}\n"
                f"Phecodes provided: {p}\n"
                f"Phecodes in dataset: {self.phecodes}\n"
            )
        if len(s_idx) == 0:
            raise ValueError(
                f"No data found for the provided stats types in {self.f.filename}\n"
                f"Stats types provided: {s}\n"
                f"Stats types in dataset: {self.stats_types}\n"
            )

        # TERM PHECODE / PHECODE TERM
        if len(t_idx) == 1 and len(p_idx) == 1:
            # Only s_idx is a 'fancy' index
            results = self.data[t_idx[0], p_idx[0], s_idx]
            columns = self.stats_types[s_idx]
            if invoker is None:
                index = None
            elif invoker == "get_phecode_stats":
                index = self.terms[t_idx]
            elif invoker == "get_term_stats":
                index = self.phecodes[p_idx]

            return pd.DataFrame(results[np.newaxis, :], index=index, columns=columns)

        # TERM PHECODES
        if len(t_idx) == 1 and len(p_idx) > 1:
            # Trying to do this with advanced indexing on 3 dimensions at once...
            # Do two separate 'fancy indexes' instead.
            results = self.data[t_idx[0], p_idx, :][:, s_idx]

            index = self.stats_types[s_idx]
            columns = self.phecodes[p_idx]

            return pd.DataFrame(results.T, index=index, columns=columns)

        # TERMS PHECODE
        if len(t_idx) > 1 and len(p_idx) == 1:
            # Trying to do this with advanced indexing on 3 dimensions at once...
            # Do two separate 'fancy indexes' instead.
            results = self.data[t_idx, p_idx[0], :][:, s_idx]

            index = self.terms[t_idx]
            columns = self.stats_types[s_idx]

            return pd.DataFrame(results, index=index, columns=columns)

        # TODO: Just need to decide on format of the return value.
        raise ValueError("Multiple terms and phecodes not currently supported")
