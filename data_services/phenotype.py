from pathlib import Path

import h5py
import numpy as np

from typing import Optional

# Create a 'dummy' profile decorator if we don't have line_profiler installed
try:
    from line_profiler import profile
except ImportError:

    def profile(func):
        return func


class PhenotypeService:
    def __init__(self, hdf5_file: Path | str):
        self._hdf = PhenotypesHDF5(hdf5_file)  # Existing implementation

    def get_cases_for_phecode(
        self, phecode, population: str | None = None, biological_sex: str | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Delegate to the underlying HDF5 file's get_cases_for_phecode method."""
        return self._hdf.get_cases_for_phecode(phecode, population, biological_sex)

    @property
    def phenotype_data(self):
        return self._hdf.phenotype_data

    @property
    def eids(self):
        return self._hdf.eids

    @property
    def phecodes(self):
        return self._hdf.phecodes

    def get_phenotypes(
        self, eids: Optional[np.ndarray] = None, phecodes: Optional[np.ndarray] = None
    ) -> np.ndarray:
        return self._hdf.get_phenotypes(eids, phecodes)

class PhenotypesHDF5:
    """Handles access to phenotype data stored in HDF5 format."""

    def __init__(self, hdf5_file: Path | str):
        self.f = h5py.File(hdf5_file, "r")

        # Jumping through hoops to fix the type checker...
        # TODO: Fix the naming of the fields
        phenotype_data = self.f["phenotype_data"]
        eids = self.f["eids"]
        phecodes = self.f["phecodes"]
        populations = self.f["populations"]
        affected_sex = self.f["phecode_sex"]
        biological_sex = self.f["affected_sex"]

        # Type checker magic
        assert isinstance(phenotype_data, h5py.Dataset)
        assert isinstance(eids, h5py.Dataset)
        assert isinstance(phecodes, h5py.Dataset)
        assert isinstance(populations, h5py.Dataset)
        assert isinstance(affected_sex, h5py.Dataset)
        assert isinstance(biological_sex, h5py.Dataset)

        # Sanity checks
        # TODO: Move these to integration tests?
        try:
            assert phenotype_data.shape[0] == eids.shape[0]
            assert biological_sex.shape[0] == eids.shape[0]
            assert populations.shape[0] == eids.shape[0]

            assert phenotype_data.shape[1] == phecodes.shape[0]
            assert affected_sex.shape[0] == phecodes.shape[0]

            assert np.all(np.diff(eids) > 0), "EIDs are expected to be sorted"

        except Exception as e:
            print(f"Error in sanity checks: {str(e)}")
            raise

        # Load the data matrix and index information as numpy arrays
        # I think this may be a massive performance hit!
        # self.phenotype_data: h5py.Dataset = phenotype_data
        self.phenotype_data: np.ndarray = phenotype_data[...].astype(int)
        self.eids: np.ndarray = eids[...].astype(int)
        self.populations: np.ndarray = populations[...].astype(str)
        self.biological_sex: np.ndarray = biological_sex[...].astype(str)

        self.phecodes: np.ndarray = phecodes[...].astype(str)
        self.affected_sex: np.ndarray = affected_sex[...].astype(str)

        # TODO: What about eid to index mapping?
        # TODO: Is this useful?
        self.phecode_to_index = {
            phecode: index for index, phecode in enumerate(self.phecodes)
        }

        # TODO: Do we need a memory-mapped matrix here too? (See GenotypeService)

        # Precompute sorted indices for faster lookups
        self._phecodes_argsort = np.argsort(self.phecodes)
        self._phecodes_sorted = self.phecodes[self._phecodes_argsort]

    def get_phenotypes(
        self, eids: Optional[np.ndarray] = None, phecodes: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Get phenotypes for a list of eids and phecodes.

        Args:
            eids: List of eids to get phenotypes for
            phecodes: List of phecodes to get phenotypes for

        Returns:
            np.ndarray: Array of phenotypes with shape (n_eids, n_phecodes) in the order of the
                provided eids and phecodes.

            Phenotype values are encoded as:
                - 0 = control
                - 1 = affected
                - 9 = excluded
                - 8 = additional sex based exclusion... TBD
        """
        if eids is None and phecodes is None:
            return self.phenotype_data

        if eids is not None:
            # Vectorized index lookup using precomputed sorted arrays
            eid_pos = np.searchsorted(self.eids, eids)
            eid_idx = np.argsort(eid_pos)

        if phecodes is not None:
            phecode_pos = np.searchsorted(self._phecodes_sorted, phecodes)
            phecode_idx = self._phecodes_argsort[phecode_pos]

        if eids is None:
            return self.phenotype_data[phecode_idx, :]

        if phecodes is None:
            return self.phenotype_data[eid_idx, :]

        return self.phenotype_data[eid_idx, :][:, phecode_idx]

    def get_population_mask(self, population):
        return self.populations == population

    def get_biological_sex_mask(self, biological_sex):
        return self.biological_sex == biological_sex

    def get_affected_sex_mask(self, affected_sex):
        return self.affected_sex == affected_sex

    @profile
    def get_cases_for_phecode(
        self, phecode, population: str | None = None, biological_sex: str | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get the cases for a given phecode, population, and biological sex."""
        phecode_index = self.phecode_to_index[phecode]

        if biological_sex is not None:
            biological_sex_mask = self.get_biological_sex_mask(biological_sex)
        else:
            biological_sex_mask = np.ones(self.biological_sex.shape, dtype=bool)

        if population is not None:
            population_mask = self.get_population_mask(population)
        else:
            population_mask = np.ones(self.populations.shape, dtype=bool)

        return (
            self.eids[population_mask & biological_sex_mask],
            self.phenotype_data[:, phecode_index][
                population_mask & biological_sex_mask
            ],
        )

    def get_case_counts_for_phecode(
        self, phecode, population: str | None = None, biological_sex: str | None = None
    ):
        eids, cases = self.get_cases_for_phecode(
            phecode, population=population, biological_sex=biological_sex
        )
        values, counts = np.unique(cases, return_counts=True)
        counts_dict = dict(zip(values, counts))

        counts_dict["affected"] = counts_dict.get(1, 0)
        counts_dict["excluded"] = counts_dict.get(9, 0)
        counts_dict["control"] = counts_dict.get(0, 0)

        return counts_dict
