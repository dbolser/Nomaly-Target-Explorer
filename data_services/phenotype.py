from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# # Create a 'dummy' profile decorator if we don't have line_profiler installed
# try:
#     from line_profiler import profile  # type: ignore
# except ImportError:

#     def profile(func):
#         return func


class PhenotypeService:
    """Service for getting phenotype data."""

    def __init__(self, hdf5_file: Path | None = None):
        self.hdf5_file = hdf5_file
        self.initialized = hdf5_file is not None

        if hdf5_file is not None:
            self._hdf = PhenotypesHDF5(hdf5_file)  # Existing implementation

    def _check_initialized(self):
        if not self.initialized:
            raise ValueError("Service not properly initialized: missing filename")

    def get_cases_for_phecode(self, phecode, ancestry: str = "EUR") -> pd.DataFrame:
        """Delegate to the underlying HDF5 file's get_cases_for_phecode method."""
        self._check_initialized()
        return self._hdf.get_cases_for_phecode(phecode, ancestry)

    def get_case_counts_for_phecode(self, phecode, ancestry: str = "EUR"):
        """Get the case counts for a given phecode and ancestry."""
        self._check_initialized()
        return self._hdf.get_case_counts_for_phecode(phecode, ancestry)

    def get_phenotypes(
        self,
        eids: np.ndarray | None = None,
        phecodes: np.ndarray | None = None,
        ancestry: str = "EUR",
    ) -> np.ndarray | h5py.Dataset:
        self._check_initialized()
        return self._hdf.get_phenotypes(eids, phecodes, ancestry)

    @property
    def phecodes(self):
        return self._hdf.phecode

    @property
    def phecode_sex(self):
        return self._hdf.disease_sex

    def get_eids(self, ancestry: str = "EUR"):
        return self._hdf.get_eids(ancestry)


class PhenotypesHDF5:
    """Handles access to phenotype data stored in HDF5 format."""

    def __init__(self, hdf5_file: Path):
        self.f = h5py.File(hdf5_file, "r")

        # TODO: Fix the naming of the fields
        phenotype_data = self.f["phenotype_data"]
        eid = self.f["eids"]
        phecode = self.f["phecodes"]
        ancestry = self.f["populations"]
        disease_sex = self.f["phecode_sex"]
        individual_sex = self.f["affected_sex"]

        # Jumping through hoops to fix the type checker...
        assert isinstance(phenotype_data, h5py.Dataset)
        assert isinstance(eid, h5py.Dataset)
        assert isinstance(phecode, h5py.Dataset)
        assert isinstance(ancestry, h5py.Dataset)
        assert isinstance(disease_sex, h5py.Dataset)
        assert isinstance(individual_sex, h5py.Dataset)

        # Sanity checks
        # TODO: Move these to integration tests?
        try:
            print(f"RUNNING SANITY CHECK ON {hdf5_file}! ONLY SEE THIS ONCE!")
            assert phenotype_data.shape[0] == eid.shape[0]
            assert individual_sex.shape[0] == eid.shape[0]
            assert ancestry.shape[0] == eid.shape[0]

            assert phenotype_data.shape[1] == phecode.shape[0]
            assert disease_sex.shape[0] == phecode.shape[0]

            # assert np.all(np.diff(eid) > 0), "EIDs are expected to be sorted"

        except Exception as e:
            print(f"Error in sanity checks: {str(e)}")
            raise

        # Load the data matrix and index information as numpy arrays
        self.phenotype_data_h5: h5py.Dataset = phenotype_data
        self.eid: np.ndarray = eid[...].astype(int)
        self.ancestry: np.ndarray = ancestry[...].astype(str)
        self.individual_sex: np.ndarray = individual_sex[...].astype(str)

        self.phecode: np.ndarray = phecode[...].astype(str)
        self.disease_sex: np.ndarray = disease_sex[...].astype(str)

        # Gah...
        npy_file = hdf5_file.with_suffix(".h5.npy")
        if not npy_file.exists():
            print("REALLY? WE'RE DOING THIS HERE? NOW?")
            np.save(npy_file, self.phenotype_data_h5)

        self.phenotype_data_mm = np.load(npy_file, mmap_mode="r")

        # Pick the memmap version
        self.phenotype_data: np.memmap = self.phenotype_data_mm

        # TODO: Is this useful?
        # TODO: What about eid to index mapping?

        print(f"BUILDING INDEX MAP FOR {hdf5_file}! ONLY SEE THIS ONCE!")
        self.phecode_to_index = {
            phecode: index for index, phecode in enumerate(self.phecode)
        }

        # TODO: Do we need a memory-mapped matrix here too? (See GenotypeService)

        # Precompute sorted indices for faster lookups
        self._phecodes_argsort = np.argsort(self.phecode)
        self._phecodes_sorted = self.phecode[self._phecodes_argsort]
        print(f"PHENOTYPE SERVICE INITIALIZED {hdf5_file}")

    def get_ancestry_mask(self, ancestry):
        """Get a mask for the given ancestry."""
        if ancestry == "EUR":
            return (self.ancestry == "EUR") | (self.ancestry == "EUR_S")
        return self.ancestry == ancestry

    def get_phenotypes(
        self,
        eids: np.ndarray | None = None,
        phecodes: np.ndarray | None = None,
        ancestry: str = "EUR",
    ) -> np.ndarray | h5py.Dataset:
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

        ancestry_mask = self.get_ancestry_mask(ancestry)

        if eids is None and phecodes is None:
            return self.phenotype_data[ancestry_mask, :]

        if eids is not None:
            # Vectorized index lookup using precomputed sorted arrays
            eid_pos = np.searchsorted(self.eid, eids)
            eid_idx = np.argsort(eid_pos)

        if phecodes is not None:
            phecode_pos = np.searchsorted(self._phecodes_sorted, phecodes)
            phecode_idx = self._phecodes_argsort[phecode_pos]

        if eids is None:
            return self.phenotype_data[ancestry_mask, :][:, phecode_idx]

        if phecodes is None:
            return self.phenotype_data[ancestry_mask, :][eid_idx, :]

        return self.phenotype_data[ancestry_mask, :][eid_idx, :][:, phecode_idx]

    def get_eids(self, ancestry: str = "EUR"):
        return self.eid[self.get_ancestry_mask(ancestry)]

    # @profile
    def get_cases_for_phecode(self, phecode, ancestry: str = "EUR") -> pd.DataFrame:
        """Get the cases for a given phecode and ancestry"""
        phecode_index = self.phecode_to_index[phecode]

        ancestry_mask = self.get_ancestry_mask(ancestry)

        return pd.DataFrame(
            {
                "eid": self.eid[ancestry_mask],
                "sex": self.individual_sex[ancestry_mask],
                "phenotype": self.phenotype_data[:, phecode_index][ancestry_mask],
            }
        )

    def get_case_counts_for_phecode(self, phecode, ancestry: str = "EUR"):
        """Get the case counts for a given phecode and ancestry."""

        ancestry_mask = self.get_ancestry_mask(ancestry)
        phecode_index = self.phecode_to_index[phecode]
        phenotypes = self.phenotype_data[:, phecode_index][ancestry_mask]

        values, counts = np.unique(phenotypes, return_counts=True)

        counts_dict = dict(zip(values, counts))

        counts_dict["affected"] = counts_dict.get(1, 0)
        counts_dict["excluded"] = counts_dict.get(9, 0)
        counts_dict["control"] = counts_dict.get(0, 0)

        return counts_dict
