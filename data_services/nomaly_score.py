from pathlib import Path

import h5py
import numpy as np


class NomalyScoreService:
    """Service for getting nomaly scores."""

    def __init__(self, hdf5_file: Path | None = None):
        """Initialize the service with the path to the HDF5 file."""
        self.hdf5_file = hdf5_file
        self.initialized = hdf5_file is not None

        if hdf5_file is not None:
            self._hdf = NomalyScoreHDF5(hdf5_file)  # Existing implementation

    def _check_initialized(self):
        if not self.initialized:
            raise ValueError("Service not properly initialized: missing filename")

    def get_scores_by_eids_unsorted(
        self, eids: np.ndarray, terms: np.ndarray | None = None
    ) -> np.ndarray:
        """Delegate to the underlying HDF5 file's get_scores_by_eids_unsorted method."""
        self._check_initialized()
        return self._hdf.get_scores_by_eids_unsorted(eids, terms)


class NomalyScoreHDF5:
    """Handles access to nomaly score data stored in HDF5 format."""

    def __init__(self, hdf5_file: Path):
        self.f = h5py.File(hdf5_file, "r")

        scores = self.f["scores"]
        eids = self.f["eid"]  # TODO: Fix the naming of the fields?
        terms = self.f["term"]  # TODO: Fix the naming of the fields?

        # Jumping through hoops to fix the type checker...
        assert isinstance(scores, h5py.Dataset)
        assert isinstance(eids, h5py.Dataset)
        assert isinstance(terms, h5py.Dataset)

        # Sanity checks
        try:
            print(f"RUNNING SANITY CHECK ON {self.f.filename}! ONLY SEE THIS ONCE!")
            assert scores.shape[0] == eids.shape[0]
            assert scores.shape[1] == terms.shape[0]
        except Exception as e:
            print(f"Error in sanity checks: {str(e)}")
            raise

        # Load the data matrix and index information
        # self.data_matrix: h5py.Dataset = score
        self.eids: np.ndarray = eids[...].astype(int)
        self.terms: np.ndarray = terms[...].astype(str)

        npy_file = hdf5_file.with_suffix(".h5.npy")
        # if not npy_file.exists():
        #     np.save(npy_file, self.data_matrix)

        # Load the data matrix 'memory-mapped'
        # self.data_matrix_mm = np.load(npy_file, mmap_mode="r")
        self.data_matrix = np.load(npy_file)

        # Swap them...
        # self.data_matrix_h5 = self.data_matrix
        # self.data_matrix = self.data_matrix_mm

        # HACK: Get sorted version of self.eids and the 'sorting indices
        self.sorting_idx = np.argsort(self.eids)
        self.sorted_eids = self.eids[self.sorting_idx]

    def get_scores_by_eid(self, eid):
        eid_mask = self.eids == eid
        eid_mask_idx = np.where(eid_mask)[0][0]
        return self.data_matrix[eid_mask_idx, :]

    def get_scores_by_term(self, term):
        term_mask = self.terms == term
        term_mask_idx = np.where(term_mask)[0][0]
        return self.data_matrix[:, term_mask_idx]

    def get_scores_by_eids(self, eids):
        eids_mask = np.isin(self.eids, eids)
        return self.data_matrix[eids_mask, :]

    def get_scores_by_terms(self, terms):
        terms_mask = np.isin(self.terms, terms)
        return self.data_matrix[:, terms_mask]

    def get_scores_by_eids_unsorted(
        self, eids: np.ndarray, terms: np.ndarray | None = None
    ) -> np.ndarray:
        """
        Get scores for a list of eids, unsorted.

        Args:
            eids: Array of eids to get scores for.
            terms: Array of terms to get scores for.
        """
        if terms is not None:
            # Convert bool to positions
            term_rows = np.where(np.isin(self.terms, terms))[0]
        else:
            term_rows = slice(None)  # Select all columns

        # Sanity check for the method to work...
        assert len(eids) == len(np.intersect1d(eids, self.eids))

        # Find where each requested eid would fit in the sorted array of 'score_eids'
        idx = np.searchsorted(self.sorted_eids, eids)

        # Use the original sorting indices to map back to the original unsorted
        # positions
        eid_rows = self.sorting_idx[idx]

        return self.data_matrix[eid_rows, term_rows]
