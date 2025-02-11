from pathlib import Path
from typing import Optional

import h5py
import numpy as np


class NomalyScoreService:
    def __init__(self, hdf5_file: Path):
        self._hdf = NomalyScoresHDF5(hdf5_file)  # Existing implementation


class NomalyScoresHDF5:
    """Handles access to nomaly score data stored in HDF5 format."""

    def __init__(self, hdf5_file: Path):
        self.f = h5py.File(hdf5_file, "r")

        # Jumping through hoops to fix the type checker...
        # TODO: Fix the naming of the fields
        try:
            nomaly_scores = self.f["scores"]
            eids = self.f["eid"]
            terms = self.f["term"]

            # Sanity checks
            assert isinstance(nomaly_scores, h5py.Dataset)
            assert isinstance(eids, h5py.Dataset)
            assert isinstance(terms, h5py.Dataset)

            # Sanity checks
            assert nomaly_scores.shape[0] == eids.shape[0]
            assert nomaly_scores.shape[1] == terms.shape[0]

        except Exception as e:
            print(f"Error in sanity checks: {str(e)}")
            raise

        # Load the data matrix and index information
        self.data_matrix: h5py.Dataset = nomaly_scores
        self.eids: np.ndarray = eids[...].astype(int)
        self.terms: np.ndarray = terms[...].astype(str)

        # FUCK! (HACK No. 1)
        npy_file = hdf5_file.with_suffix(".h5.npy")
        if not npy_file.exists():
            np.save(npy_file, self.data_matrix)

        # Load the data matrix 'memory-mapped'
        self.data_matrix_mm = np.load(npy_file, mmap_mode="r")

        # Swap them...
        self.data_matrix_h5 = self.data_matrix
        self.data_matrix = self.data_matrix_mm

        # DOUBLE FUCK!! (HACK No. 2)
        # Get sorted version of self.eids and the indices that would sort it
        self.sorted_idx = np.argsort(self.eids)
        self.sorted_eids = self.eids[self.sorted_idx]

    def get_scores_by_eid(self, eid):
        mask_row = self.eids == eid
        mask_row_indices = np.where(mask_row)[0][0]
        return self.data_matrix[mask_row_indices, :]

    def get_scores_by_term(self, term):
        mask_column = self.terms == term
        mask_column_indices = np.where(mask_column)[0][0]
        return self.data_matrix[:, mask_column_indices]

    def get_scores_by_eids(self, eids):
        mask_rows = np.isin(self.eids, eids)
        return self.data_matrix[mask_rows, :]

    def get_scores_by_terms(self, terms):
        mask_columns = np.isin(self.terms, terms)
        return self.data_matrix[:, mask_columns]

    def get_scores_by_eids_unsorted(
        self, eids: np.ndarray, terms: Optional[np.ndarray] = None
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
        eid_rows = self.sorted_idx[idx]

        return self.data_matrix[eid_rows, term_rows]
