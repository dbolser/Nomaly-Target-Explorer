from pathlib import Path

import h5py
import numpy as np


class NomalyScoreService:
    def __init__(self, hdf5_file: Path | str):
        self._hdf = NomalyScoresHDF5(hdf5_file)  # Existing implementation


class NomalyScoresHDF5:
    """Handles access to nomaly score data stored in HDF5 format."""

    def __init__(self, hdf5_file):
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

    def get_score_by_eid(self, eid):
        mask_row = self.eids == eid
        mask_row_indices = np.where(mask_row)[0][0]
        return self.data_matrix[mask_row_indices, :]

    def get_score_by_term(self, term):
        mask_column = self.terms == term
        mask_column_indices = np.where(mask_column)[0][0]
        return self.data_matrix[:, mask_column_indices]

    def get_scores_by_eids(self, eids):
        mask_rows = np.isin(self.eids, eids)
        return self.data_matrix[mask_rows, :]

    def get_scores_by_terms(self, terms):
        mask_columns = np.isin(self.terms, terms)
        return self.data_matrix[:, mask_columns]
