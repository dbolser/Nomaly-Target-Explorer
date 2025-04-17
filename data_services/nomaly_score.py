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
        self.data_matrix_h5: h5py.Dataset = scores
        self.eids: np.ndarray = eids[...].astype(int)
        self.terms: np.ndarray = terms[...].astype(str)

        npy_file = hdf5_file.with_suffix(".h5.npy")
        if not npy_file.exists():
            print("REALLY? WE'RE DOING THIS HERE NOW?")
            np.save(npy_file, self.data_matrix_h5)

        self.data_matrix_mm = np.load(npy_file, mmap_mode="r")

        # Pick the memmap version
        self.data_matrix: np.memmap = self.data_matrix_mm

        print(f"SORTING STUFF IN {self.f.filename}")

        # NOTE: We get sorted version of self.eids and the 'sorting indices so
        # we don't have to assume anything about the ordering of the eids
        self.sorting_idx = np.argsort(self.eids)
        self.sorted_eids = self.eids[self.sorting_idx]

        print(f"INITIALIZED {self.f.filename}")

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
            # Convert bool to positions for the requested terms. Using `np.where` ensures we
            # always have an *array* of indices which can safely be passed to ``np.ix_``.
            term_rows = np.where(np.isin(self.terms, terms))[0]
        else:
            # Select all columns
            term_rows = slice(None)

        # Sanity check – all requested eids must be present in the file.
        intersected = np.intersect1d(eids, self.eids)
        assert len(eids) == len(intersected)

        # Find the row indices of the requested eids, taking advantage of the
        # pre‑computed sorted array for speed.
        idx = np.searchsorted(self.sorted_eids, eids)
        eid_rows = self.sorting_idx[idx]

        # ---------------------------------------------------------------------
        # Fancy indexing – we want the result to *always* be 2‑dimensional so
        # that downstream code doesn't have to worry about corner cases such as
        # a single eid or a single term being requested (these would normally
        # lead NumPy to squeeze the corresponding dimension).
        # ---------------------------------------------------------------------
        if isinstance(term_rows, slice):
            # Selecting *all* columns or with a slice already preserves the 2D
            # shape when we pass an integer/array for the rows first.
            result: np.ndarray = self.data_matrix[eid_rows, term_rows]
        else:
            # When both rows and columns are arrays we need ``np.ix_`` to keep
            # the outer (Cartesian) product which guarantees a 2‑D result even
            # if either array has length 1 or 0.
            result = self.data_matrix[np.ix_(eid_rows, term_rows)]

        # If NumPy still managed to squeeze the array into 1‑D (this can happen
        # when one of the dimensions is length 0 and the other is length >0),
        # reshape it explicitly so that callers always get (n_eids, n_terms).
        if result.ndim == 1:
            result = result.reshape(len(eid_rows), -1)

        return result
