import logging
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from data_services.utils import deprecated

logger = logging.getLogger(__name__)


class GenotypeService:
    """Service for getting genotype data."""

    def __init__(self, hdf5_file: Path | None = None):
        self.hdf5_file = hdf5_file
        self.initialized = hdf5_file is not None

        if hdf5_file is not None:
            self._hdf = GenotypesHDF5(hdf5_file)  # Existing implementation

    def _check_initialized(self):
        if not self.initialized:
            raise ValueError("Service not properly initialized: missing filename")

    def get_genotypes(
        self,
        eids: np.ndarray | None = None,
        vids: np.ndarray | None = None,
        nomaly_ids: bool = False,
    ) -> np.ndarray:
        """Delegate to the underlying HDF5 file's method."""
        self._check_initialized()
        return self._hdf.get_genotypes(eids=eids, vids=vids, nomaly_ids=nomaly_ids)

    def get_genotype_counts_and_freqs(
        self,
        eids: np.ndarray | None = None,
        vids: np.ndarray | None = None,
        nomaly_ids: bool = False,
        ancestry: str = "EUR",
    ) -> pd.DataFrame:
        """Delegate to the underlying HDF5 file's method."""
        self._check_initialized()
        return self._hdf.get_genotype_counts_and_freqs(
            eids=eids, vids=vids, nomaly_ids=nomaly_ids, ancestry=ancestry
        )

    @deprecated("Use get_genotype_counts_and_freqs instead.")
    def get_variant_counts(
        self, nomaly_variant_id: str, ancestry: str = "EUR"
    ) -> dict[str, int]:
        """Delegate to the underlying HDF5 file's get_variant_counts method."""
        self._check_initialized()
        return self._hdf.get_variant_counts(nomaly_variant_id, ancestry)

    @property
    def individual(self):
        """Delegate to the underlying HDF5 file's individual property."""
        self._check_initialized()
        return self._hdf.individual

    @property
    def plink_variant_ids(self) -> np.ndarray:
        """Delegate to the underlying HDF5 file's genotype_variant_id property."""
        self._check_initialized()
        return self._hdf.plink_variant_id

    @property
    def nomaly_variant_ids(self) -> np.ndarray:
        """Delegate to the underlying HDF5 file's genotype_variant_id property."""
        self._check_initialized()
        return self._hdf.nomaly_variant_id


class GenotypesHDF5:
    """Handles access to genotype data stored in HDF5 format.

    This class expects both an HDF5 file and a corresponding .npy file
    containing the same genotype matrix data. The .npy file is used for
    memory-mapped access to improve performance.

    Args:
        hdf5_file (str): Path to the HDF5 file. A corresponding .npy file
            must exist at {hdf5_file}.npy
    """

    def __init__(self, hdf5_file: Path):
        self.f = h5py.File(hdf5_file, "r")

        genotype_matrix = self.f["genotypes"]

        # COLUMNS
        individual = self.f["eid"]
        individual_sex = self.f["sex"]
        ancestry = self.f["ancestry"]

        # ROWS
        plink_variant_id = self.f["plink_variant_id"]
        nomaly_variant_id = self.f["nomaly_variant_id"]
        reference_allele = self.f["REF"]
        alternate_allele = self.f["ALT"]
        rsid = self.f["rsID"]

        # NOTE: We used to explicitly store genotype counts in the HDF5 file
        # genotype_counts = self.f["genotype_counts"]

        # Jumping through hoops to fix the type checker...
        assert isinstance(genotype_matrix, h5py.Dataset)

        assert isinstance(individual, h5py.Dataset)
        assert isinstance(individual_sex, h5py.Dataset)
        assert isinstance(ancestry, h5py.Dataset)

        assert isinstance(plink_variant_id, h5py.Dataset)
        assert isinstance(nomaly_variant_id, h5py.Dataset)
        assert isinstance(reference_allele, h5py.Dataset)
        assert isinstance(alternate_allele, h5py.Dataset)
        assert isinstance(rsid, h5py.Dataset)

        # Sanity checks
        try:
            print(f"RUNNING SANITY CHECK ON {self.f.filename}")

            assert genotype_matrix.shape[0] == plink_variant_id.shape[0]
            assert genotype_matrix.shape[1] == individual.shape[0]

            assert individual.shape[0] == individual_sex.shape[0]
            assert individual.shape[0] == ancestry.shape[0]

            assert plink_variant_id.shape[0] == nomaly_variant_id.shape[0]
            assert plink_variant_id.shape[0] == reference_allele.shape[0]
            assert plink_variant_id.shape[0] == alternate_allele.shape[0]
            assert plink_variant_id.shape[0] == rsid.shape[0]
        except Exception as e:
            print(f"Error in sanity checks: {str(e)}")
            raise

        # Load the data matrix and index information as numpy arrays

        # self.genotype_matrix: h5py.Dataset = genotype_matrix  # NOTE: SEE BELOW!
        self.individual: np.ndarray = individual[...].astype(int)
        self.individual_sex: np.ndarray = individual_sex[...].astype(str)
        self.ancestry: np.ndarray = ancestry[...].astype(str)

        self.plink_variant_id: np.ndarray = plink_variant_id[...].astype(str)
        self.nomaly_variant_id: np.ndarray = nomaly_variant_id[...].astype(str)
        self.rsid: np.ndarray = rsid[...].astype(str)

        # NOTE:Load genotype_matrix as memmap!
        # Load NPY file # TODO: What is mode?
        npy_file = hdf5_file.with_suffix(".npy")
        self.genotype_matrix = np.memmap(
            npy_file, dtype=np.int8, mode="r", shape=genotype_matrix.shape
        )

        print(f"SORTING STUFF IN {self.f.filename}")

        # Precompute sorted indices for faster row / column lookups

        self._individual_argsort = np.argsort(self.individual)
        self._individual_sorted = self.individual[self._individual_argsort]

        self._plink_variant_argsort = np.argsort(self.plink_variant_id)
        self._plink_variant_sorted = self.plink_variant_id[self._plink_variant_argsort]

        # TODO: Could we use a mapping here instead?
        self._nomaly_variant_argsort = np.argsort(self.nomaly_variant_id)
        self._nomaly_variant_sorted = self.nomaly_variant_id[
            self._nomaly_variant_argsort
        ]

        print(f"INITIALIZED {self.f.filename}")

    def get_genotypes(
        self,
        eids: np.ndarray | None = None,
        vids: np.ndarray | None = None,
        nomaly_ids: bool = False,
    ) -> np.ndarray:
        """
        Get genotypes for a list of eids and variant IDs.

        NOTE:
            We don't filter by ancestry here. Use a list of EIDs, or fetch the
            ancestry separately if you want to do that.

        Args:
            eids:
                An optional list of eids to get genotypes for
            vids:
                An optional list of variant IDs to get genotypes for
            nomaly_ids:
                If True, the variant IDs will be 'Nomaly Format' variant IDs. If
                False, the variant IDs will be 'Genotype Format' variant IDs.

        Returns:
            np.ndarray: Array of genotypes with shape (n_variants, n_eids) in
            the order of the
                provided eids and variant IDs.

            Genotype values are encoded as:
                - 0 = homozygous ref
                - 1 = heterozygous
                - 2 = homozygous alt
                - -1 = missing
        """

        if eids is None and vids is None:
            return self.genotype_matrix

        if eids is not None:
            # Vectorized index lookup using precomputed sorted arrays
            eid_pos = np.searchsorted(self._individual_sorted, eids)
            if np.any(eid_pos == len(self._individual_sorted)):
                raise IndexError(
                    f"Individual {eids} not found in genotype file. Please check your individual IDs."
                )
            eid_idx = self._individual_argsort[eid_pos]

        if vids is not None:
            if nomaly_ids:
                vid_pos = np.searchsorted(self._nomaly_variant_sorted, vids)
                if np.any(vid_pos == len(self._nomaly_variant_sorted)):
                    raise IndexError(
                        f"Variant {vids} not found in genotype file. Please check your (nomaly) variant IDs."
                    )
                vid_idx = self._nomaly_variant_argsort[vid_pos]
            else:
                vid_pos = np.searchsorted(self._plink_variant_sorted, vids)
                if np.any(vid_pos == len(self._plink_variant_sorted)):
                    raise IndexError(
                        f"Variant {vids} not found in genotype file. Please check your (plink) variant IDs."
                    )
                vid_idx = self._plink_variant_argsort[vid_pos]

        # It was looking so clean until here...

        if eids is not None:
            eid_bad = self.individual[eid_idx] != eids
            if np.any(eid_bad):
                raise IndexError(
                    f"Individual {eids} not found in genotype file. Please check your individual IDs."
                )

            # NOTE: One option would be to return a matrix with -9s for the
            # missing values (corresponding to non-existant eids)...
            # matrix = self.genotype_matrix[vid_idx, :].copy()
            # matrix[:, eid_bad] = -9
            # return matrix

        if vids is not None:
            if nomaly_ids:
                vid_bad = self.nomaly_variant_id[vid_idx] != vids
            else:
                vid_bad = self.plink_variant_id[vid_idx] != vids

            if np.any(vid_bad):
                raise IndexError(
                    f"Variant {vids} not found in genotype file. Please check your variant IDs."
                )
            # NOTE: One option would be to return a matrix with -9s for the
            # missing values (corresponding to non-existant vids)...
            # matrix = self.genotype_matrix[vid_idx, :].copy()
            # matrix[vid_bad, :] = -9
            # return matrix

        if eids is None:
            return self.genotype_matrix[vid_idx, :]

        if vids is None:
            return self.genotype_matrix[:, eid_idx]

        return self.genotype_matrix[vid_idx, :][:, eid_idx]

    def get_genotype_counts_and_freqs(
        self,
        eids: np.ndarray | None = None,
        vids: np.ndarray | None = None,
        nomaly_ids: bool = False,
        ancestry: str = "EUR",
    ) -> pd.DataFrame:
        """Get counts of genotypes for a specific NOMALY variant.

        NOTE: The ancestry param is a bit of a hack...
        """

        # This is rather annoying...
        if eids is None:
            ancestry_mask = self.ancestry == ancestry
        else:
            ancestry_mask = np.ones(eids.shape, dtype=bool)

        genotypes = self.get_genotypes(eids=eids, vids=vids, nomaly_ids=nomaly_ids)

        # Get the individual genotype masks
        het_genotype_mask = genotypes == 1
        ref_genotype_mask = genotypes == 0
        alt_genotype_mask = genotypes == 2
        mis_genotype_mask = genotypes == -9

        # Get the counts
        het_count = np.sum(het_genotype_mask & ancestry_mask)
        ref_count = np.sum(ref_genotype_mask & ancestry_mask)
        alt_count = np.sum(alt_genotype_mask & ancestry_mask)
        mis_count = np.sum(mis_genotype_mask & ancestry_mask)

        # Get the frequencies
        total = het_count + ref_count + alt_count + mis_count
        het_freq = het_count / total
        ref_freq = ref_count / total
        alt_freq = alt_count / total
        mis_freq = mis_count / total

        # Return the counts and frequencies as a dataframe
        return pd.DataFrame(
            {
                "variant_id": vids,
                "het_count": het_count,
                "ref_count": ref_count,
                "alt_count": alt_count,
                "mis_count": mis_count,
                "het_freq": het_freq,
                "ref_freq": ref_freq,
                "alt_freq": alt_freq,
                "mis_freq": mis_freq,
            }
        )

    @deprecated(
        "Use get_genotype_counts instead. This will be removed in a future version."
    )
    def get_variant_counts(
        self,
        nomaly_variant_id: str,
        ancestry: str = "EUR",
    ) -> dict[str, int]:
        """Get counts of genotypes for a specific NOMALY variant

        Args:
            nomaly_variant_id: The NOMALY format variant ID to get counts for
            ancestry: Ancestry to filter by (must match values in ancestry dataset)

        Returns:
            Dictionary with counts for each genotype category: {
                'total': int,
                'homozygous_ref': int,  # Count of 0s
                'heterozygous': int,    # Count of 1s
                'homozygous_alt': int,  # Count of 2s
                'missing': int,         # Count of -1s
            }
        """

        counts = self.get_genotype_counts_and_freqs(
            vids=np.array([nomaly_variant_id]),
            nomaly_ids=True,
            ancestry=ancestry,
        )
        counts = counts.iloc[0]

        # Count occurrences of each genotype value
        return {
            "total": counts.sum(),
            "homozygous_ref": counts[0],
            "heterozygous": counts[1],
            "homozygous_alt": counts[2],
            "missing": counts[3],
        }
