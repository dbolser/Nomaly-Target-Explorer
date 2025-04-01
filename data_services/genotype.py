import logging
import re
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from data_services.utils import deprecated

logger = logging.getLogger(__name__)

chr_str = re.compile("chr", re.IGNORECASE)


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

    def get_genotypes(self, eids=None, vids=None, nomaly_ids=False) -> np.ndarray:
        """Delegate to the underlying HDF5 file's get_genotypes method."""
        self._check_initialized()
        return self._hdf.get_genotypes(eids=eids, vids=vids, nomaly_ids=nomaly_ids)

    def get_variant_counts(
        self,
        nomaly_variant_id: str,
        individual_sex: str | None = None,
        ancestry: str | None = None,
    ) -> dict[str, int]:
        """Delegate to the underlying HDF5 file's get_variant_counts method."""
        self._check_initialized()
        return self._hdf.get_variant_counts(nomaly_variant_id, individual_sex, ancestry)

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

        # # TODO: Hopefully we never need to use this...
        # self.allele_flipped_in_genotype_file_relative_to_nomaly_variant_id = (
        #     self._calculate_allele_flips(self.plink_variant_id, self.nomaly_variant_id)
        # )

        print(f"INITIALIZED {self.f.filename}")

    # TODO: Add sex and ancestry here?
    def get_genotypes(
        self,
        eids: Optional[np.ndarray] = None,
        vids: Optional[np.ndarray] = None,
        nomaly_ids: bool = False,
    ) -> np.ndarray:
        """
        Get genotypes for a list of eids and variant IDs.

        Args:
            eids: List of eids to get genotypes for
            vids: List of variant IDs to get genotypes for.
            genotype_ids:
                If True, the variant IDs will be 'Genotype Format' variant IDs.
                If False, the variant IDs will be 'Nomaly Format' variant IDs.

        Returns:
            np.ndarray: Array of genotypes with shape (n_variants, n_eids) in the order of the
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

    @deprecated("This method will be removed in a future version.")
    def _flip_alleles(self, variant: str) -> str:
        """
        Flip the ref and alt alleles in a variant ID.

        Args:
            variant (str): Variant ID in format "CHR:POS:REF:ALT" or "CHR_POS_REF/ALT"

        Returns:
            str: Flipped variant ID in same format as input
        """
        chrom, pos, ref, alt = variant.split(":")
        return f"{chrom}:{pos}:{alt}:{ref}"

    @deprecated("This method will be removed in a future version.")
    def _standardize_variant_format(self, variant: str) -> str:
        """
        Convert any variant format to CHR:POS:REF:ALT format.

        Handles formats:
        - CHR_POS_REF/ALT (e.g., "8_6870776_C/T")
        - CHR_POS_REF_ALT (e.g., "8_6870776_C_T")
        - CHR:POS:REF:ALT (e.g., "8:6870776:C:T")

        Args:
            variant (str): Variant in any supported format

        Returns:
            str|None: Standardized variant string or None if invalid format
        """
        try:
            chrom, pos, ref, alt = re.split(r"[_/:]", variant)
            # Validate components
            if not chrom or not pos or not ref or not alt:
                raise ValueError(f"Missing component in variant: {variant}")

            if not pos.isdigit():
                raise ValueError(f"Position must be numeric: {pos}")

            # Clean chromosome format (e.g., 'chr8' -> '8')
            chrom = chr_str.sub("", chrom)

            # Standardize to colon format for genotype lookup (why?)
            return f"{chrom}:{pos}:{ref}:{alt}"

        except Exception as e:
            raise ValueError(f"Error standardizing variant format {variant}: {str(e)}")

    def get_variant_counts(
        self,
        nomaly_variant_id: str,
        individual_sex: str | None = None,
        ancestry: str | None = None,
    ) -> dict[str, int]:
        """Get counts of genotypes for a specific NOMALY variant

        Optionally filtered by sex and/or ancestry.

        Args:
            nomaly_variant_id: The NOMALY format variant ID to get counts for
            sex: Optional sex to filter by ('M' or 'F')
            ancestry: Optional ancestry to filter by (must match values in ancestry dataset)

        Returns:
            Dictionary with counts for each genotype category: {
                'homozygous_ref': int,  # Count of 0s 'heterozygous': int,    #
                Count of 1s 'homozygous_alt': int,  # Count of 2s 'missing': int
                # Count of -1s
            }
        """
        # Get genotypes for this variant using nomaly_variant_id
        genotypes = self.get_genotypes(
            vids=np.array([nomaly_variant_id]), nomaly_ids=True
        )

        # Get ancestry mask
        if ancestry is not None:
            sample_mask = self.ancestry == ancestry
        else:
            sample_mask = np.ones(len(self.individual), dtype=bool)
        filtered_genotypes = genotypes[:, sample_mask]

        # TODO: Think about how to keep the variant dimension here... For now,
        # lets just give up. The plan is to use all in one phewas function
        # anyway... isn't it? Thinking about it, the phewas function does one
        # variant across all phenotypes. Here we want to do all variants across
        # one phenotype... SHould be simple eh?

        filtered_genotypes = filtered_genotypes[0]

        # Count occurrences of each genotype value
        return {
            "total": len(filtered_genotypes),
            "homozygous_ref": int(np.sum(filtered_genotypes == 0)),
            "heterozygous": int(np.sum(filtered_genotypes == 1)),
            "homozygous_alt": int(np.sum(filtered_genotypes == 2)),
            "missing": int(np.sum(filtered_genotypes == -9)),  ### FUCK!!!!
            # "missing": int(np.sum(filtered_genotypes == -1)),
        }

    @deprecated("This method will be removed in a future version.")
    def query_variants(self, variant: str) -> np.ndarray:
        """
        Query genotypes for a variant, trying flipped alleles if original not found.

        Args:
            variant (str): Variant ID in any supported format

        Returns:
            np.ndarray: Array of genotypes (shape: (n_variants, n_eids))
        """
        try:
            std_variant = self._standardize_variant_format(variant)

            # Try original variant
            genotypes = self._query_variants_internal(std_variant)
            if genotypes is not None:
                logger.debug(f"Found variant {variant} in original orientation")
                return genotypes

            # Try flipped alleles
            flipped = self._flip_alleles(std_variant)
            genotypes = self._query_variants_internal(flipped)
            if genotypes is not None:
                logger.debug(f"Found variant {variant} in flipped orientation")
                # Because we flipped the alleles, we need to flip the genotypes...
                genotypes[genotypes != -1] = 2 - genotypes[genotypes != -1]
                return genotypes

            logger.warning(f"Variant not found in either orientation: {variant}")
            return np.array([])

        except Exception as e:
            logger.error(f"Error in query_variants: {str(e)}")
            return np.array([])

    @deprecated("This method will be removed in a future version.")
    def _query_variants_internal(self, variant: str) -> np.ndarray | None:
        """Internal method to query a single variant."""
        try:
            variant_mask = self._single_variant_mask(variant)
            if variant_mask is None:
                return None

            # Find matching variant
            selected_variant_indices = np.where(variant_mask)[0]
            if len(selected_variant_indices) == 0:
                return None

            # Query the submatrix
            try:
                submatrix = self.genotype_matrix[selected_variant_indices, :]
                if submatrix.size == 0:
                    return None
                return submatrix
            except Exception as e:
                logger.error(f"Error accessing genotype matrix: {str(e)}")
                return None

        except Exception as e:
            logger.error(f"Error in _query_variants_internal: {str(e)}")
            return None

    @deprecated("This method will be removed in a future version.")
    def _single_variant_mask(self, varsel: str) -> np.ndarray | None:
        """Create mask for single variant."""
        try:
            mask = self.plink_variant_id == varsel
            if mask.sum() == 0:
                return None
            return mask
        except Exception as e:
            logger.error(f"Error in _single_variant_mask: {str(e)}")
            raise

    @deprecated("This method will be removed in a future version.")
    def query_variantID_genotypes(
        self, variant: str
    ) -> tuple[np.ndarray, np.ndarray] | None:
        """
        Get genotype data for a variant.

        Args:
            variant (str): Variant ID in format "CHR:POS:REF:ALT"

        Returns:
            tuple: (sorted_eids, sorted_genotypes) or None if error
        """
        try:
            # Get genotype data
            genotype_eids = self.individual
            genotypes = self.query_variants(variant)
            if genotypes is None or len(genotypes) == 0:
                logger.warning(f"No genotype data found for variant {variant}")
                return None

            # Get first row of genotypes (for single variant)
            genotypes = genotypes[0]
            if len(genotypes) != len(genotype_eids):
                logger.warning(
                    f"Mismatch between genotypes ({len(genotypes)}) and IDs ({len(genotype_eids)})"
                )
                return None

            # Sort the data
            sorted_indices = np.argsort(genotype_eids)
            sorted_genotype_eids = genotype_eids[sorted_indices]
            sorted_genotypes = genotypes[sorted_indices]

            return sorted_genotype_eids, sorted_genotypes

        except Exception as e:
            logger.error(
                f"Error in get_genotype_data for variant (query_variantID_genotypes) {variant}: {str(e)}"
            )
            return None
