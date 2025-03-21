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
    def __init__(self, hdf5_file: Path | str):
        # Convert string to Path if needed
        if isinstance(hdf5_file, str):
            hdf5_file = Path(hdf5_file)
        self._hdf = GenotypesHDF5(hdf5_file)  # Existing implementation

    @property
    def individual(self):
        """Delegate to the underlying HDF5 file's individual property."""
        return self._hdf.individual

    def query_variantID_genotypes(
        self, variant: str
    ) -> tuple[np.ndarray, np.ndarray] | None:
        """Delegate to the underlying HDF5 file's query_variantID_genotypes method."""
        return self._hdf.query_variantID_genotypes(variant)

    def get_genotypes(self, eids=None, vids=None, nomaly_ids=False) -> np.ndarray:
        """Delegate to the underlying HDF5 file's get_genotypes method."""
        return self._hdf.get_genotypes(eids=eids, vids=vids, nomaly_ids=nomaly_ids)


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

        # TODO: Fix the naming of the fields!
        genotype_matrix = self.f["genotype_matrix"]
        individual = self.f["fam"]
        biological_sex = self.f["sex"]
        ancestry = self.f["ancestry"]

        genotype_variant_id = self.f["bim"]
        nomaly_variant_id = self.f["nomaly_variant_id"]
        plink_variant_id = self.f["plink_variant_id"]

        # genotype_counts = self.f["genotype_counts"]

        # Jumping through hoops to fix the type checker...
        assert isinstance(genotype_matrix, h5py.Dataset)
        assert isinstance(individual, h5py.Dataset)
        assert isinstance(biological_sex, h5py.Dataset)
        assert isinstance(ancestry, h5py.Dataset)
        assert isinstance(genotype_variant_id, h5py.Dataset)
        assert isinstance(nomaly_variant_id, h5py.Dataset)
        assert isinstance(plink_variant_id, h5py.Dataset)

        # Sanity checks
        try:
            print(f"RUNNING SANITY CHECK ON {self.f.filename}")
            # assert isinstance(genotype_counts, h5py.Dataset)

            assert genotype_matrix.shape[0] == genotype_variant_id.shape[0]
            assert genotype_matrix.shape[1] == individual.shape[0]

            assert biological_sex.shape[0] == individual.shape[0]
            assert ancestry.shape[0] == individual.shape[0]

            assert nomaly_variant_id.shape[0] == genotype_variant_id.shape[0]

            # assert genotype_counts.shape[0] == genotype_variant_id.shape[0]

        except Exception as e:
            print(f"Error in sanity checks: {str(e)}")
            raise

        # Load the data matrix and index information as numpy arrays
        # self.genotype_matrix: h5py.Dataset = genotype_matrix
        self.individual: np.ndarray = individual[...].astype(int)
        self.biological_sex: np.ndarray = biological_sex[...].astype(str)
        self.ancestry: np.ndarray = ancestry[...].astype(str)

        self.genotype_variant_id: np.ndarray = genotype_variant_id[...].astype(str)
        self.nomaly_variant_id: np.ndarray = nomaly_variant_id[...].astype(str)
        self.plink_variant_id: np.ndarray = plink_variant_id[...].astype(str)

        # self.genotype_counts: np.ndarray = genotype_counts[...].astype(int)

        # TODO: Convert genotype matrix to np.memmap here? (Currently we
        # assume that the corresponding .npy file already exists.)
        genotype_matrix_np_path = f"{self.f.file.filename}.npy"

        # Load genotype_matrix as memmap!
        self.genotype_matrix: np.ndarray = np.load(
            genotype_matrix_np_path,
            mmap_mode="r",  # TODO: What is mode?
        )

        # More sanity checks

        assert np.all(np.diff(self.individual) > 0), "EIDs are expected to be sorted"

        # # TODO: Move these to integration tests?
        # try:
        #     assert self.genotype_matrix_mm.shape == self.genotype_matrix.shape
        #     assert self.genotype_matrix_mm.dtype == self.genotype_matrix.dtype

        #     # A bit random, but hey...
        #     assert np.all(
        #         self.genotype_matrix[0:10, :] == self.genotype_matrix_mm[0:10, :]
        #     )
        #     assert np.all(
        #         self.genotype_matrix[:, 0:10] == self.genotype_matrix_mm[:, 0:10]
        #     )
        # except Exception as e:
        #     print(
        #         f"Error in sanity checks, HDF5 genotype_matrix != memmap: {str(e)}"
        #     )
        #     raise

        # # Switch to mm (for testing)
        # self.genotype_matrix_h5 = self.genotype_matrix
        # self.genotype_matrix = self.genotype_matrix_mm

        # Above we try asserting that eid's are sorted, but the following is
        # fast enough, I think it's probably preferable to not care and do
        # this anyway. It's a small up-front cost for a lot of downstream
        # benefit, and praying that we sorted everything as expectd is a
        # burden to carry.

        # Precompute sorted indices for faster lookups
        self._individual_argsort = np.argsort(self.individual)
        self._individual_sorted = self.individual[self._individual_argsort]

        self._nomaly_variant_argsort = np.argsort(self.nomaly_variant_id)
        self._nomaly_variant_sorted = self.nomaly_variant_id[
            self._nomaly_variant_argsort
        ]

        self._genotype_variant_argsort = np.argsort(self.genotype_variant_id)
        self._genotype_variant_sorted = self.genotype_variant_id[
            self._genotype_variant_argsort
        ]

        # Hopefully we never need to use this...
        self.allele_flipped_in_genotype_file_relative_to_nomaly_variant_id = (
            self._calculate_allele_flips(
                self.nomaly_variant_id, self.genotype_variant_id
            )
        )

    @staticmethod
    def _calculate_allele_flips(
        nomaly_variant_id: np.ndarray,
        genotype_variant_id: np.ndarray,
    ) -> np.ndarray:
        """
        Returns an array indicating whether the alleles are flipped in the
        genotype file relative to the nominal variant ID.

        True = alleles are flipped, e.g.
          - 3:122257322:A:G (genotype file)            or 3:122257322:A/G (plink
            file) vs.
          - 3_122257322_G/A (Nomaly variant ID... ish) or 3_122257322_G_A
            (Nomaly variant ID).

            Isn't it fun to be alive?

        False = alleles are not flipped.

        None = There is no Nomaly variant ID for this particular variant in the
        genotype file.

        """

        result = np.full(len(nomaly_variant_id), None)

        for i in range(len(nomaly_variant_id)):
            if nomaly_variant_id[i] == "Missing":
                continue

            ref_match = nomaly_variant_id[i][-3] == genotype_variant_id[i][-3]
            alt_match = nomaly_variant_id[i][-1] == genotype_variant_id[i][-1]

            if ref_match and alt_match:
                result[i] = False  # alleles are not flipped
            elif not ref_match and not alt_match:
                result[i] = True  # alleles are flipped
            else:
                # Fuck knows!
                raise ValueError(
                    f"Variant {nomaly_variant_id[i]} has different alleles in the genotype file and the nominal variant ID"
                )

        return result

    # Example class property method:
    @property
    def num_variants(self) -> int:
        """
        Returns the number of variants in the genotype matrix.

        Returns:
            int: The number of variants
        """
        return self.genotype_matrix.shape[0]

    # TODO: Add sex and ancestry here!
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
            eid_idx = self._individual_argsort[eid_pos]

        if vids is not None:
            if nomaly_ids:
                vid_pos = np.searchsorted(self._nomaly_variant_sorted, vids)
                vid_idx = self._nomaly_variant_argsort[vid_pos]
            else:
                vid_pos = np.searchsorted(self._genotype_variant_sorted, vids)
                vid_idx = self._genotype_variant_argsort[vid_pos]

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

    @deprecated("This method will be removed in a future version.")
    def get_variant_counts(
        self,
        nomaly_variant_id: str,
        sex: Optional[str] = None,
        ancestry: Optional[str] = None,
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

        # Get variant index
        variant_idx = np.where(self.nomaly_variant_id == nomaly_variant_id)

        if len(variant_idx[0]) == 0:
            raise ValueError(f"Variant {nomaly_variant_id} not found")

        if len(variant_idx[0]) != 1:
            raise ValueError(
                f"Expected exactly one variant index for {nomaly_variant_id}, got {len(variant_idx[0])}"
            )

        # Get genotypes for this variant
        genotypes = self.genotype_matrix[variant_idx[0], :][0]

        # Create sample mask starting with all True
        sample_mask = np.ones(len(self.individual), dtype=bool)

        # Apply sex filter if specified
        if sex is not None:
            sex_mask = self.biological_sex == sex
            sample_mask &= sex_mask

        # Apply ancestry filter if specified
        if ancestry is not None:
            ancestry_data = self.ancestry
            sample_mask &= ancestry_data == ancestry

        # Apply sample mask to genotypes
        filtered_genotypes = genotypes[sample_mask]

        # Count occurrences of each genotype value
        return {
            "total": len(filtered_genotypes),
            "homozygous_ref": int(np.sum(filtered_genotypes == 0)),
            "heterozygous": int(np.sum(filtered_genotypes == 1)),
            "homozygous_alt": int(np.sum(filtered_genotypes == 2)),
            "missing": int(np.sum(filtered_genotypes == -1)),
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

    def _single_variant_mask(self, varsel: str) -> np.ndarray | None:
        """Create mask for single variant."""
        try:
            mask = self.genotype_variant_id == varsel
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
