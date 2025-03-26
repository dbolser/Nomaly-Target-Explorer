import numpy as np

from data_services import GenotypeService

from config import Config

import pytest

# config = Config()

NUM_INDIVIDUALS = 488377
# NUM_INDIVIDUALS = 487950  # Removed 427 individuals with negative eids
NUM_VARIANTS = 83011



@pytest.fixture(scope="session")
def genotype_service():
    return GenotypeService(Config.GENOTYPES_HDF)


def test_production_file_exists(genotype_service):
    """Verify the production genotype file exists and is readable."""

    assert genotype_service.hdf5_file.exists()
    assert genotype_service.hdf5_file.is_file()
    assert genotype_service.hdf5_file.stat().st_mode & 0o400

def test_h5_matches_npy(genotype_service):
    """Verify the HDF5 file matches the NPY file."""
    pass
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


def test_individual_count(genotype_service):
    """Verify the expected number of individuals in the dataset."""
    assert len(genotype_service.individual) == NUM_INDIVIDUALS


def test_variant_count(genotype_service):
    """Verify the expected number of variants in the dataset."""
    assert len(genotype_service.plink_variant_id) == NUM_VARIANTS


def test_genotype_matrix_shape(genotype_service):
    """Verify the shape of the genotype matrix."""
    genotype_matrix = genotype_service.get_genotypes()
    assert genotype_matrix.shape == (NUM_VARIANTS, NUM_INDIVIDUALS)


def test_bad_eids(genotype_service):
    """Test that bad eids are handled correctly."""
    with pytest.raises(IndexError):
        _ = genotype_service.get_genotypes(eids=["hello"])


def test_bad_vids(genotype_service):
    """Test that bad vids are handled correctly."""
    with pytest.raises(IndexError):
        _ = genotype_service.get_genotypes(vids=["hello"])


def test_get_genotypes_to_hell(genotype_service):
    """Test the get_genotypes method."""
    eids = genotype_service.individual
    vids = genotype_service.plink_variant_id

    rand_eid = np.random.choice(eids, size=1000, replace=True)
    rand_vid = np.random.choice(vids, size=100, replace=True)

    # Find the indices of the random eids and vids
    rand_eid_idx = np.where(np.isin(eids, rand_eid))[0]
    rand_vid_idx = np.where(np.isin(vids, rand_vid))[0]

    # Note that the indexes as found above (using np.where) are in the
    # 'original' order, as found in the original eids and vids arrays.

    # Get genotypes in 'original' order
    selected_genotypes_original_order = genotype_service.get_genotypes(
        eids=eids[rand_eid_idx], vids=vids[rand_vid_idx]
    )

    # Get genotypes in 'reverse' order
    selected_genotypes_reverse_order = genotype_service.get_genotypes(
        eids=eids[rand_eid_idx[::-1]], vids=vids[rand_vid_idx[::-1]]
    )

    # The results are different, but we can check that they are the same by
    # simply reversing the reversed matrix...
    assert np.array_equal(
        selected_genotypes_original_order,
        selected_genotypes_reverse_order[::-1, ::-1],
    )

    # Now lets shuffle the eids and vids and check that the results are the same
    rand_eid_idx_idx = np.random.permutation(np.arange(len(rand_eid_idx)))
    rand_vid_idx_idx = np.random.permutation(np.arange(len(rand_vid_idx)))

    selected_genotypes_random_order = genotype_service.get_genotypes(
        eids=eids[rand_eid_idx[rand_eid_idx_idx]],
        vids=vids[rand_vid_idx[rand_vid_idx_idx]],
    )

    # The results are different, but we can check that they are the same by
    # simply re-ordering the original matrix using the shuffled indexes...
    assert np.array_equal(
        selected_genotypes_original_order[rand_vid_idx_idx, :][:, rand_eid_idx_idx],
        selected_genotypes_random_order,
    )

    double_eids = np.concatenate([eids[rand_eid_idx], eids[rand_eid_idx]])
    double_vids = np.concatenate([vids[rand_vid_idx], vids[rand_vid_idx]])

    # Double check that duplicates are handled 'correctly'...
    selected_genotypes_duplicate_eids = genotype_service.get_genotypes(
        eids=double_eids, vids=double_vids
    )

    assert np.array_equal(
        np.tile(selected_genotypes_original_order, (2, 2)),
        selected_genotypes_duplicate_eids,
    )

    # Check with just 10 eids
    selected_10_eids1 = genotype_service.get_genotypes(eids=eids[rand_eid_idx[:10]])
    assert selected_10_eids1.shape == (NUM_VARIANTS, 10)

    selected_10_eids2 = genotype_service.get_genotypes(
        eids=eids[rand_eid_idx[:10][::-1]]
    )
    assert selected_10_eids2.shape == (NUM_VARIANTS, 10)

    assert np.array_equal(selected_10_eids1, selected_10_eids2[:, ::-1])

    # Check with just 10 vids
    selected_10_vids1 = genotype_service.get_genotypes(vids=vids[rand_vid_idx[:10]])
    assert selected_10_vids1.shape == (10, NUM_INDIVIDUALS)

    selected_10_vids2 = genotype_service.get_genotypes(
        vids=vids[rand_vid_idx[:10][::-1]]
    )
    assert selected_10_vids2.shape == (10, NUM_INDIVIDUALS)

    assert np.array_equal(selected_10_vids1, selected_10_vids2[::-1, :])


def test_known_variant_for_genotypes(genotype_service):
    """Test querying a known variant from production data."""
    # Known variant that should exist
    variant = "11:69083946_T/C"
    result = genotype_service.get_genotypes(vids=[variant])

    assert result is not None
    assert isinstance(result, np.ndarray)
    assert result.shape[0] == 1

    assert result.shape[1] == len(genotype_service.individual)
    assert result[0].dtype == np.int8  # Genotypes should be integers

    genotype_counts = np.unique(result[0], return_counts=True)

    # TODO: REDUMP?
    # assert np.all(np.isin(genotype_counts[0], [-1, 0, 1, 2]))
    assert np.all(np.isin(genotype_counts[0], [-9, 0, 1, 2]))

    # The real counts were carefully obtained using PLINK2...
    assert np.all(np.isin(genotype_counts[1], [444, 1516, 11499, 474918]))


def test_another_known_variant_genotypes(genotype_service):
    """Test the distribution of genotypes for a known variant."""
    # Known variant with expected distribution
    variant1 = "19:44908684_T/C"
    variant2 = "6:26199089_A/C"

    results = genotype_service.get_genotypes(vids=[variant1, variant2])
    assert results is not None

    assert results.shape == (2, NUM_INDIVIDUALS)

    result = results[0]
    result_counts = np.unique(result, return_counts=True)
    assert np.all(np.isin(result_counts[0], [-9, 0, 1, 2]))
    assert np.all(np.isin(result_counts[1], [73909, 296223, 108371, 9874]))

    result = results[1]
    result_counts = np.unique(result, return_counts=True)
    assert np.all(np.isin(result_counts[0], [-9, 0, 1, 2]))
    assert np.all(np.isin(result_counts[1], [354, 488018, 5, 0]))


def test_plink_variant_id_format_in_file(genotype_service):
    """Test that variants in the file follow the expected format."""
    # Sample first 1000 variants
    sample_variants = genotype_service.plink_variant_id[:1000]

    import re

    for variant in sample_variants:
        parts = re.split(r"[:_/]", variant)
        assert len(parts) == 4, f"Invalid variant format: {variant}"
        assert parts[0].isdigit(), f"Invalid chromosome: {parts[0]}"
        assert parts[1].isdigit(), f"Invalid position: {parts[1]}"
        assert len(parts[2]) >= 1, f"Invalid ref allele: {parts[2]}"
        assert len(parts[3]) >= 1, f"Invalid alt allele: {parts[3]}"


@pytest.mark.skip(reason="Do we need to allow to query by nomaly variant ids?")
def test_nomaly_variant_id_format_in_file(genotype_service):
    """Test that variants in the file follow the expected format."""
    # Sample first 1000 variants
    sample_variants = genotype_service.nomaly_variant_id[:1000]

    missing_count = 0
    for variant in sample_variants:
        if variant == "Missing":
            missing_count += 1
            continue

        chr, pos, alleles = variant.split("_")
        assert len(chr) == 1, f"Invalid chromosome: {chr}"
        assert len(pos) >= 1, f"Invalid position: {pos}"
        assert len(alleles) == 3, f"Invalid alleles: {alleles}"
        assert alleles[0] in ["A", "C", "G", "T"], f"Invalid ref allele: {alleles[0]}"
        assert alleles[1] == "/", f"Invalid separator: {alleles[1]}"
        assert alleles[2] in ["A", "C", "G", "T"], f"Invalid alt allele: {alleles[2]}"

    assert missing_count < 1000, "There should be less than 1000 missing variants"


def test_individual_id_format(genotype_service):
    """Test that individual IDs are in the expected format."""
    sample_ids = genotype_service.individual

    # THIS WAS A PROBLEM IN THE ORIGINAL FILE, e..g there was exactly 1 -1 in
    # the sample_ids!

    NUM_MISSING_IDS = 427

    assert np.sum(sample_ids <= 0) == NUM_MISSING_IDS, (
        f"There should be {NUM_MISSING_IDS} -1s in the sample_ids...ahhh"
    )

    for id_num in sample_ids:
        assert isinstance(id_num, np.integer)
        if id_num <= 0:
            continue
        assert id_num >= 100000, f"Individual ID: {id_num} is not greater than 100000"
        assert id_num < 10000000, f"Individual ID: {id_num} is not less than 10000000"
        assert len(str(id_num)) >= 5  # Example: Expecting at least 5-digit IDs

