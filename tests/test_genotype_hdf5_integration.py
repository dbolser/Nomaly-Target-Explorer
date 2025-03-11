import numpy as np

from data_services.genotype import GenotypesHDF5

from config import Config

import pytest

# config = Config()

NUM_INDIVIDUALS = 488377
NUM_INDIVIDUALS = 487950  # Removed 427 individuals with negative eids
NUM_VARIANTS = 83011



@pytest.fixture(scope="session")
def nomaly_genotype(integration_app):
    return integration_app.extensions["nomaly_services"].genotype._hdf


def test_production_file_exists(nomaly_genotype):
    """Verify the production genotype file exists and is readable."""


    assert nomaly_genotype is not None
    assert hasattr(nomaly_genotype, "f")
    assert nomaly_genotype.hdf5_file.exists()
    assert nomaly_genotype.hdf5_file.is_file()
    assert nomaly_genotype.hdf5_file.stat().st_mode & 0o400


def test_known_variant_query(nomaly_genotype):
    """Test querying a known variant from production data."""
    # Known variant that should exist
    variant = "11:69083946:T:C"
    result = nomaly_genotype.query_variants(variant)

    assert result is not None
    assert isinstance(result, np.ndarray)
    assert result.shape[1] == len(
        nomaly_genotype.individual
    )  # Should match number of individuals
    assert result[0].dtype == np.int8  # Genotypes should be integers
    assert all(g in [-1, 0, 1, 2] for g in result[0])  # Valid genotype values


def test_individual_count(nomaly_genotype):
    """Verify the expected number of individuals in the dataset."""
    assert len(nomaly_genotype.individual) == NUM_INDIVIDUALS


def test_variant_count(nomaly_genotype):
    """Verify the expected number of variants in the dataset."""
    assert len(nomaly_genotype.genotype_variant_id) == NUM_VARIANTS
    assert len(nomaly_genotype.nomaly_variant_id) == NUM_VARIANTS


def test_genotype_matrix_shape(nomaly_genotype):
    """Verify the shape of the genotype matrix."""
    assert nomaly_genotype.genotype_matrix.shape == (NUM_VARIANTS, NUM_INDIVIDUALS)


def test_get_genotypes_to_hell(nomaly_genotype):
    """Test the get_genotypes method."""
    eids = nomaly_genotype.individual
    vids = nomaly_genotype.genotype_variant_id

    rand_eid = np.random.choice(eids, size=1000, replace=True)
    rand_vid = np.random.choice(vids, size=100, replace=True)

    # Find the indices of the random eids and vids
    rand_eid_idx = np.where(np.isin(eids, rand_eid))[0]
    rand_vid_idx = np.where(np.isin(vids, rand_vid))[0]

    # Note that the indexes as found above (using np.where) are in the
    # 'original' order, as found in the original eids and vids arrays.

    # Get genotypes in 'original' order
    selected_genotypes_original_order = nomaly_genotype.get_genotypes(
        eids=eids[rand_eid_idx], vids=vids[rand_vid_idx]
    )

    # Get genotypes in 'reverse' order
    selected_genotypes_reverse_order = nomaly_genotype.get_genotypes(
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

    selected_genotypes_random_order = nomaly_genotype.get_genotypes(
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
    selected_genotypes_duplicate_eids = nomaly_genotype.get_genotypes(
        eids=double_eids, vids=double_vids
    )

    assert np.array_equal(
        np.tile(selected_genotypes_original_order, (2, 2)),
        selected_genotypes_duplicate_eids,
    )

    # Check with just 10 eids
    selected_10_eids1 = nomaly_genotype.get_genotypes(eids=eids[rand_eid_idx[:10]])
    assert selected_10_eids1.shape == (NUM_VARIANTS, 10)

    selected_10_eids2 = nomaly_genotype.get_genotypes(
        eids=eids[rand_eid_idx[:10][::-1]]
    )
    assert selected_10_eids2.shape == (NUM_VARIANTS, 10)

    assert np.array_equal(selected_10_eids1, selected_10_eids2[:, ::-1])

    # Check with just 10 vids
    selected_10_vids1 = nomaly_genotype.get_genotypes(vids=vids[rand_vid_idx[:10]])
    assert selected_10_vids1.shape == (10, NUM_INDIVIDUALS)

    selected_10_vids2 = nomaly_genotype.get_genotypes(
        vids=vids[rand_vid_idx[:10][::-1]]
    )
    assert selected_10_vids2.shape == (10, NUM_INDIVIDUALS)

    assert np.array_equal(selected_10_vids1, selected_10_vids2[::-1, :])

    # Check with garbage eids and vids
    # Currently we're fucked...
    garbage = np.array(["foo", "bar", "baz"])

    with pytest.raises(IndexError):
        _ = nomaly_genotype.get_genotypes(eids=garbage)

    # This is like the least best result we could expect... other than what we currently do...
    # assert selected_garbage_eids.shape == (NUM_VARIANTS, 0)

    with pytest.raises(IndexError):
        _ = nomaly_genotype.get_genotypes(vids=garbage)

    # This is like the least best result we could expect... other than what we currently do...
    # assert selected_garbage_vids.shape == (0, NUM_INDIVIDUALS)


def test_known_variant_genotype_distribution(nomaly_genotype):
    """Test the distribution of genotypes for a known variant."""
    # Known variant with expected distribution
    variant = "19:44908684:C:T"
    result = nomaly_genotype.query_variants(variant)
    assert result is not None

    genotypes = result[0]

    expected_missing_______ = 73846
    expected_ref_homozygous = 9863
    expected_heterozygous__ = 108276
    expected_alt_homozygous = 295965

    assert np.sum(genotypes == -1) == expected_missing_______
    assert np.sum(genotypes == 0) == expected_ref_homozygous
    assert np.sum(genotypes == 1) == expected_heterozygous__
    assert np.sum(genotypes == 2) == expected_alt_homozygous

    variant = "6:26199089:A:C"
    result = nomaly_genotype.query_variants(variant)
    assert result is not None

    genotypes = result[0]

    expected_missing_______ = 354
    expected_ref_homozygous = 487591
    expected_heterozygous__ = 5
    expected_alt_homozygous = 0

    assert np.sum(genotypes == -1) == expected_missing_______
    assert np.sum(genotypes == 0) == expected_ref_homozygous
    assert np.sum(genotypes == 1) == expected_heterozygous__
    assert np.sum(genotypes == 2) == expected_alt_homozygous


def test_genotype_variant_id_format_in_file(nomaly_genotype):
    """Test that variants in the file follow the expected format."""
    # Sample first 1000 variants
    sample_variants = nomaly_genotype.genotype_variant_id[:1000]

    for variant in sample_variants:
        parts = variant.split(":")
        assert len(parts) == 4, f"Invalid variant format: {variant}"
        assert parts[0].isdigit(), f"Invalid chromosome: {parts[0]}"
        assert parts[1].isdigit(), f"Invalid position: {parts[1]}"
        assert len(parts[2]) >= 1, f"Invalid ref allele: {parts[2]}"
        assert len(parts[3]) >= 1, f"Invalid alt allele: {parts[3]}"


def test_nomaly_variant_id_format_in_file(nomaly_genotype):
    """Test that variants in the file follow the expected format."""
    # Sample first 1000 variants
    sample_variants = nomaly_genotype.nomaly_variant_id[:1000]

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


def test_individual_id_format(nomaly_genotype):
    """Test that individual IDs are in the expected format."""
    sample_ids = nomaly_genotype.individual[:1000]

    # THIS WAS A PROBLEM IN THE ORIGINAL FILE, e..g there was exactly 1 -1 in
    # the sample_ids!

    NUM_MISSING_IDS = 0

    assert np.sum(sample_ids == -1) == NUM_MISSING_IDS, (
        f"There should be {NUM_MISSING_IDS} -1s in the sample_ids...ahhh"
    )

    assert np.sum(sample_ids < 0) == NUM_MISSING_IDS, (
        f"There should be {NUM_MISSING_IDS} -1s in the sample_ids...ahhh"
    )

    for id_num in sample_ids:
        assert isinstance(id_num, np.integer)
        if id_num == -1:
            continue
        assert id_num >= 100000, f"Individual ID: {id_num} is not greater than 100000"
        assert id_num < 10000000, f"Individual ID: {id_num} is not less than 10000000"
        assert len(str(id_num)) >= 5  # Example: Expecting at least 5-digit IDs


def test_missing_data_handling(nomaly_genotype):
    """Test handling of missing genotype data."""
    # Find a variant with some missing data (genotype = -1)
    variant = "11:69083946:T:C"
    result = nomaly_genotype.query_variants(variant)
    assert result is not None

    genotypes = result[0]
    missing_count = np.sum(genotypes == -1)
    assert missing_count < len(genotypes) * 0.1  # Less than 10% missing


def test_flipped_allele_query(nomaly_genotype):
    """Test querying variants with flipped alleles."""
    # Original variant
    variant = "8_6870776_C/T"
    # Expected flipped variant in file
    flipped = "8_6870776_T/C"

    # First verify the original variant isn't found
    result_original = nomaly_genotype.query_variants(variant)
    assert result_original is not None, (
        "Either original or flipped variant should be found"
    )

    # The query should automatically try the flipped version
    genotypes = result_original[0]
    assert all(g in [-1, 0, 1, 2] for g in genotypes)

    # Verify the flipped variant is also found
    result_flipped = nomaly_genotype.query_variants(flipped)
    assert result_flipped is not None, "Flipped variant should be found"


def test_variant_format_standardization(nomaly_genotype):
    """Test variant format standardization."""
    geno = nomaly_genotype

    test_cases = [
        # Format: (input, expected_output)
        ("8_6870776_C/T", "8:6870776:C:T"),
        ("8_6870776_C_T", "8:6870776:C:T"),
        ("8:6870776:C:T", "8:6870776:C:T"),  # Already standard
        # With chr prefix
        ("chr8_6870776_C/T", "8:6870776:C:T"),
        ("Chr8_6870776_C_T", "8:6870776:C:T"),
        ("CHR8:6870776:C:T", "8:6870776:C:T"),
    ]

    for input_variant, expected in test_cases:
        result = geno._standardize_variant_format(input_variant)
        assert result == expected, f"Failed for {input_variant}"

    test_cases = [
        # Invalid formats
        "invalid_format",
        "8_6870776",  # Missing alleles
        "8_pos_C/T",  # Invalid position
        "",  # Empty string
        "8_6870776_C",  # Missing alt
    ]

    # Test that invalid formats raise ValueError
    for invalid_variant in test_cases:
        with pytest.raises(ValueError):
            geno._standardize_variant_format(invalid_variant)


def test_query_with_different_formats(nomaly_genotype):
    """Test querying variants in different formats."""
    # Same variant in different formats
    variants = [
        "8_6870776_C/T",  # underscore format
        "8:6870776:C:T",  # colon format
        "chr8_6870776_C/T",  # with chr prefix
    ]

    # All should return the same data
    results = [nomaly_genotype.query_variants(v) for v in variants]
    results = [r for r in results if r is not None]

    # Check all results are identical
    assert len(results) > 0, "No variants found"
    for r in results[1:]:
        np.testing.assert_array_equal(
            results[0], r, "Different formats of same variant returned different data"
        )


def test_genotype_flipping(nomaly_genotype):
    """Test that genotypes are correctly flipped when alleles are flipped."""
    # Example variants where we know one is the flipped version of the other
    # These are placeholders - replace with actual variants from the data
    variant_ref = (
        "19:44908684:C:T"  # Replace with variant where we know genotype=0 means CC
    )
    variant_flipped = (
        "19:44908684:T:C"  # Same variant but flipped, so genotype=0 should mean TT
    )

    # Get genotypes for both orientations
    result_ref = nomaly_genotype.query_variants(variant_ref)
    result_flipped = nomaly_genotype.query_variants(variant_flipped)

    assert result_ref is not None, f"Variant {variant_ref} not found"
    assert result_flipped is not None, f"Variant {variant_flipped} not found"

    genotypes_ref = result_ref[0]
    genotypes_flipped = result_flipped[0]

    # Verify that homozygous ref (0) in one is homozygous alt (2) in the other
    ref_homozygous_mask = genotypes_ref == 0
    assert all(genotypes_flipped[ref_homozygous_mask] == 2), (
        "Homozygous ref not flipped to homozygous alt"
    )

    # Verify that homozygous alt (2) in one is homozygous ref (0) in the other
    alt_homozygous_mask = genotypes_ref == 2
    assert all(genotypes_flipped[alt_homozygous_mask] == 0), (
        "Homozygous alt not flipped to homozygous ref"
    )

    # Verify that heterozygous (1) stays heterozygous
    het_mask = genotypes_ref == 1
    assert all(genotypes_flipped[het_mask] == 1), (
        "Heterozygous genotypes changed during flipping"
    )

    # Verify that missing (-1) stays missing
    missing_mask = genotypes_ref == -1
    assert all(genotypes_flipped[missing_mask] == -1), (
        "Missing genotypes changed during flipping"
    )


def test_variant_counts(nomaly_genotype):
    """Test the variant counts for specific variants."""
    # Test a variant with known distribution
    variant = "19_44908684_T/C"
    counts = nomaly_genotype.get_variant_counts(variant)

    # Check the counts match expected values
    assert counts["missing"] == 73846
    assert counts["homozygous_ref"] == 9863
    assert counts["heterozygous"] == 108276
    assert counts["homozygous_alt"] == 295965
    assert counts["total"] == 487950

    # Test another variant with different distribution
    variant = "6_26199089_A/C"
    counts = nomaly_genotype.get_variant_counts(variant)

    assert counts["missing"] == 354
    assert counts["homozygous_ref"] == 0
    assert counts["heterozygous"] == 5
    assert counts["homozygous_alt"] == 487591
    assert counts["total"] == 487950

    # Test with ancestry filter
    counts_eur = nomaly_genotype.get_variant_counts(variant, ancestry="EUR")
    assert counts_eur["missing"] == 319
    assert counts_eur["homozygous_ref"] == 0
    assert counts_eur["heterozygous"] == 5
    assert counts_eur["homozygous_alt"] == 449099
    assert counts_eur["total"] == 449423

    assert sum(counts_eur.values()) <= sum(counts.values())

    # Test with sex filter
    counts_female = nomaly_genotype.get_variant_counts(variant, sex="F")
    assert isinstance(counts_female["missing"], int)
    assert isinstance(counts_female["homozygous_ref"], int)
    assert isinstance(counts_female["heterozygous"], int)
    assert isinstance(counts_female["homozygous_alt"], int)
    assert isinstance(counts_female["total"], int)
    assert sum(counts_female.values()) <= sum(counts.values())
