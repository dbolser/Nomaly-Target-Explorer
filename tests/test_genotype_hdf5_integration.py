import numpy as np

from blueprints.nomaly import GenotypeHDF5

from config import Config

import pytest

config = Config()
nomaly_genotype = GenotypeHDF5(config.GENOTYPES_H5)

NUM_INDIVIDUALS = 488377
NUM_INDIVIDUALS = 487950  # Removed 427 individuals with negative eids
NUM_VARIANTS = 83011


def test_production_file_exists():
    """Verify the production genotype file exists and is readable."""
    assert nomaly_genotype is not None
    assert hasattr(nomaly_genotype, "f")
    assert "genotype_matrix" in nomaly_genotype.f


def test_known_variant_query():
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


def test_individual_count():
    """Verify the expected number of individuals in the dataset."""
    assert len(nomaly_genotype.individual) == NUM_INDIVIDUALS


def test_variant_count():
    """Verify the expected number of variants in the dataset."""
    assert len(nomaly_genotype.variants) == NUM_VARIANTS


def test_genotype_matrix_shape():
    """Verify the shape of the genotype matrix."""
    assert nomaly_genotype.genotype_matrix.shape == (NUM_VARIANTS, NUM_INDIVIDUALS)


def test_known_variant_genotype_distribution():
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


def test_variant_format_in_file():
    """Test that variants in the file follow the expected format."""
    # Sample first 1000 variants
    sample_variants = nomaly_genotype.variants[:1000]

    for variant in sample_variants:
        parts = variant.split(":")
        assert len(parts) == 4, f"Invalid variant format: {variant}"
        assert parts[0].isdigit(), f"Invalid chromosome: {parts[0]}"
        assert parts[1].isdigit(), f"Invalid position: {parts[1]}"
        assert len(parts[2]) >= 1, f"Invalid ref allele: {parts[2]}"
        assert len(parts[3]) >= 1, f"Invalid alt allele: {parts[3]}"


def test_individual_id_format():
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


def test_missing_data_handling():
    """Test handling of missing genotype data."""
    # Find a variant with some missing data (genotype = -1)
    variant = "11:69083946:T:C"
    result = nomaly_genotype.query_variants(variant)
    assert result is not None

    genotypes = result[0]
    missing_count = np.sum(genotypes == -1)
    assert missing_count < len(genotypes) * 0.1  # Less than 10% missing


def test_flipped_allele_query():
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


def test_variant_format_standardization():
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


def test_query_with_different_formats():
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


def test_genotype_flipping():
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


def test_variant_counts():
    """Test the variant counts."""
    counts = nomaly_genotype.get_variant_counts()

    # Verify we have the expected number of variants
    assert len(counts) == NUM_VARIANTS

    # Test specific variant counts using DataFrame indexing
    variant1 = "19:44908684:C:T"
    assert counts.loc[variant1, "heterozygous"] == 108371
    assert counts.loc[variant1, "homozygous_alt"] == 296223
    assert counts.loc[variant1, "homozygous_ref"] == 9874
    assert counts.loc[variant1, "missing"] == 73909

    variant2 = "11:69083946:T:C"
    assert counts.loc[variant2, "heterozygous"] == 11499
    assert counts.loc[variant2, "homozygous_alt"] == 474918
    assert counts.loc[variant2, "homozygous_ref"] == 1516
    assert counts.loc[variant2, "missing"] == 444
