import pytest
import numpy as np
import h5py
import os
import tempfile

from blueprints.nomaly import GenotypeHDF5


@pytest.fixture
def mock_genotype_file():
    """Create a temporary mock HDF5 file with test data."""
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create required datasets
            f.create_dataset(
                "genotype_matrix",
                data=np.array(
                    [
                        [0, 1, 2],  # Genotypes for first variant
                        [1, 0, 2],  # Genotypes for second variant
                    ]
                ),
            )
            f.create_dataset("fam", data=np.array([1001, 1002, 1003]))  # Individual IDs
            f.create_dataset(
                "bim",
                data=np.array(
                    [
                        b"1:100:A:T",  # First variant
                        b"2:200:C:G",  # Second variant
                    ]
                ),
            )
        yield tmp.name
    os.unlink(tmp.name)


def test_init_with_valid_file(mock_genotype_file):
    """Test initialization with a valid HDF5 file."""
    geno = GenotypeHDF5(mock_genotype_file)
    assert geno is not None
    assert hasattr(geno, "f")
    assert hasattr(geno, "genotype_matrix")
    assert hasattr(geno, "individual")
    assert hasattr(geno, "variants")


def test_init_with_nonexistent_file():
    """Test initialization with a non-existent file."""
    with pytest.raises(FileNotFoundError):
        GenotypeHDF5("nonexistent.h5")


def test_init_with_invalid_file():
    """Test initialization with an invalid HDF5 file."""
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create file without required datasets
            f.create_dataset("some_other_data", data=[1, 2, 3])

        with pytest.raises(KeyError):
            GenotypeHDF5(tmp.name)

        os.unlink(tmp.name)


def test_query_variants_single(mock_genotype_file):
    """Test querying a single variant."""
    geno = GenotypeHDF5(mock_genotype_file)
    result = geno.query_variants("1:100:A:T")
    assert result is not None
    assert isinstance(result, np.ndarray)
    assert result.shape == (1, 3)  # One variant, three individuals
    np.testing.assert_array_equal(result[0], [0, 1, 2])


def test_query_variants_multiple(mock_genotype_file):
    """Test querying multiple variants."""
    geno = GenotypeHDF5(mock_genotype_file)
    result = [geno.query_variants(variant) for variant in ["1:100:A:T", "2:200:C:G"]]

    for r in result:
        assert r is not None

    # Make a 2D array of results (ugly...)
    result = [x[0] for x in result]
    result = np.array(result)

    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 3)  # Two variants, three individuals
    np.testing.assert_array_equal(result, [[0, 1, 2], [1, 0, 2]])


def test_query_variants_nonexistent(mock_genotype_file):
    """Test querying a non-existent variant."""
    geno = GenotypeHDF5(mock_genotype_file)
    result = geno.query_variants("3:300:G:C")
    assert result is None


def test_query_variants_invalid_format(mock_genotype_file):
    """Test querying with invalid variant format."""
    geno = GenotypeHDF5(mock_genotype_file)
    result = geno.query_variants("invalid_format")
    assert result is None


def test_query_variants_empty_input(mock_genotype_file):
    """Test querying with empty input."""
    geno = GenotypeHDF5(mock_genotype_file)
    result = geno.query_variants("")
    assert result is None
    result = geno.query_variants([])
    assert result is None


def test_single_variant_mask(mock_genotype_file):
    """Test creation of single variant mask."""
    geno = GenotypeHDF5(mock_genotype_file)
    mask = geno._single_variant_mask("1:100:A:T")
    assert isinstance(mask, np.ndarray)
    assert mask.dtype == bool
    assert np.sum(mask) == 1
    assert mask[0]  # First variant should match


def test_large_dataset_handling(mock_genotype_file):
    """Test handling of larger datasets."""
    with h5py.File(mock_genotype_file, "w") as f:
        # Create larger datasets
        f.create_dataset("genotype_matrix", data=np.random.randint(0, 3, (1000, 1000)))
        f.create_dataset("fam", data=np.arange(1000))
        f.create_dataset("bim", data=[f"{i}:100:A:T".encode() for i in range(1000)])

    geno = GenotypeHDF5(mock_genotype_file)
    result = geno.query_variants("0:100:A:T")
    assert result is not None
    assert result.shape == (1, 1000)


def test_error_handling_corrupted_data(mock_genotype_file):
    """Test handling of corrupted data."""
    with h5py.File(mock_genotype_file, "w") as f:
        # Create mismatched datasets
        f.create_dataset("genotype_matrix", data=np.array([[0, 1]]))  # 2 individuals
        f.create_dataset("fam", data=np.array([1001, 1002, 1003]))  # 3 individuals
        f.create_dataset("bim", data=np.array([b"1:100:A:T"]))

    geno = GenotypeHDF5(mock_genotype_file)
    # This should still initialize, but querying might fail
    result = geno.query_variants("1:100:A:T")
    # The implementation should either return None or raise an error
    # depending on how we want to handle this case
    assert result is None or isinstance(result, np.ndarray)


@pytest.mark.parametrize(
    "variant_id",
    [
        "1:100:A:T",  # Valid format
        "1:100:A:T:extra",  # Too many parts
        "1:100:A",  # Too few parts
        "",  # Empty string
        None,  # None
        123,  # Wrong type
    ],
)
def test_variant_format_validation(mock_genotype_file, variant_id):
    """Test validation of various variant ID formats."""
    geno = GenotypeHDF5(mock_genotype_file)
    if variant_id == "1:100:A:T":
        result = geno.query_variants(variant_id)
        assert result is not None
    else:
        result = geno.query_variants(variant_id)
        assert result is None
