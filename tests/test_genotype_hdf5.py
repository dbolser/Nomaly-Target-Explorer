import pytest
import numpy as np
import h5py
import os
import tempfile
from pathlib import Path
from data_services.genotype import GenotypesHDF5


def test_init_with_valid_file(mock_genotype_hdf5_file_with_npy):
    """Test initialization with a valid HDF5 file."""
    geno = GenotypesHDF5(mock_genotype_hdf5_file_with_npy)
    assert geno is not None
    assert hasattr(geno, "f")
    assert hasattr(geno, "genotype_matrix")
    assert hasattr(geno, "genotype_matrix_mm")
    assert hasattr(geno, "genotype_matrix_h5")
    assert hasattr(geno, "individual")
    assert hasattr(geno, "genotype_variant_id")
    assert hasattr(geno, "nomaly_variant_id")
    assert hasattr(geno, "ancestry")
    assert hasattr(geno, "biological_sex")


def test_init_with_nonexistent_file():
    """Test initialization with a non-existent file."""
    with pytest.raises(FileNotFoundError):
        GenotypesHDF5("nonexistent.h5")


def test_init_with_invalid_file():
    """Test initialization with an invalid HDF5 file."""
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create file without required datasets
            f.create_dataset("some_other_data", data=[1, 2, 3])

        with pytest.raises(KeyError):
            GenotypesHDF5(Path(tmp.name))

        os.unlink(tmp.name)


def test_query_variants_single(mock_genotype_hdf5_file_with_npy):
    """Test querying a single variant."""
    geno = GenotypesHDF5(mock_genotype_hdf5_file_with_npy)
    result = geno.query_variants("1:100:A:T")
    assert result is not None
    assert isinstance(result, np.ndarray)
    assert result.shape == (1, 4)  # One variant, four individuals
    np.testing.assert_array_equal(result[0], [0, 1, 1, 2])


def test_query_variants_multiple(mock_genotype_hdf5_file_with_npy):
    """Test querying multiple variants."""
    geno = GenotypesHDF5(mock_genotype_hdf5_file_with_npy)
    result = [geno.query_variants(variant) for variant in ["1:100:A:T", "2:200:C:G"]]

    for r in result:
        assert r is not None

    # Make a 2D array of results (ugly...)
    result = [x[0] for x in result]
    result = np.array(result)

    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 4)  # Two variants, four individuals
    np.testing.assert_array_equal(result, [[0, 1, 1, 2], [1, 0, 0, 1]])


def test_query_variants_nonexistent(mock_genotype_hdf5_file_with_npy):
    """Test querying a non-existent variant."""
    geno = GenotypesHDF5(mock_genotype_hdf5_file_with_npy)
    result = geno.query_variants("3:300:G:C")
    assert result.size == 0


def test_query_variants_invalid_format(mock_genotype_hdf5_file_with_npy):
    """Test querying with invalid variant format."""
    geno = GenotypesHDF5(mock_genotype_hdf5_file_with_npy)
    result = geno.query_variants("invalid_format")
    assert result.size == 0


def test_query_variants_empty_input(mock_genotype_hdf5_file_with_npy):
    """Test querying with empty input."""
    geno = GenotypesHDF5(mock_genotype_hdf5_file_with_npy)
    result = geno.query_variants("")
    assert result.size == 0


def test_single_variant_mask(mock_genotype_hdf5_file_with_npy):
    """Test creation of single variant mask."""
    geno = GenotypesHDF5(mock_genotype_hdf5_file_with_npy)
    mask = geno._single_variant_mask("1:100:A:T")
    assert isinstance(mask, np.ndarray)
    assert mask.dtype == bool
    assert np.sum(mask) == 1
    assert mask[0]  # First variant should match


def test_large_dataset_handling():
    """Test handling of larger datasets."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create larger datasets
            genotype_data = np.random.randint(0, 3, (1000, 1000), dtype=np.int8)
            f.create_dataset("genotype_matrix", data=genotype_data)
            f.create_dataset("fam", data=np.arange(1000))
            f.create_dataset(
                "sex", data=np.random.choice(["M", "F"], 1000).astype(np.string_)
            )
            f.create_dataset(
                "ancestry",
                data=np.random.choice(["EUR", "SAS"], 1000).astype(np.string_),
            )
            f.create_dataset("bim", data=[f"{i}:100:A:T".encode() for i in range(1000)])
            f.create_dataset(
                "plink_variant_id", data=[f"{i}:100:A:T".encode() for i in range(1000)]
            )
            f.create_dataset(
                "nomaly_variant_id", data=[f"{i}:100:A:T".encode() for i in range(1000)]
            )

        # Create the corresponding .npy file
        np.save(f"{tmp.name}.npy", genotype_data)

        try:
            geno = GenotypesHDF5(Path(tmp.name))
            result = geno.query_variants("0:100:A:T")
            assert result is not None
            assert result.shape == (1, 1000)
        finally:
            # Clean up the .npy file
            if Path(f"{tmp.name}.npy").exists():
                Path(f"{tmp.name}.npy").unlink()


def test_error_handling_corrupted_data():
    """Test handling of corrupted data with mismatched dimensions."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create intentionally mismatched datasets
            f.create_dataset(
                "genotype_matrix", data=np.array([[0, 1]], dtype=np.int8)
            )  # 1x2 matrix
            f.create_dataset("fam", data=np.array([1001, 1002, 1003]))  # 3 individuals
            f.create_dataset("sex", data=np.array(["M", "F", "M"]).astype(np.string_))
            f.create_dataset(
                "ancestry", data=np.array(["EUR", "EUR", "SAS"]).astype(np.string_)
            )
            f.create_dataset("bim", data=np.array([b"1:100:A:T"]))  # 1 variant
            f.create_dataset("nomaly_variant_id", data=np.array([b"1:100:A:T"]))
            f.create_dataset("plink_variant_id", data=np.array([b"1:100:A:T"]))
        with pytest.raises(AssertionError) as excinfo:
            _ = GenotypesHDF5(tmp.name)
            assert "Error in sanity checks" in str(excinfo.value)


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
def test_variant_format_validation(mock_genotype_hdf5_file_with_npy, variant_id):
    """Test validation of various variant ID formats."""
    geno = GenotypesHDF5(mock_genotype_hdf5_file_with_npy)
    if variant_id == "1:100:A:T":
        result = geno.query_variants(variant_id)
        assert result.size > 0
    else:
        result = geno.query_variants(variant_id)
        assert result.size == 0
