"""
Tests for the GenotypesHDF5 class.

Testing Strategy:
---------------
This file contains only basic initialization tests for the GenotypesHDF5 class.
Functional tests have been moved to test_genotype_hdf5_integration.py, which
uses real data to ensure more reliable testing of the actual functionality.

The separation allows us to: 1. Keep lightweight initialization tests for basic
setup and error handling 2. Use comprehensive integration tests for functional
verification 3. Reduce maintenance burden by eliminating duplicate tests 4.
Avoid issues with mock data that may not accurately represent real data
"""

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
    # assert hasattr(geno, "genotype_matrix_mm")
    # assert hasattr(geno, "genotype_matrix_h5")
    assert hasattr(geno, "individual")
    assert hasattr(geno, "plink_variant_id")
    assert hasattr(geno, "nomaly_variant_id")
    assert hasattr(geno, "ancestry")
    assert hasattr(geno, "individual_sex")


def test_init_with_nonexistent_file():
    """Test initialization with a non-existent file."""
    with pytest.raises(FileNotFoundError):
        GenotypesHDF5(Path("nonexistent.h5"))


def test_init_with_invalid_file():
    """Test initialization with an invalid HDF5 file."""
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create file without required datasets
            f.create_dataset("some_other_data", data=[1, 2, 3])

        with pytest.raises(KeyError):
            GenotypesHDF5(Path(tmp.name))

        os.unlink(tmp.name)

# The following functional tests have been removed as they are covered by the comprehensive
# integration tests in test_genotype_hdf5_integration.py
