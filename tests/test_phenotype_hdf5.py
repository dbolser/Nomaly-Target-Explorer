"""This is an DATA TEST, where we look at the real production data."""

import pytest
import numpy as np
import h5py
import tempfile

from data_services.phenotype import PhenotypeService

from config import Config


# By setting the scope to "module", we only read the HDF5 file once for all the
# following tests.. I think.
@pytest.fixture(scope="module")
def production_phenotype_service():
    """Fixture to create a PhenotypeService instance."""
    return PhenotypeService(hdf5_file=Config.PHENOTYPES_HDF)


# TODO: Move this to a test for the HDF5 file itself
def test_get_cases_invalid_data():
    """Test handling of invalid data in the HDF5 file."""
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        with h5py.File(tmp.name, "w") as f:
            # Create mismatched data
            f.create_dataset(
                "phenotype_data",
                data=np.array([[0], [0]], dtype=np.int8),  # Only 2 individuals
            )
            f.create_dataset("eids", data=np.array([101001, 101002, 101003]))
            f.create_dataset("phecodes", data=np.array([b"777.7"]))
            f.create_dataset("sex", data=np.array([b"M", b"F", b"M"]))
            f.create_dataset("populations", data=np.array([b"EUR", b"EUR", b"SAS"]))

            with pytest.raises(KeyError):
                _ = PhenotypeService(tmp.name)

