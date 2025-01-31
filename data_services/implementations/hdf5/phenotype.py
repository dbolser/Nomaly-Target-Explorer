from os import PathLike
from typing import Tuple
import numpy as np

from data_services.interfaces.phenotype import PhenotypeService
from blueprints.nomaly import PhenotypesHDF5


class HDF5PhenotypeService(PhenotypeService):
    def __init__(self, hdf5_file: PathLike):
        self._hdf = PhenotypesHDF5(hdf5_file)  # Existing implementation

    def get_cases_for_phecode(
        self, phecode: str, sex: str | None = None, population: str | None = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        return self._hdf.get_cases_for_phecode(phecode, population)
