from abc import ABC, abstractmethod
from typing import Tuple
import numpy as np


class PhenotypeService(ABC):
    @abstractmethod
    def get_cases_for_phecode(
        self, phecode: str, sex: str | None = None, population: str | None = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Get cases and controls for a phecode.

        Returns:
            Tuple[np.ndarray, np.ndarray]: (eids, case_status)
            where case_status is 1 for cases, 0 for controls
        """
        pass
