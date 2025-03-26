"""Data services package providing access to genomic and phenotypic data.

NOTE(S):
    - Without init imports (below), users would need to do:
      from data_services.stats import StatsService

    - With init imports, they can do:
      from data_services import StatsService

    - We could also set up 'package level logging' here.
"""

from .registry import ServiceRegistry

from .genotype import GenotypeService
from .phenotype import PhenotypeService
from .nomaly_data import NomalyDataService
from .nomaly_score import NomalyScoreService

from .stats import StatsRegistry, StatsService

__all__ = [
    "ServiceRegistry",
    "GenotypeService",
    "PhenotypeService",
    "NomalyDataService",
    "NomalyScoreService",
    "StatsRegistry",
    "StatsService",
]
