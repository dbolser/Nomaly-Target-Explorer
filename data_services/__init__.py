"""Data services package providing access to genomic and phenotypic data.

NOTE(S):
    - Without init imports (below), users would need to do:
      from data_services.stats import StatsService

    - With init imports, they can do:
      from data_services import StatsService

    - We could also set up 'package level logging' here.
"""

from .genotype import GenotypeService
from .nomaly_score import NomalyScoreService
from .phenotype import PhenotypeService
from .registry import ServiceRegistry
from .stats import StatsService, StatsRegistry

__all__ = [
    "ServiceRegistry",
    "GenotypeService",
    "PhenotypeService",
    "NomalyScoreService",
    "StatsService",
    "StatsRegistry",
]
