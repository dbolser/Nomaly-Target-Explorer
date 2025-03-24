from typing import Optional

from data_services.genotype import GenotypeService
from data_services.nomaly_data import NomalyDataService
from data_services.nomaly_score import NomalyScoreService
from data_services.phenotype import PhenotypeService
from data_services.stats import StatsRegistry, StatsService


class ServiceRegistry:
    def __init__(self, app=None):
        self.genotype = GenotypeService()
        self.phenotype = PhenotypeService()
        self.stats_registry: StatsRegistry = StatsRegistry()
        self.nomaly_score: Optional[NomalyScoreService] = None
        self.nomaly_data: Optional[NomalyDataService] = None

        # If TESTING is True, we don't want to initialise the data services
        if app is not None and not app.config.get("TESTING"):
            self.init_app(app)

    def init_app(self, app):
        """Initialize services with Flask app"""
        if not hasattr(app, "extensions"):
            app.extensions = {}

        app.extensions["nomaly_services"] = self

        # Initialize services from config
        self.genotype = GenotypeService(app.config.get("GENOTYPES_H5"))
        self.nomaly_data = NomalyDataService(app.config.get("NOMALY_VARIANTS_PATH"))
        self.nomaly_score = NomalyScoreService(app.config.get("NOMALY_SCORES_H5"))
        self.phenotype = PhenotypeService(app.config.get("PHENOTYPES_H5"))
        # self.stats = StatsService(app.config.get("STATS_H5"))

        # Initialize the stats registry with the stats selector from config
        if app.config.get("STATS_SELECTOR"):
            self.stats_registry = StatsRegistry(app.config.get("STATS_SELECTOR"))

    def init_from_config(self, config):
        self.genotype = GenotypeService(config.GENOTYPES_H5)
        self.nomaly_data = NomalyDataService(config.NOMALY_VARIANTS_PATH)
        self.nomaly_score = NomalyScoreService(config.NOMALY_SCORES_H5)
        self.phenotype = PhenotypeService(config.PHENOTYPES_H5)
        # self.stats = StatsService(config.STATS_H5)
        self.stats_registry = StatsRegistry(config.STATS_SELECTOR)

    @classmethod
    def from_config(cls, config):
        registry = cls()
        registry.init_from_config(config)

        return registry
