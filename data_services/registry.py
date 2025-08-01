from config import Config

from data_services.genotype import GenotypeService
from data_services.nomaly_data import NomalyDataService
from data_services.nomaly_score import NomalyScoreService
from data_services.phenotype import PhenotypeService
from data_services.stats import StatsRegistry


class ServiceRegistry:
    def __init__(self, app=None):
        self.genotype = GenotypeService()
        self.phenotype = PhenotypeService()
        self.nomaly_data = NomalyDataService()
        self.nomaly_score = NomalyScoreService()

        self.stats_registry = StatsRegistry()

        # If TESTING is True, we don't want to initialise the data services
        if app is not None and not app.config.get("TESTING"):
            self.init_app(app)

    def init_app(self, app):
        """Initialize services with Flask app"""
        if not hasattr(app, "extensions"):
            app.extensions = {}

        app.extensions["nomaly_services"] = self

        # Initialize services from app config
        self.genotype = GenotypeService(app.config.get("GENOTYPES_HDF"))
        self.phenotype = PhenotypeService(app.config.get("PHENOTYPES_HDF"))
        self.nomaly_data = NomalyDataService(
            app.config.get("NOMALY_VARIANT_MAPPING_PATH")
        )
        self.nomaly_score = NomalyScoreService(app.config.get("NOMALY_SCORES_H5"))

        # Initialize the stats registry with the stats selector from config
        if app.config.get("STATS_SELECTOR"):
            self.stats_registry = StatsRegistry(app.config.get("STATS_SELECTOR"))

    def init_from_config(self, config: Config):
        self.genotype = GenotypeService(config.GENOTYPES_HDF)
        self.phenotype = PhenotypeService(config.PHENOTYPES_HDF)
        self.nomaly_data = NomalyDataService(config.NOMALY_VARIANT_MAPPING_PATH)
        self.nomaly_score = NomalyScoreService(config.NOMALY_SCORES_H5)

        self.stats_registry = StatsRegistry(config.STATS_SELECTOR)

    @classmethod
    def from_config(cls, config):
        registry = cls()
        registry.init_from_config(config)

        return registry
