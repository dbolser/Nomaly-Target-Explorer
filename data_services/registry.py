from data_services.genotype import GenotypeService
from data_services.phenotype import PhenotypeService
from data_services.nomaly_score import NomalyScoreService
from data_services.stats import StatsService


class ServiceRegistry:
    def __init__(self, app=None):
        self.genotype = None
        self.phenotype = None
        self.stats = None
        self.nomaly_score = None

        # Don't ask...
        self.stats_v2 = None
        self.nomaly_score_v2 = None

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
        self.phenotype = PhenotypeService(app.config.get("PHENOTYPES_H5"))
        self.stats = StatsService(app.config.get("STATS_H5"))
        self.nomaly_score = NomalyScoreService(app.config.get("NOMALY_SCORES_H5"))

        # Please don't ask
        self.stats_v2 = StatsService(app.config.get("STATS_H5_V2"))
        self.nomaly_score_v2 = NomalyScoreService(app.config.get("NOMALY_SCORES_H5_V2"))

    def init_from_config(self, config):
        self.genotype = GenotypeService(config.GENOTYPES_H5)
        self.phenotype = PhenotypeService(config.PHENOTYPES_H5)
        self.stats = StatsService(config.STATS_H5)
        self.nomaly_score = NomalyScoreService(config.NOMALY_SCORES_H5)

    @classmethod
    def from_config(cls, config):
        registry = cls()
        registry.init_from_config(config)

        return registry
