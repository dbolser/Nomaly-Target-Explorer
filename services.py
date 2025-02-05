from blueprints.nomaly import StatsHDF5
from data_services.genotype import GenotypeService
from data_services.phenotype import PhenotypeService
from data_services.nomaly_score import NomalyScoreService


class ServiceRegistry:
    def __init__(self, app=None):
        self.genotype = None
        self.phenotype = None
        self.stats = None
        self.stats_v2 = None
        self.nomaly_score = None
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
        self.nomaly_score = NomalyScoreService(app.config.get("NOMALY_SCORES_H5"))
        self.stats = StatsHDF5(app.config.get("STATS_H5"))
        self.stats_v2 = StatsHDF5(app.config.get("STATS_H5_V2"))
