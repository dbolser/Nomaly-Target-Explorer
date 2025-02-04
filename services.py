from blueprints.nomaly import StatsHDF5
from data_services.genotype import GenotypeService
from data_services.phenotype import PhenotypeService


class ServiceRegistry:
    def __init__(self, app=None):
        self.genotype: GenotypeService | None = None
        self.phenotype: PhenotypeService | None = None
        self.stats: StatsHDF5 | None = None
        self.stats_v2: StatsHDF5 | None = None

        if app is not None:
            self.init_app(app)

    def init_app(self, app):
        """Initialize services with Flask app"""
        if not hasattr(app, "extensions"):
            app.extensions = {}

        app.extensions["nomaly_services"] = self

        # Initialize services from config
        self.genotype = GenotypeService(app.config.get("GENOTYPES_H5"))
        self.phenotype = PhenotypeService(app.config.get("PHENOTYPES_H5"))
        self.stats = StatsHDF5(app.config.get("STATS_H5"))
        self.stats_v2 = StatsHDF5(app.config.get("STATS_H5_V2"))


services = ServiceRegistry()
