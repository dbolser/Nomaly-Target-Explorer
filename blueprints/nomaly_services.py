from blueprints.nomaly import GenotypeHDF5, StatsHDF5
from config import Config

from services import services as new_services


class NomalyServices:
    def __init__(self):
        self.genotype = None
        self.phenotype = None
        self.stats = None
        self.stats_v2 = None
        self._initialized = False

    def init_app(self, app):
        """Initialize with Flask app (exactly once)"""
        if self._initialized:
            return

        # Initialize real services
        config = app.config
        self.genotype = GenotypeHDF5(config.get("GENOTYPES_H5", Config.GENOTYPES_H5))
        self.phenotype = new_services.phenotype
        self.stats = StatsHDF5(config.get("STATS_H5", Config.STATS_H5))
        self.stats_v2 = StatsHDF5(config.get("STATS_H5_V2", Config.STATS_H5_V2))

        self._initialized = True

        if not hasattr(app, "extensions"):
            app.extensions = {}
        app.extensions["nomaly_services"] = self

    def reset(self):
        """Reset the service state - useful for testing"""
        self._initialized = False
        self.genotype = None
        self.phenotype = None
        self.stats = None
        self.stats_v2 = None


# Single instance - create it but don't initialize until explicitly requested
services = NomalyServices()
