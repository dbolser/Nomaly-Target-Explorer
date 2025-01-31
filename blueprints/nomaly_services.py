import os

from blueprints.nomaly import GenotypeHDF5, PhenotypesHDF5, StatsHDF5
from config import Config


class NomalyServices:
    def __init__(self):
        self.genotype = None
        self.phenotype = None
        self.stats = None
        self.stats_v2 = None

    def init_app(self, app):
        """Initialize with Flask app (exactly once)"""
        if hasattr(self, "_initialized") and self._initialized:
            return

        self.genotype = GenotypeHDF5(Config.GENOTYPES_H5)
        self.phenotype = PhenotypesHDF5(Config.PHENOTYPES_H5)
        self.stats = StatsHDF5(Config.STATS_H5)
        self.stats_v2 = StatsHDF5(Config.STATS_H5_V2)
        self._initialized = True

        if not hasattr(app, "extensions"):
            app.extensions = {}
        app.extensions["nomaly_services"] = self


# Single instance
if not os.environ.get("TESTING"):
    services = NomalyServices()
