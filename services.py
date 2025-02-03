from typing import Optional

from blueprints.nomaly import StatsHDF5
from config import Config
from data_services.genotype import GenotypeService
from data_services.phenotype import PhenotypeService


class ServiceRegistry:
    def __init__(self):
        self.genotype: GenotypeService | None = None
        self.phenotype: PhenotypeService | None = None
        self.stats: StatsHDF5 | None = None
        self.stats_v2: StatsHDF5 | None = None
        self._initialized = False

    @classmethod
    def create_from_config(cls, config_dict):
        """Create and return a new initialized instance from config dict."""
        instance = cls()
        instance.init_from_config(config_dict)
        return instance

    def init_from_app(self, app):
        """Initialize with Flask app (maintaining backward compatibility)."""
        self.init_from_config(app.config)

        if not hasattr(app, "extensions"):
            app.extensions = {}
        app.extensions["nomaly_services"] = self

    def init_from_config(self, config: Optional[dict] = None):
        """Initialize services directly from config"""
        if self._initialized:
            return

        # Use provided config or fall back to default Config
        config = config or Config.__dict__

        # Initialize services
        self.genotype = GenotypeService(config.get("GENOTYPES_H5", Config.GENOTYPES_H5))
        self.phenotype = PhenotypeService(
            config.get("PHENOTYPES_H5", Config.PHENOTYPES_H5)
        )
        self.stats = StatsHDF5(config.get("STATS_H5", Config.STATS_H5))
        self.stats_v2 = StatsHDF5(config.get("STATS_H5_V2", Config.STATS_H5_V2))
        self._initialized = True


# Global instance - keeping for backward compatibility
services = ServiceRegistry()
