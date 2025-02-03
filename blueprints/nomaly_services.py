from blueprints.nomaly import StatsHDF5
from config import Config
from typing import Optional

class NomalyServices:
    def __init__(self):
        self.genotype = None
        self.phenotype = None
        self.stats = None
        self.stats_v2 = None
        self._initialized = False

    def init_from_config(self, config: Optional[dict] = None):
        """Initialize services directly from config"""
        if self._initialized:
            return

        # Use provided config or fall back to default Config
        config = config or Config.__dict__
        
        # Initialize services
        self.stats = StatsHDF5(config.get("STATS_H5", Config.STATS_H5))
        self.stats_v2 = StatsHDF5(config.get("STATS_H5_V2", Config.STATS_H5_V2))
        self._initialized = True

    def init_app(self, app):
        """Initialize with Flask app (maintaining backward compatibility)"""
        self.init_from_config(app.config)

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
