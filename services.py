from data_services.interfaces.phenotype import PhenotypeService
from data_services.implementations.hdf5.phenotype import HDF5PhenotypeService


class ServiceRegistry:
    def __init__(self):
        self._phenotype_service: PhenotypeService | None = None

    def init_app(self, app):
        if hasattr(self, "_initialized"):
            return

        self._phenotype_service = HDF5PhenotypeService(app.config["PHENOTYPES_H5"])
        self._initialized = True

    @property
    def phenotype(self) -> PhenotypeService:
        if self._phenotype_service is None:
            raise RuntimeError("Phenotype service not initialized")
        return self._phenotype_service

    @phenotype.setter
    def phenotype(self, service: PhenotypeService):
        self._phenotype_service = service


# Global instance
services = ServiceRegistry()
