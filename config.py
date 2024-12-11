import os
from pathlib import Path

from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()


class Config:
    """Base configuration"""

    # Flask settings
    HOST = "0.0.0.0"
    PORT = 8756
    DEBUG = False

    # Session settings
    SESSION_TYPE = "filesystem"
    SECRET_KEY = os.getenv("FLASK_SECRET_KEY")

    # Database settings
    MYSQL_HOST = os.getenv("MYSQL_HOST")
    MYSQL_USER = os.getenv("MYSQL_USER")
    MYSQL_PASSWORD = os.getenv("MYSQL_PASSWORD")
    MYSQL_DB = os.getenv("MYSQL_DB")

    # Application paths
    NOMALY_RESULTS_DIR = Path("/data/general/UKBB/Run-v2/DatabaseInputs")
    PHEWAS_PHENO_DIR = Path("/data/clu/ukbb/by_variant")
    UKBB_PHENO_DIR = Path("/data/general/UKBB/Phenotypes")
    GWAS_PHENO_DIR = Path("/data/clu/ukbb/by_pheno")

    PHENOTYPES_H5 = UKBB_PHENO_DIR / "phecode_cases_excludes.h5"
    GENOTYPES_H5 = NOMALY_RESULTS_DIR / "genotypes.h5"

class DevelopmentConfig(Config):
    """Development configuration"""
    DEBUG = True


class ProductionConfig(Config):
    """Production configuration"""
    pass


class TestingConfig(Config):
    """Testing configuration"""
    TESTING = True
    DEBUG = True


# Map environment names to config objects
config = {
    "development": DevelopmentConfig,
    "production": ProductionConfig,
    "testing": TestingConfig,
    "default": DevelopmentConfig,
}
