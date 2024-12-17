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

    # Caching directories
    PHEWAS_PHENO_DIR = Path("/data/clu/ukbb/by_variant")
    GWAS_PHENO_DIR = Path("/data/clu/ukbb/by_pheno")

    # Phenotype directories and files
    UKBB_PHENO_DIR = Path("/data/general/UKBB/Phenotypes")
    PHENOTYPES_H5 = UKBB_PHENO_DIR / "phecode_cases_excludes.h5"

    # Nomaly results directories (V1)
    NOMALY_RESULTS_DIR = Path("/data/general/UKBB/Run-v1/DatabaseInputs")

    # Genotype files
    GENOTYPES_H5 = NOMALY_RESULTS_DIR / "genotypes_with_counts.h5"

    # Stats and scores files (V1)
    STATS_H5 = NOMALY_RESULTS_DIR / "stats.h5"
    # SCORES_H5 = NOMALY_RESULTS_DIR / "float16_scores.h5"

    # Nomaly results directories (V2)
    NOMALY_RESULTS_DIR_V2 = Path("/data/general/UKBB/Run-v2/DatabaseInputs")
    STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats.h5"
    # SCORES_H5_V2 = NOMALY_RESULTS_DIR_V2 / "float16_scores.h5"


    RESOURCE_DATA_DIR = Path("/data/general/Data/")


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
