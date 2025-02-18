import os
from pathlib import Path

from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()


class Config:
    """Base configuration"""

    # FLASK SETTINGS

    # Session settings (common across all environments)
    SESSION_TYPE = "filesystem"
    SESSION_PROTECTION = "basic"
    SESSION_COOKIE_HTTPONLY = True
    SESSION_COOKIE_SECURE = True  # HTTPS only by default
    SESSION_COOKIE_SAMESITE = "Lax"

    # NOMALY SETTINGS

    # Database settings from environment variables
    MYSQL_HOST = os.getenv("MYSQL_HOST")
    MYSQL_PORT = os.getenv("MYSQL_PORT", 3306)
    MYSQL_USER = os.getenv("MYSQL_USER")
    MYSQL_PASSWORD = os.getenv("MYSQL_PASSWORD")
    MYSQL_DB = os.getenv("MYSQL_DB")

    # Tools directories
    SOURCE_PLINK_GENOME = "/data/clu/ukbb/genotypes_nomaly"
    # SOURCE_PLINK_GENOME = "/data/clu/ukbb/genotypes_nomaly_eur"
    PLINK_BINARY = "/data/clu/ukbb/plink"

    # Important mapping file
    NOMALY_VARIANTS_PATH = Path("/data/clu/ukbb/nomaly_variants.tsv")

    # Phenotype directories and files
    UKBB_PHENO_DIR = Path("/data/general/UKBB/Phenotypes")
    # PHENOTYPES_H5 = UKBB_PHENO_DIR / "phecode_cases_excludes.h5"
    PHENOTYPES_H5 = UKBB_PHENO_DIR / "phecode_cases_excludes_with_metadata.h5"
    # PHENOTYPES_PKL = UKBB_PHENO_DIR / "phecode_cases_excludes.pkl"

    # Genotype directories and files
    UKBB_GENOTYPES_DIR = Path("/data/general/UKBB/Genotypes/GRCh38")
    # GENOTYPES_H5 = UKBB_GENOTYPES_DIR / "genotypes_with_counts.h5"
    GENOTYPES_H5 = UKBB_GENOTYPES_DIR / "genotypes_with_metadata.h5"

    # Nomaly results directories and files (V1)
    NOMALY_RESULTS_DIR = Path("/data/general/UKBB/Run-v1/DatabaseInputs")
    NOMALY_SCORES_H5 = NOMALY_RESULTS_DIR / "float16_scores.h5"
    # STATS_H5 = NOMALY_RESULTS_DIR / "stats.h5"
    # STATS_H5 = NOMALY_RESULTS_DIR / "stats-fixed.h5"
    # STATS_H5 = NOMALY_RESULTS_DIR / "stats-EUR-2025-02-07.h5"
    STATS_H5 = NOMALY_RESULTS_DIR / "stats-All-2025-02-10.h5"

    # Nomaly results directories and files (V2)
    NOMALY_RESULTS_DIR_V2 = Path("/data/general/UKBB/Run-v2/DatabaseInputs")
    NOMALY_SCORES_H5_V2 = NOMALY_RESULTS_DIR_V2 / "float16_scores.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-fixed.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-fixed-EUR.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-EUR-2025-02-07.h5"
    STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-All-2025-02-10.h5"

    RESOURCE_DATA_DIR = Path("/data/general/Data/")

    # Caching directoriesd
    # GWAS_PHENO_DIR = Path("/data/clu/ukbb/by_pheno")
    # PHEWAS_PHENO_DIR = Path("/data/clu/ukbb/by_variant")
    # PHECODE_TERM_DIR = Path("/data/clu/ukbb/by_phecode_term")
    # VARIANT_SCORES_DIR = Path("/data/clu/ukbb/variant_scores")

    # Caching directories
    GWAS_PHENO_DIR = Path("/data/personal/danbolser/ukbb/cache/by_pheno")
    PHEWAS_PHENO_DIR = Path("/data/personal/danbolser/ukbb/cache/by_variant")
    PHECODE_TERM_DIR = Path("/data/personal/danbolser/ukbb/cache/by_phecode_term")
    VARIANT_SCORES_DIR = Path("/data/personal/danbolser/ukbb/cache/variant_scores")


class DevelopmentConfig(Config):
    """Development configuration"""

    TESTING = True
    DEBUG = True

    # Session/Cookie settings
    SESSION_COOKIE_NAME = "session-dev"
    SESSION_COOKIE_SECURE = False  # Allow HTTP for development
    SESSION_COOKIE_DOMAIN = None  # Restrict to same domain
    SECRET_KEY = os.getenv("DEV_SECRET_KEY", "development-secret-key")

    # Other development settings
    WTF_CSRF_ENABLED = True


class TestingConfig(Config):
    """Testing configuration"""

    TESTING = True
    DEBUG = False

    # Session/Cookie settings
    SESSION_COOKIE_NAME = "session-test"
    SESSION_COOKIE_SECURE = False  # Allow HTTP for testing
    SESSION_COOKIE_DOMAIN = None  # Restrict to same domain
    SECRET_KEY = os.getenv("TEST_SECRET_KEY", "testing-secret-key")

    # Test database settings
    MYSQL_HOST = "localhost"
    MYSQL_PORT = 3306
    MYSQL_USER = "testuser"
    MYSQL_PASSWORD = "testpass"
    MYSQL_DB = "testdb"

    # Other test settings
    WTF_CSRF_ENABLED = False
    LOGIN_DISABLED = False


class ProductionConfig(Config):
    """Production configuration"""

    TESTING = False
    DEBUG = False

    SESSION_COOKIE_NAME = "session"
    SESSION_COOKIE_DOMAIN = None  # Set to your domain in deployment
    SECRET_KEY = os.getenv("PROD_SECRET_KEY")  # Must be set in production

    # Security settings
    WTF_CSRF_ENABLED = True


class StagingConfig(Config):
    """Staging configuration"""

    SESSION_COOKIE_NAME = "session-staging"
    SESSION_COOKIE_DOMAIN = None  # Set to your staging domain in deployment
    SECRET_KEY = os.getenv("STAGING_SECRET_KEY")  # Must be set in staging


# Configuration dictionary
config = {
    "default": DevelopmentConfig,
    "development": DevelopmentConfig,
    "testing": TestingConfig,
    "production": ProductionConfig,
    "staging": StagingConfig,
}
