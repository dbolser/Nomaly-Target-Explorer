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
    PLINK_BINARY = "/data/clu/ukbb/plink"

    # Important mapping files and other junk
    NOMALY_VARIANT_MAPPING_PATH = Path("/data/clu/ukbb/nomaly_variants.tsv")
    NOMALY_VARIANT_SCORES_PATH = Path("/data/clu/ukbb/variantscores.tsv")
    NOMALY_VARIANT_SCORES_PATH = Path("/data/clu/ukbb/genotype_counts_with_vs.tsv")
    NOMALY_VARIANT2GENE_PATH = Path("/data/clu/ukbb/variant2gene.tsv")

    DATA_ROOT = Path("/data/analysis/UKBB")

    # Phenotype directories and files
    UKBB_PHENO_DIR = DATA_ROOT / "Phenotypes"
    # PHENOTYPES_H5 = UKBB_PHENO_DIR / "phecode_cases_excludes.h5"
    PHENOTYPES_HDF = UKBB_PHENO_DIR / "phecode_cases_excludes_with_metadata.h5"
    # PHENOTYPES_PKL = UKBB_PHENO_DIR / "phecode_cases_excludes.pkl"

    # Genotype directories and files
    UKBB_GENOTYPES_DIR = DATA_ROOT / "Genotypes/GRCh38"

    GENOTYPES_BED = UKBB_GENOTYPES_DIR / "genotypes-ukbb.bed"
    GENOTYPES_BIM = UKBB_GENOTYPES_DIR / "genotypes-ukbb.bim"
    GENOTYPES_FAM = UKBB_GENOTYPES_DIR / "genotypes-ukbb.fam"

    GENOTYPES_HDF = UKBB_GENOTYPES_DIR / "genotypes-ukbb.h5"
    GENOTYPES_NPY = UKBB_GENOTYPES_DIR / "genotypes-ukbb.npy"

    # Nomaly results directories and files (V1)
    NOMALY_RESULTS_DIR = DATA_ROOT / "Run-v1/DatabaseInputs"
    NOMALY_SCORES_H5 = NOMALY_RESULTS_DIR / "float32_scores.h5"
    # STATS_H5 = NOMALY_RESULTS_DIR / "stats.h5"
    # STATS_H5 = NOMALY_RESULTS_DIR / "stats-fixed.h5"
    # STATS_H5 = NOMALY_RESULTS_DIR / "stats-EUR-2025-02-07.h5"
    STATS_H5 = NOMALY_RESULTS_DIR / "stats-All-2025-02-10.h5"

    # TODO: Need to implement this in VP code
    # Nomaly results directories and files (V2)
    NOMALY_RESULTS_DIR_V2 = DATA_ROOT / "Run-v2/DatabaseInputs"
    NOMALY_SCORES_H5_V2 = NOMALY_RESULTS_DIR_V2 / "float32_scores.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-fixed.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-fixed-EUR.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-EUR-2025-02-07.h5"
    STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-All-2025-02-10.h5"

    STATS_SELECTOR = {
        "Run-v1": {
            "AFR": DATA_ROOT / "Run-v1/DatabaseInputs/AFR-sexmatch.h5",
            "EUR": DATA_ROOT / "Run-v1/DatabaseInputs/EUR-sexmatch.h5",
            "SAS": DATA_ROOT / "Run-v1/DatabaseInputs/SAS-sexmatch.h5",
            "EAS": DATA_ROOT / "Run-v1/DatabaseInputs/EAS-sexmatch.h5",
            "ALL": DATA_ROOT / "Run-v1/DatabaseInputs/ALL-sexmatch.h5",
        },
        "Run-v2": {
            "AFR": DATA_ROOT / "Run-v2/DatabaseInputs/AFR-sexmatch.h5",
            "EUR": DATA_ROOT / "Run-v2/DatabaseInputs/EUR-sexmatch.h5",
            "SAS": DATA_ROOT / "Run-v2/DatabaseInputs/SAS-sexmatch.h5",
            "EAS": DATA_ROOT / "Run-v2/DatabaseInputs/EAS-sexmatch.h5",
            "ALL": DATA_ROOT / "Run-v2/DatabaseInputs/ALL-sexmatch.h5",
        },
    }

    PHAROS_DATA_DIR = Path("/data/public/Pharos")

    # Caching directories
    GWAS_PHENO_DIR = Path("/data/danbolser/ukbb/cache/by_pheno")
    PHEWAS_PHENO_DIR = Path("/data/danbolser/ukbb/cache/by_variant")
    PHECODE_TERM_DIR = Path("/data/danbolser/ukbb/cache/by_phecode_term")
    VARIANT_SCORES_DIR = Path("/data/danbolser/ukbb/cache/variant_scores")


class DevelopmentConfig(Config):
    """Development configuration"""

    TESTING = False
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

    # Mock file paths for testing - these should be overridden in tests
    # GENOTYPES_HDF = "MOCK_PATH_OVERRIDE_IN_TESTS"
    # PHENOTYPES_H5 = "MOCK_PATH_OVERRIDE_IN_TESTS"
    # NOMALY_SCORES_H5 = "MOCK_PATH_OVERRIDE_IN_TESTS"
    # STATS_H5 = "MOCK_PATH_OVERRIDE_IN_TESTS"

    # Other test settings
    WTF_CSRF_ENABLED = False
    LOGIN_DISABLED = False


class ProductionConfig(Config):
    """Production configuration"""

    TESTING = False
    DEBUG = False

    SESSION_COOKIE_NAME = "session-live"
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
