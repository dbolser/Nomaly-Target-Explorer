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
    SECRET_KEY = os.environ.get("SECRET_KEY", "your-default-secret-key")

    # Database settings
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
    STATS_H5 = NOMALY_RESULTS_DIR / "stats-EUR-2025-02-07.h5"

    # Nomaly results directories and files (V2)
    NOMALY_RESULTS_DIR_V2 = Path("/data/general/UKBB/Run-v2/DatabaseInputs")
    NOMALY_SCORES_H5_V2 = NOMALY_RESULTS_DIR_V2 / "float16_scores.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-fixed.h5"
    # STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-fixed-EUR.h5"
    STATS_H5_V2 = NOMALY_RESULTS_DIR_V2 / "stats-EUR-2025-02-07.h5"

    RESOURCE_DATA_DIR = Path("/data/general/Data/")

    # Caching directoriesd
    GWAS_PHENO_DIR = Path("/data/clu/ukbb/by_pheno")
    PHEWAS_PHENO_DIR = Path("/data/clu/ukbb/by_variant")
    PHECODE_TERM_DIR = Path("/data/clu/ukbb/by_phecode_term")
    VARIANT_SCORES_DIR = Path("/data/clu/ukbb/variant_scores")

    # Caching directories
    # GWAS_PHENO_DIR = Path("/data/personal/danbolser/ukbb/eur_gwas")
    # PHEWAS_PHENO_DIR = Path("/data/personal/danbolser/ukbb/eur_phewas")
    # PHECODE_TERM_DIR = Path("/data/personal/danbolser/ukbb/by_phecode_term")


class DevelopmentConfig(Config):
    """Development configuration"""

    DEBUG = True


class ProductionConfig(Config):
    """Production configuration"""

    WTF_CSRF_ENABLED = True


config = {
    "default": DevelopmentConfig,
    "development": DevelopmentConfig,
    "production": ProductionConfig,
}
