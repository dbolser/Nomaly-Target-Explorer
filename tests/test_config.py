from pathlib import Path
from config import (
    Config,
    DevelopmentConfig,
    TestingConfig,
    ProductionConfig,
    StagingConfig,
    config as config_dict,
)

import pytest

from tests.data_utils import production_data_available

pytestmark = pytest.mark.skipif(
    not production_data_available(),
    reason="requires production configuration files",
)


def test_base_paths_exist():
    """Test that all base paths that don't change with population exist."""

    config = Config()

    # Test resource directories
    assert config.NOMALY_VARIANT_MAPPING_PATH.exists(), (
        "NOMALY_VARIANT_MAPPING_PATH does not exist"
    )
    assert config.NOMALY_VARIANT_SCORES_PATH.exists(), (
        "NOMALY_VARIANT_SCORES_PATH does not exist"
    )
    assert config.NOMALY_VARIANT2GENE_PATH.exists(), (
        "NOMALY_VARIANT2GENE_PATH does not exist"
    )
    assert config.NOMALY_RESULTS_DIR.exists(), "NOMALY_RESULTS_DIR does not exist"
    assert config.UKBB_PHENO_DIR.exists(), "UKBB_PHENO_DIR does not exist"
    assert config.UKBB_GENOTYPES_DIR.exists(), "UKBB_GENOTYPES_DIR does not exist"
    assert config.PHENOTYPES_HDF.exists(), "PHENOTYPES_HDF does not exist"
    assert config.GENOTYPES_BED.exists(), "GENOTYPES_BED does not exist"
    assert config.GENOTYPES_BIM.exists(), "GENOTYPES_BIM does not exist"
    assert config.GENOTYPES_FAM.exists(), "GENOTYPES_FAM does not exist"
    assert config.GENOTYPES_HDF.exists(), "GENOTYPES_HDF does not exist"
    assert config.GENOTYPES_NPY.exists(), "GENOTYPES_NPY does not exist"


def test_pharros_data_dir_exists():
    """Test that the Pharos data directory exists."""
    config = Config()
    assert config.PHAROS_DATA_DIR.exists(), "PHAROS_DATA_DIR does not exist"


def test_cache_dirs_exist():
    """Test that all cache directories exist."""
    config = Config()
    assert config.GWAS_PHENO_DIR.exists(), "GWAS_PHENO_DIR does not exist"
    assert config.PHEWAS_PHENO_DIR.exists(), "PHEWAS_PHENO_DIR does not exist"
    assert config.PHECODE_TERM_DIR.exists(), "PHECODE_TERM_DIR does not exist"
    assert config.VARIANT_SCORES_DIR.exists(), "VARIANT_SCORES_DIR does not exist"


@pytest.mark.parametrize(
    "env, expected_class",
    [
        ("default", DevelopmentConfig),
        ("development", DevelopmentConfig),
        ("testing", TestingConfig),
        ("production", ProductionConfig),
        ("staging", StagingConfig),
    ],
)
def test_config_dictionary(env, expected_class):
    """Test that the config dictionary maps environments to the correct classes."""
    assert env in config_dict
    assert config_dict[env] == expected_class
