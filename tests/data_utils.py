from config import Config


def production_data_available():
    """Return True if critical production data files are present."""
    required_paths = [
        Config.PHENOTYPES_HDF,
        Config.GENOTYPES_HDF,
        Config.STATS_H5,
    ]
    return all(path.exists() for path in required_paths)
