import json
import logging
from pathlib import Path
from typing import Dict, List

from config import Config


class CacheMiss(Exception):
    """Exception raised when cache miss occurs."""

    pass


logger = logging.getLogger(__name__)

PHECODE_TERM_DIR: Path = Config.PHECODE_TERM_DIR


def get_cache_path(phecode: str, term: str, ancestry: str) -> Path:
    """Get cache path."""
    return PHECODE_TERM_DIR / f"ancestry_{ancestry}_phecode_{phecode}_term_{term}.json"


def load_cache_results(phecode: str, term: str, ancestry: str) -> Dict:
    """Load results from cache."""

    cache_path = get_cache_path(phecode, term, ancestry)
    logger.debug(f"Checking cache at: {cache_path}")

    if cache_path.exists():
        try:
            with open(cache_path, "r") as f:
                cached_data = json.load(f)
            return cached_data
        except Exception as e:
            raise Exception(f"Error loading cached results: {e}")

    logger.info(f"No cache found at {cache_path}")
    raise CacheMiss(f"No cache found at {cache_path}")


def save_cache_results(phecode: str, term: str, ancestry: str, data: List[Dict]):
    """Save results to cache."""

    cache_path = get_cache_path(phecode, term, ancestry)
    try:
        with open(cache_path, "w") as f:
            json.dump(
                {
                    "phecode": phecode,
                    "term": term,
                    "ancestry": ancestry,
                    "data": data,
                },
                f,
                indent=2,
            )
        logger.info(
            f"Saved results to cache for phecode {phecode}, term {term}, ancestry {ancestry}"
        )
    except Exception as e:
        logger.error(f"Error saving results to cache: {e}")
