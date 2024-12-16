import os
from typing import Dict, List
import logging
import json
from datetime import datetime

logger = logging.getLogger(__name__)

CACHE_DIR = "/data/clu/ukbb/by_phecode_term/"


def ensure_cache_dir():
    """Ensure cache directory exists"""
    os.makedirs(CACHE_DIR, exist_ok=True)


def get_cache_path(phecode: str, term: str) -> str:
    """Get path for cached results"""
    return os.path.join(CACHE_DIR, f"phecode_{phecode}_term_{term}.json")


def delete_cache(phecode: str, term: str) -> bool:
    """Delete cached results if they exist"""
    cache_path = get_cache_path(phecode, term)
    if os.path.exists(cache_path):
        try:
            os.remove(cache_path)
            logger.info(f"Deleted cache for phecode {phecode}, term {term}")
            return True
        except Exception as e:
            logger.error(f"Error deleting cache: {e}")
    return False


def load_cached_results(phecode: str, term: str, flush: bool = False) -> Dict | None:
    """Load cached results if they exist and flush is not requested"""
    if flush:
        logger.info(f"Flush requested for phecode {phecode}, term {term}")
        delete_cache(phecode, term)
        return None

    cache_path = get_cache_path(phecode, term)
    logger.info(f"Checking cache at: {cache_path}")
    
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "r") as f:
                cached_data = json.load(f)
            logger.info(f"Successfully loaded cache for phecode {phecode}, term {term}")
            return cached_data
        except Exception as e:
            logger.error(f"Error loading cached results: {e}")
            # If there's an error reading the cache, delete it
            delete_cache(phecode, term)
            return None
    
    logger.info(f"No cache file found at {cache_path}")
    return None


def save_results(phecode: str, term: str, data: List[Dict]):
    """Save results to cache"""
    ensure_cache_dir()
    cache_path = get_cache_path(phecode, term)
    try:
        with open(cache_path, "w") as f:
            json.dump(
                {
                    "timestamp": datetime.now().isoformat(),
                    "phecode": phecode,
                    "term": term,
                    "data": data,
                },
                f,
                indent=2,
            )
        logger.info(f"Saved results to cache for phecode {phecode}, term {term}")
    except Exception as e:
        logger.error(f"Error saving results to cache: {e}")
