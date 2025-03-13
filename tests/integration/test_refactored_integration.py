import json
import pytest

from blueprints.prioritisation_by_nomaly_scores import (
    get_top_variants,
    convert_numpy_types,
    check_json_safety,
)


def get_expected_output(phecode: str, term: str):
    # These are just copied over from the cache directory!
    from config import Config

    CACHE_DIR = Config.VARIANT_SCORES_DIR

    with open(f"{CACHE_DIR}/variant_prioritization_{phecode}_{term}.json", "r") as f:
        return json.load(f)


@pytest.mark.integration
@pytest.mark.parametrize(
    "phecode_term",
    [
        ("250.2", "UP:UPA00246"),
        ("282", "GO:0004586"),
        ("705.3", "CC:MESH:C028632"),
        ("256", "GO:0003941"),
    ],
)
def test_refactored_json_output(integration_app, phecode_term):
    """Integration test for the refactored get_top_variants function.

    This test calls get_top_variants with real (or near-real) service
    implementations and verifies that the final JSON output contains the
    expected keys.

    This serves as a snapshot test to ensure refactoring preserves the original
    output behavior.
    """
    phecode, term = phecode_term

    expected_output = get_expected_output(phecode, term)

    services = integration_app.extensions["nomaly_services"]

    phenotype_service = services.phenotype
    genotype_service = services.genotype
    nomaly_scores_service = services.nomaly_score
    stats_service = services.stats

    # Force recomputation by setting no_cache=True
    result = get_top_variants(
        phecode,
        term,
        phenotype_service,
        genotype_service,
        nomaly_scores_service,
        stats_service,
        no_cache=True,
    )

    result = convert_numpy_types(result)

    # Define essential keys that should be present in the output JSON
    assert result == expected_output

    # Optionally, perform further structure/value tests or compare with a snapshot.
