import json
import pytest
from unittest.mock import patch

from blueprints.prioritisation_by_nomaly_scores import (
    get_top_variants as get_top_variants_original,
    convert_numpy_types,
)

from blueprints.prioritisation_by_nomaly_scores_refactored import (
    get_top_variants_refactored as get_top_variants_refactored,
)


def get_expected_output(phecode: str, term: str, run_version: str, ancestry: str):
    from config import Config

    CACHE_DIR = Config.VARIANT_SCORES_DIR
    cache_file = (
        f"variant_prioritization_{phecode}_{term}_{run_version}_{ancestry}.json"
    )

    with open(CACHE_DIR / cache_file, "r") as f:
        return json.load(f)


# @pytest.mark.integration
@pytest.mark.parametrize(
    "phecode_term",
    [
        ("250.2", "UP:UPA00246"),
        ("282", "GO:0004586"),
        ("705.3", "CC:MESH:C028632"),
        ("256", "GO:0003941"),
    ],
)
def test_refactored_matches_original(integration_app, phecode_term):
    """Test that the refactored get_top_variants function produces the same output as the original."""
    phecode, term = phecode_term

    run_version = "Run-v1"
    ancestry = "EUR"

    services = integration_app.extensions["nomaly_services"]

    phenotype_service = services.phenotype
    genotype_service = services.genotype
    nomaly_scores_service = services.nomaly_score
    stats_service = services.stats_registry.get(run_version, ancestry)

    # Force recomputation by setting no_cache=True
    original_result = get_top_variants_original(
        phecode,
        term,
        phenotype_service,
        genotype_service,
        nomaly_scores_service,
        stats_service,
        run_version=run_version,
        ancestry=ancestry,
        no_cache=True,
    )

    refactored_result = get_top_variants_refactored(
        phecode,
        term,
        phenotype_service,
        genotype_service,
        nomaly_scores_service,
        stats_service,
        run_version=run_version,
        ancestry=ancestry,
        no_cache=True,
    )

    # Convert to standard Python types for comparison
    original_result = convert_numpy_types(original_result)
    refactored_result = convert_numpy_types(refactored_result)

    # Use deepdiff to find differences between the two results
    # This will provide a detailed report of any differences
    diff = {}
    try:
        from deepdiff import DeepDiff
        import pprint

        diff = DeepDiff(
            original_result,
            refactored_result,
            ignore_order=True,  # Ignore list order differences
            significant_digits=6,  # Allow small floating point differences
        )
    except ImportError:
        # Fallback if deepdiff is not available
        print("DeepDiff not available, using simple comparison")
        assert original_result == refactored_result, "Results don't match"
        return

    # If there are differences, print them for debugging
    if diff:
        print("\nDifferences between original and refactored results:")
        pprint.pprint(diff)

        # Check if differences are only in floating point values with small differences
        # This helps identify if the differences are just due to floating point precision
        only_numeric_diffs = True
        for key in diff:
            if key == "values_changed":
                for change_key, change_data in diff[key].items():
                    if not (
                        isinstance(change_data.get("old_value"), (int, float))
                        and isinstance(change_data.get("new_value"), (int, float))
                    ):
                        only_numeric_diffs = False
                        break

                    if abs(change_data["old_value"] - change_data["new_value"]) > 1e-6:
                        only_numeric_diffs = False
                        break
            else:
                only_numeric_diffs = False
                break

        if only_numeric_diffs:
            print("Only small floating point differences detected - this is acceptable")
        else:
            # If there are significant differences, fail the test
            assert not diff, (
                "Refactored function produces different output than original"
            )


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
    run_version = "Run-v1"
    ancestry = "EUR"

    expected_output = get_expected_output(phecode, term, run_version, ancestry)

    services = integration_app.extensions["nomaly_services"]

    phenotype_service = services.phenotype
    genotype_service = services.genotype
    nomaly_scores_service = services.nomaly_score
    stats_service = services.stats_registry.get(run_version, ancestry)

    # NOTE: We need to mock the call to 'save_results_to_cache' to keep the test
    # accurate (we don't want it to overwrite the original cache)
    with patch(
        "blueprints.prioritisation_by_nomaly_scores_refactored.save_results_to_cache"
    ) as mock_save_results_to_cache:
        mock_save_results_to_cache.return_value = expected_output

        # Force recomputation by setting no_cache=True
        result = get_top_variants_refactored(
            phecode,
            term,
            phenotype_service,
            genotype_service,
            nomaly_scores_service,
            stats_service,
            run_version=run_version,
            ancestry=ancestry,
            no_cache=True,
        )

    result = convert_numpy_types(result)

    # Use deepdiff for detailed comparison
    from deepdiff import DeepDiff

    diff = DeepDiff(
        expected_output,
        result,
        ignore_order=True,
        significant_digits=6,
    )

    assert not diff, "Results don't match expected output"
