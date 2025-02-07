import numpy as np
import pandas as pd
import pytest

# Known test cases - these values will need to be updated with actual production data
KNOWN_PHECODES = ["290.11", "250.2"]
KNOWN_BUGGY_TERMS = ["GO:0016403", "GO:0036265", "GO:0048712"]
KNOWN_TERMS = ["GO:0030800"]


def test_production_file_exists(integration_app):
    """Verify the production phenotype file exists and is readable."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = services.stats
        assert stats_service is not None
        assert hasattr(stats_service, "_hdf")
        assert hasattr(stats_service._hdf, "data")


def test_known_phecode_query(integration_app):
    """Test querying a known phecode from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = services.stats
        for phecode in KNOWN_PHECODES:
            single_phecode_data = stats_service._hdf.get_stats_by_phecode(
                phecode, statstype=None
            )
            assert isinstance(single_phecode_data, pd.DataFrame)
            assert len(single_phecode_data) > 0


def test_known_phecode_query_specific_stat(integration_app):
    """Test querying a known phecode from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = services.stats
        for phecode in KNOWN_PHECODES:
            single_phecode_data = stats_service._hdf.get_stats_by_phecode(
                phecode, statstype="metric1_pvalue"
            )
            assert isinstance(single_phecode_data, np.ndarray)
            assert len(single_phecode_data) > 0


def test_known_phecode_query_specific_stats(integration_app):
    """Test querying a known phecode from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = services.stats
        for phecode in KNOWN_PHECODES:
            single_phecode_data = stats_service._hdf.get_stats_by_phecode(
                phecode, statstype=["metric1_pvalue", "mwu_pvalue", "mcc_pvalue"]
            )
            assert isinstance(single_phecode_data, np.ndarray)
            assert len(single_phecode_data) > 0


def test_get_stats_by_term_phecode(integration_app):
    """Test querying a known phecode from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = services.stats

        for phecode in KNOWN_PHECODES:
            for term in KNOWN_TERMS:
                single_term_data = stats_service._hdf.get_stats_by_term_phecode(
                    term, phecode
                )
                assert isinstance(single_term_data, dict)
                assert len(single_term_data) > 0


def test_sanity_check_a_specific_stat(integration_app):
    """Test querying specific known values from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = services.stats._hdf

        # Test GO:0030800 - 250.2
        stats = stats_service.get_stats_by_term_phecode("GO:0030800", "250.2")
        assert stats["mwu_pvalue"] == pytest.approx(6.35231e-13, rel=1e-5)
        assert stats["metric1_pvalue"] == pytest.approx(0.00446315, rel=1e-5)

        # Test GO:0036265 - 250.2
        # stats = stats_service.get_stats_by_term_phecode("GO:0036265", "250.2")
        # assert stats["mwu_pvalue"] == pytest.approx(0.0587734, rel=1e-5)
        # assert stats["metric1_pvalue"] == pytest.approx(0.116089, rel=1e-5)

        # Test GO:0030800 - 290.11
        stats = stats_service.get_stats_by_term_phecode("GO:0030800", "290.11")
        assert stats["mwu_pvalue"] == pytest.approx(0.0770663, rel=1e-5)
        assert stats["metric1_pvalue"] == pytest.approx(0.254069, rel=1e-5)

        # Test GO:0036265 - 290.11
        # stats = stats_service.get_stats_by_term_phecode("GO:0036265", "290.11")
        # assert stats["mwu_pvalue"] == pytest.approx(0.817572, rel=1e-5)
        # assert stats["metric1_pvalue"] is None
