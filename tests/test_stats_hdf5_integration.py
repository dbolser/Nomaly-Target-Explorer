import numpy as np
import pandas as pd
import pytest

# Known test cases - these values will need to be updated with actual production data

SAMPLE_COUNT_ALL = [467305, 486145]  # v1 / v2
SAMPLE_COUNT_EUR = 449423

KNOWN_PHECODES = ["290.11", "250.2"]
KNOWN_BUGGY_TERMS = ["GO:0016403", "GO:0036265", "GO:0048712"]
KNOWN_TERMS = ["GO:0030800"]

# Just some random columns we expect...
EXPECTED_COLUMNS = [
    "num_rp",
    "num_rn",
    "mwu_pvalue",
    "tti_pvalue",
    "metric1_pvalue",
]


@pytest.mark.parametrize("version", ["stats", "stats_v2"])
def test_production_file_exists(integration_app, version):
    """Verify the production phenotype file exists and is readable."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = getattr(services, version)
        assert stats_service is not None
        assert hasattr(stats_service, "_hdf")
        assert hasattr(stats_service._hdf, "data")


@pytest.mark.parametrize("version", ["stats", "stats_v2"])
def test_get_stats_by_phecode_1(integration_app, version):
    """Test querying a known phecode from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = getattr(services, version)
        for phecode in KNOWN_PHECODES:
            data = stats_service._hdf.get_stats_by_phecode(phecode, statstype=None)
            assert isinstance(data, pd.DataFrame)
            assert len(data) > 0

            assert all(col in data.columns for col in EXPECTED_COLUMNS)


@pytest.mark.skip(reason="NEED TO UPDATE.")
@pytest.mark.parametrize("version", ["stats", "stats_v2"])
def test_counts_for_phecode_008(integration_app, version):
    """Test the counts for a phecode with zero excludes."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = getattr(services, version)
        data = stats_service._hdf.get_stats_by_phecode("008", statstype=None)

        # NOTE: Hard coded to EUR for now...
        # assert np.all(data["num_rp"] + data["num_rn"] == SAMPLE_COUNT_EUR)

        # Value for ALL
        assert np.all(data["num_rp"] + data["num_rn"] == SAMPLE_COUNT_ALL[0]) or np.all(
            data["num_rp"] + data["num_rn"] == SAMPLE_COUNT_ALL[1]
        )


@pytest.mark.parametrize("version", ["stats", "stats_v2"])
def test_get_stats_by_term_phecode(integration_app, version):
    """Test querying a known phecode from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = getattr(services, version)

        for phecode in KNOWN_PHECODES:
            for term in KNOWN_TERMS:
                single_term_data = stats_service._hdf.get_stats_by_term_phecode(
                    term, phecode
                )
                assert isinstance(single_term_data, dict)
                assert len(single_term_data) > 0


@pytest.mark.parametrize("version", ["stats", "stats_v2"])
def test_buggy_terms_are_missing(integration_app, version):
    """Test that buggy terms are missing from the data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = getattr(services, version)

        for term in KNOWN_BUGGY_TERMS:
            with pytest.raises(IndexError):
                _ = stats_service._hdf.get_stats_by_term_phecode(term, "250.2")


# TODO: Make a 'sane range' consistent across versions / ancestries?
@pytest.mark.skip(reason="NEED TO CHECK STATS DATA!!!.")
def test_sanity_check_a_specific_stat(integration_app):
    """Test querying specific known values from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = services.stats._hdf

        # TODO: Integrate stats service 'handler' for version and populations here!
        stats = stats_service.get_stats_by_term_phecode("GO:0030800", "250.2")

        # SELECT
        #   num_rp,num_rn,mwu_pvalue,tti_pvalue,metric1_pvalue,roc_stats_mcc_pvalue,roc_stats_yjs_pvalue
        # FROM
        #   stats_test_denormalised_run_v1_eur
        # WHERE
        #   phecode = '250.2' AND term = 'GO:0030800';

        # +---------+----------+-------------+-------------+----------------+----------------------+----------------------+
        # | num_rp  | num_rn   | mwu_pvalue  | tti_pvalue  | metric1_pvalue | roc_stats_mcc_pvalue | roc_stats_yjs_pvalue |
        # +---------+----------+-------------+-------------+----------------+----------------------+----------------------+
        # | 37073.0 | 408655.0 | 1.48228e-32 | 0.000793399 | 0.368464       | 6.17771e-08          | 1.23201e-07          |
        # +---------+----------+-------------+-------------+----------------+----------------------+----------------------+

        # assert stats["num_rp"] == 37073
        # assert stats["num_rn"] == 408655
        # assert stats["mwu_pvalue"] == pytest.approx(1.48228e-32, rel=1e-5)
        # assert stats["tti_pvalue"] == pytest.approx(0.000793399, rel=1e-5)
        # assert stats["metric1_pvalue"] == pytest.approx(0.368464, rel=1e-5)
        # assert stats["roc_stats_mcc_pvalue"] == pytest.approx(6.17771e-08, rel=1e-5)
        # assert stats["roc_stats_yjs_pvalue"] == pytest.approx(1.23201e-07, rel=1e-5)

        # DEES ' ALL' not 'EUR'
        # assert stats["mwu_pvalue"] == pytest.approx(6.35231e-13, rel=1e-5)
        # assert stats["metric1_pvalue"] == pytest.approx(0.00446315, rel=1e-5)

        stats = stats_service.get_stats_by_term_phecode("GO:0030800", "290.11")

        # +--------+----------+------------+------------+----------------+----------------------+----------------------+
        # | num_rp | num_rn   | mwu_pvalue | tti_pvalue | metric1_pvalue | roc_stats_mcc_pvalue | roc_stats_yjs_pvalue |
        # +--------+----------+------------+------------+----------------+----------------------+----------------------+
        # | 4044.0 | 421145.0 | 0.334302   | 0.24151    | 0.339152       | 0.86028              | 0.86028              |
        # +--------+----------+------------+------------+----------------+----------------------+----------------------+

        assert stats["num_rp"] == 4044
        assert stats["num_rn"] == 421145
        assert stats["mwu_pvalue"] == pytest.approx(0.334302, rel=1e-5)
        assert stats["tti_pvalue"] == pytest.approx(0.24151, rel=1e-5)
        assert stats["metric1_pvalue"] == pytest.approx(0.339152, rel=1e-5)
        assert stats["roc_stats_mcc_pvalue"] == pytest.approx(0.86028, rel=1e-5)
        assert stats["roc_stats_yjs_pvalue"] == pytest.approx(0.86028, rel=1e-5)

        # DEES ' ALL' not 'EUR'
        # assert stats["mwu_pvalue"] == pytest.approx(0.0770663, rel=1e-5)
        # assert stats["metric1_pvalue"] == pytest.approx(0.254069, rel=1e-5)


# Below we just call the stats service using different values for statstype...
# I don't think these vesions are used anywhere.


@pytest.mark.parametrize("version", ["stats", "stats_v2"])
def test_get_stats_by_phecode_2(integration_app, version):
    """Test querying a known phecode from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = getattr(services, version)
        for phecode in KNOWN_PHECODES:
            data = stats_service._hdf.get_stats_by_phecode(
                phecode, statstype="metric1_pvalue"
            )
            assert isinstance(data, np.ndarray)
            assert len(data) > 0


@pytest.mark.parametrize("version", ["stats", "stats_v2"])
def test_get_stats_by_phecode_3(integration_app, version):
    """Test querying a known phecode from production data."""
    with integration_app.app_context():
        services = integration_app.extensions["nomaly_services"]
        stats_service = getattr(services, version)
        for phecode in KNOWN_PHECODES:
            data = stats_service._hdf.get_stats_by_phecode(
                phecode, statstype=EXPECTED_COLUMNS
            )
            assert isinstance(data, np.ndarray)
            assert len(data) > 0
