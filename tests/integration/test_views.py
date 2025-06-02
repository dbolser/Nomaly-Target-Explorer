import json
import time
import pytest

def test_search_route_unauthenticated(integration_app_client):
    """Test a random route redirects to index when not authenticated."""
    response = integration_app_client.get("/search")
    assert response.status_code == 302  # Redirect to login
    assert response.location == "/login?next=http://localhost.localdomain/search"


def test_index_route_unauthenticated(integration_app_client):
    """Test the index route when authenticated."""
    response = integration_app_client.get("/")
    assert response.status_code == 200
    assert b"welcome" in response.data.lower()


def test_login_route(integration_app_client, test_admin):
    """Test login functionality."""
    # Test GET request
    response = integration_app_client.get("/login")
    assert response.status_code == 200
    assert b"login" in response.data.lower()

    # Test POST request with valid credentials
    response = integration_app_client.post(
        "/login",
        data={"username": test_admin["username"], "password": test_admin["password"]},
    )
    assert response.status_code == 302  # Should be a redirect
    assert response.location == "/"  # Should redirect to index

    response = integration_app_client.get("/logout")
    assert response.status_code == 302  # Should be a redirect
    assert response.location == "/"  # Should redirect to index

    # Test POST request with invalid credentials
    response = integration_app_client.post(
        "/login",
        data={"username": test_admin["username"], "password": "wr0nG_pa55woRd"},
    )
    assert response.status_code == 200  # Should stay on login page
    assert b"invalid username or password" in response.data.lower()


def test_search_route(auth_integration_app_client):
    """Test the disease search functionality."""
    response = auth_integration_app_client.get("/diseasesearch?q=diabetes")
    assert response.status_code == 200
    data = json.loads(response.data)
    assert isinstance(data, list)

    # Verify the structure of returned data
    if len(data) > 0:
        assert all(
            key in data[0]
            for key in [
                "phecode",
                "description",
                "sex",
                "affected",
                "excluded",
                "phecode_exclude",
                "phecode_group",
            ]
        )


def test_phecode_route(auth_integration_app_client):
    """Test the phecode detail route."""
    test_phecode = "250.2"  # Example phecode
    response = auth_integration_app_client.get(f"/phecode/{test_phecode}")
    assert response.status_code == 200
    assert bytes(test_phecode, "utf-8") in response.data


def test_disease_sets1_structure(auth_integration_app_client):
    """Test the structure and content of the page1 route."""
    response = auth_integration_app_client.get("/disease-sets/set1")
    assert response.status_code == 200

    # Convert response data to string for easier testing
    html = response.data.decode("utf-8")

    # Test main heading
    assert '<h1 class="text-center">Selected diseases</h1>' in html

    # Test category headings (H2s)
    assert '<h2 class="mt-5 mb-4 ps-2">Skin Conditions</h2>' in html

    # Test result sections (H3s)
    assert 'Results for "Hidradenitis"' in html

    # Test description text under each search
    assert "Any Phecode or ICD10 descriptions that contains the word" in html

    # Test dynamic content loading
    assert "results-list-1" in html  # First results list
    assert "results-list-2" in html  # Second results list

    # Test JavaScript initialization
    assert "searchData(" in html
    assert "async function searchData(query, listIndex)" in html


def test_disease_sets2_structure(auth_integration_app_client):
    """Test the structure and content of the page1 route."""
    response = auth_integration_app_client.get("/disease-sets/set2")
    assert response.status_code == 200

    # Convert response data to string for easier testing
    html = response.data.decode("utf-8")

    # Test main heading
    assert '<h1 class="text-center">Selected diseases</h1>' in html

    # Test category headings (H2s)
    assert '<h2 class="mt-5 mb-4 ps-2">Women&#39;s Health</h2>' in html

    # Test result sections (H3s)
    assert 'Results for "Polycystic Ovarian"' in html

    # Test description text under each search
    assert "Any Phecode or ICD10 descriptions that contains the word" in html

    # Test dynamic content loading
    assert "results-list-1" in html  # First results list
    assert "results-list-2" in html  # Second results list

    # Test JavaScript initialization
    assert "searchData(" in html
    assert "async function searchData(query, listIndex)" in html


def test_phecode_term_structure(auth_integration_app_client):
    """Test the structure and content of a specific phecode term page."""
    # phecode = "649.1"
    phecode = "635.2"
    term = "GO:0035235"

    response = auth_integration_app_client.get(f"/phecode/{phecode}/term/{term}")
    assert response.status_code == 200

    html = response.data.decode("utf-8")

    # Test phecode section
    assert "<h3>Disease: " in html
    assert f"/phecode/{phecode}" in html
    assert f"PheCode {phecode}" in html

    assert "Sex: Female" in html
    assert "<span>Affected: <strong>" in html
    assert "Excluded: " in html
    assert "Disease: pregnancy complications" in html

    # Test term section
    assert "ionotropic glutamate receptor signaling pathway" in html
    assert f"Term: {term}" in html
    assert "23 domains" in html
    assert "27 genes" in html
    assert "GRID1" in html
    assert "GRIK4" in html

    # Test table container exists
    assert '<div id="tableContainer"' in html
    assert '<table id="resultsTable"' in html

    # Test JavaScript initialization
    assert f'const phecode = "{phecode}";' in html
    assert f'const term = "{term}";' in html


def test_phecode_page_structure(auth_integration_app_client):
    """Test the structure and content of a specific phecode page."""
    phecode = "649.1"

    response = auth_integration_app_client.get(f"/phecode/{phecode}")
    assert response.status_code == 200

    html = response.data.decode("utf-8")

    # Test main content
    assert "Diabetes or abnormal glucose tolerance complicating pregnancy" in html
    assert "<span>Phecode: " in html
    assert "Sex: Female" in html
    assert "Population: EUR" in html
    assert "Affected: <strong>" in html
    assert "Control: " in html
    assert "Excluded: " in html
    assert "Disease: pregnancy complications" in html

    # Test Nomaly Results section
    assert "Per Term Nomaly Results" in html
    assert '<table id="resultsTable"' in html

    # Test JavaScript initialization
    assert f'const phecode = "{phecode}";' in html
    # assert 'const runbatch = "Run v1";' in html


def test_phecode_page_with_gwas(auth_integration_app_client):
    """Test the phecode page with GWAS functionality enabled."""
    phecode = "649.1"

    # response = auth_integration_app_client.get(f"/phecode/{phecode}?gwas=1")
    # It's now on all the time
    response = auth_integration_app_client.get(f"/phecode/{phecode}")
    assert response.status_code == 200

    html = response.data.decode("utf-8")

    # Test GWAS button exists
    assert (
        '<button id="runTaskButton" class="btn btn-primary">Run GWAS</button>' in html
    )

    # Test GWAS results containers exist
    assert '<div id="taskResult"' in html
    assert '<div id="GWASTableContainer"' in html
    assert '<table id="variantTable"' in html

    # Test GWAS table headers
    assert "<th>Variant</th>" in html
    assert "<th>Gene</th>" in html
    assert "<th>RSID</th>" in html
    assert "<th>F_A</th>" in html
    assert "<th>F_U</th>" in html
    assert "<th>OR</th>" in html
    assert "<th>P</th>" in html


def test_phecode_nomaly_stats(auth_integration_app_client):
    """Test the Nomaly stats endpoint."""
    phecode = "649.1"

    response = auth_integration_app_client.post(f"/nomaly-stats/{phecode}")
    assert response.status_code == 200

    data = json.loads(response.data)

    # Check response structure
    # assert "qqplot" in data
    assert "data" in data
    assert "columns" in data
    assert "columnNames" in data
    assert "defaultColumns" in data
    # assert "numColumns" in data

    # Check expected columns exist
    expected_columns = [
        "minrank",
        "term",
        "name",
        # "domain",
        "mwu_pvalue",
        "mcc_pvalue",
        "yjs_pvalue",
        "lrp_pvalue",
        "metric1_pvalue",
        "lrn_protective_pvalue",
    ]
    assert all(col in data["columns"] for col in expected_columns)

    # Check data records if any exist
    if len(data["data"]) > 0:
        first_record = data["data"][0]
        # Check required fields in first record
        assert "term" in first_record
        assert "name" in first_record
        # assert "domain" in first_record
        assert "minrank" in first_record

        # Check p-value formatting
        for pval_field in ["mwu_pvalue", "metric1_pvalue", "mcc_pvalue"]:
            if pval_field in first_record:
                pval = first_record[pval_field]
                if pval != "nan":
                    # Should be in scientific notation
                    assert "e" in pval.lower()

        # Check data limits
        assert len(data["data"]) <= 1050  # Should be limited to 1000 entries


def test_variant_page_structure(auth_integration_app_client):
    """Test the structure and content of a specific variant page."""
    variant = "17_80117714_G_A"

    response = auth_integration_app_client.get(f"/variant/{variant}")
    assert response.status_code == 200

    html = response.data.decode("utf-8")

    # Test page title and headers
    # assert "<h3>Variant Details</h3>" in html
    assert '<h5 class="card-title">Basic Information' in html

    # Test variant information
    assert "<strong>Variant ID:</strong>" in html
    assert '<span id="variant-display">17_80117714_G/A</span>' in html
    assert "<strong>Chromosome:</strong> 17" in html
    assert "<strong>Position:</strong> 80117714" in html

    # NOTE THAT IT'S NOT FLIPPED!?!
    assert "<strong>Genotyping Ref Allele:</strong> G" in html
    assert "<strong>Genotyping Alt Allele:</strong> A" in html

    # Test PheWAS button exists
    assert (
        '<button id="runPheWASButton" class="btn btn-primary">PheWAS Task</button>'
        in html
    )

    # Test containers for results exist
    assert '<div id="phewasResult"' in html
    assert '<div id="phewasTableContainer"' in html

    # We need to click the button to get the results!

    # assert "PheWAS identified 327 phecodes with association p<0.05" in html


def test_variant_phewas_results(auth_integration_app_client):
    """Test the PheWAS results endpoint for a variant."""
    variant = "17_80117714_G_A"

    # First check initial state
    response = auth_integration_app_client.get(f"/phewas-result/{variant}")
    assert response.status_code == 200

    initial_data = json.loads(response.data)
    assert initial_data["result"] == "Processing..."
    assert initial_data["associations"] == []

    # Trigger PheWAS analysis
    response = auth_integration_app_client.post(f"/run-phewas/{variant}?flush=1")
    assert response.status_code == 202
    assert json.loads(response.data)["status"] == "Task started"

    # Poll for results with timeout
    max_wait = 60  # Maximum wait time in seconds
    start_time = time.time()
    data = None

    while time.time() - start_time < max_wait:
        response = auth_integration_app_client.get(f"/phewas-result/{variant}")
        assert response.status_code == 200
        data = json.loads(response.data)

        # If we have results or a failure, break
        if data["result"] != "Processing...":
            break

        time.sleep(1)  # Wait 1 second before next poll

    # Ensure we didn't timeout
    assert time.time() - start_time < max_wait, "Timed out waiting for PheWAS results"
    assert data is not None
    assert data["result"] != "Processing..."

    if "Failed" not in data["result"]:
        # Check associations data structure
        if len(data["associations"]) > 0:
            first_assoc = data["associations"][0]
            expected_fields = [
                "Phecode",
                "Sex",
                "Description",
                "Group",
                "P",
                "OR",
                "Counts",
                "RefAF",
                "AltAF",
            ]
            assert all(field in first_assoc for field in expected_fields)


def test_variant_id_formats(auth_integration_app_client):
    """Test different variant ID format handling."""
    variants = [
        "17_80117714_G_A",  # Underscore format
    ]

    for variant in variants:
        response = auth_integration_app_client.get(f"/variant/{variant}")
        assert response.status_code == 200
        html = response.data.decode("utf-8")
        assert "<strong>Chromosome:</strong> 17" in html
        assert "<strong>Position:</strong> 80117714" in html


def test_phecode_gwas_pvalues(auth_integration_app_client):
    """Test that GWAS P-values are present in phecode page for a specific case."""
    phecode = "561"

    # First check the page with GWAS enabled
    response = auth_integration_app_client.get(f"/phecode/{phecode}")
    assert response.status_code == 200
    html = response.data.decode("utf-8")
    assert (
        '<button id="runTaskButton" class="btn btn-primary">Run GWAS</button>' in html
    )

    # Now simulate clicking the GWAS button by calling the run-task endpoint
    response = auth_integration_app_client.post(f"/run-task/{phecode}/1")
    # response = auth_integration_app_client.post(f"/run-task/{phecode}/0")
    assert response.status_code == 200
    data = json.loads(response.data)

    # Verify GWAS completed successfully
    assert data["status"] == "completed"
    assert "associations" in data
    assert len(data["associations"]) == 83011

    # Check P-values in GWAS results
    for assoc in data["associations"]:
        assert "P" in assoc
        assert assoc["P"] is not None
        p_value = float(assoc["P"])
        assert p_value >= 0
        assert p_value <= 1  # P-values should be between 0 and 1

    # The first record should have the lowest P-value
    first_assoc = data["associations"][0]
    assert first_assoc["Variant"] == "4:73576632_C/G"
    assert first_assoc["P"] == 0.000002986
    assert first_assoc["OR"] == 0.9197
    assert first_assoc["F_U"] == 0.05884
    assert first_assoc["F_A"] == 0.05437


def test_phecode_term_gwas_pvalues(auth_integration_app_client):
    """Test that GWAS P-values are present in phecode term page for a specific case."""
    phecode = "561"
    term = "HP:0002242"

    # Get the variant detail data
    response = auth_integration_app_client.get(
        f"/phecode/{phecode}/term/{term}/tableVariantDetail?flush=0"
    )
    assert response.status_code == 200
    data = json.loads(response.data)

    # Verify we have data
    assert "data" in data
    assert len(data["data"]) > 0

    # Check that GWAS P-values are present and not empty in at least one record
    has_pvalue = False
    for record in data["data"]:
        if record.get("GWAS_P") and record["GWAS_P"] != "":
            has_pvalue = True
            # Convert to float to ensure it's a valid number
            p_value = float(record["GWAS_P"])
            assert p_value >= 0
            assert p_value <= 1  # P-values should be between 0 and 1

    assert has_pvalue, "No GWAS P-values found in any records"

    # Find the record that has Variant = 3_37014530_T/C
    for record in data["data"]:
        if record["Variant"] == "3:37014530_T/C":
            assert record["GWAS_P"] == 2.34e-3
            assert record["Variant"] == "3:37014530_T/C"
            assert record["GWAS_P"] == 2.34e-3
            assert record["OR"] == 7.81
            assert record["vs11"] == 1.0521
            assert record["Gene"] == "MLH1"
            assert record["HMM Score"] == "0.54"
            assert record["Classification"] == "4.2"


@pytest.fixture
def expected_gwas_data():
    """FUCK FUCK FUCK FUCK FUCK. </passive aggressive comment>"""
    import pandas as pd

    return pd.DataFrame(
        [
            {
                "Variant": "19:44908684_C/T",
                "Gene": "APOE",
                "RSID": "rs429358",
                "F_A": 0.3835,
                "F_U": 0.1517,
                "OR": 3.479,
                "P": 0.0,
            },
            {
                "Variant": "19:44905910_C/G",
                "Gene": None,
                "RSID": "rs440446",
                "F_A": 0.2718,
                "F_U": 0.3584,
                "OR": 0.6683,
                "P": 2.007e-54,
            },
            {
                "Variant": "19:44908822_T/C",
                "Gene": "APOE",
                "RSID": "rs7412",
                "F_A": 0.0392,
                "F_U": 0.08032,
                "OR": 0.4672,
                "P": 1.087e-44,
            },
            {
                "Variant": "19:44813331_A/G",
                "Gene": None,
                "RSID": "rs28399654",
                "F_A": 0.01968,
                "F_U": 0.03323,
                "OR": 0.584,
                "P": 3.583e-13,
            },
            {
                "Variant": "19:44793549_T/C",
                "Gene": None,
                "RSID": "Affx-16018598",
                "F_A": 0.02189,
                "F_U": 0.03554,
                "OR": 0.6074,
                "P": 1.996e-12,
            },
        ]
    )


def test_phecode_gwas_pvalues_explicitly(
    auth_integration_app_client, expected_gwas_data
):
    """Test that GWAS P-values are present in phecode term page for a specific case."""
    phecode = "290.11"

    # Now simulate clicking the GWAS button by calling the run-task endpoint
    # response = auth_integration_app_client.post(f"/run-task/{phecode}/1")
    response = auth_integration_app_client.post(f"/run-task/{phecode}/0")
    assert response.status_code == 200
    data = json.loads(response.data)

    # Verify GWAS completed successfully
    assert data["status"] == "completed"
    assert "associations" in data
    assert len(data["associations"]) == 83011

    import pandas as pd
    import numpy as np

    associations_df = pd.DataFrame(data["associations"])

    # Look for the expected rows in the associations dataframe
    for index, row in expected_gwas_data.iterrows():
        assert row["Variant"] in associations_df["Variant"].values

        assoc = associations_df[associations_df["Variant"] == row["Variant"]]
        assert np.isclose(assoc["P"].values[0], row["P"])
        assert assoc["OR"].values[0] == row["OR"]
        assert assoc["F_A"].values[0] == row["F_A"]
        assert assoc["F_U"].values[0] == row["F_U"]
        assert assoc["Gene"].values[0] == row["Gene"]
        # The table has a link
        assert row["RSID"] in assoc["RSID"].values[0]
