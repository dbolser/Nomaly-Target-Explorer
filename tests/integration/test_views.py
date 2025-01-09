import json
import time

# The client is created in conftest.py
# from conftest import client, auth_client


def test_index_route_unauthenticated(client):
    """Test the index route redirects to login when not authenticated."""
    response = client.get("/anything")
    assert response.status_code == 302  # Redirect to login
    assert "/login" in response.location


def test_index_route_authenticated(auth_client):
    """Test the index route when authenticated."""
    response = auth_client.get("/")
    assert response.status_code == 200
    assert b"welcome" in response.data.lower()


def test_login_route(client, test_admin):
    """Test login functionality."""
    # Test GET request
    response = client.get("/login")
    assert response.status_code == 200
    assert b"login" in response.data.lower()

    # Test POST request with valid credentials
    response = client.post(
        "/login",
        data={"username": test_admin["username"], "password": test_admin["password"]},
    )
    assert response.status_code == 302  # Should be a redirect
    assert response.location == "/"  # Should redirect to index

    response = client.get("/logout")
    assert response.status_code == 302  # Should be a redirect
    assert response.location == "/"  # Should redirect to index

    # Test POST request with invalid credentials
    response = client.post(
        "/login",
        data={"username": test_admin["username"], "password": "wr0nG_pa55woRd"},
    )
    assert response.status_code == 200  # Should stay on login page
    assert b"incorrect username or password" in response.data.lower()


def test_search_route(auth_client):
    """Test the disease search functionality."""
    response = auth_client.get("/diseasesearch?q=diabetes")
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


def test_phecode_route(auth_client):
    """Test the phecode detail route."""
    test_phecode = "250.2"  # Example phecode
    response = auth_client.get(f"/phecode/{test_phecode}")
    assert response.status_code == 200
    assert bytes(test_phecode, "utf-8") in response.data


def test_page1_structure(auth_client):
    """Test the structure and content of the page1 route."""
    response = auth_client.get("/page1")
    assert response.status_code == 200

    # Convert response data to string for easier testing
    html = response.data.decode("utf-8")

    # Test main heading
    assert '<h1 class="text-center">Selected diseases</h1>' in html

    # Test category headings (H2s)
    assert '<h2 class="mt-5 mb-4 ps-2">Skin Conditions</h2>' in html
    assert '<h2 class="mt-5 mb-4 ps-2">Women&#39;s Health</h2>' in html

    # Test result sections (H3s)
    assert 'Results for "Hidradenitis"' in html
    assert 'Results for "Polycystic Ovarian"' in html

    # Test description text under each search
    assert "Any Phecode or ICD10 descriptions that contains the word" in html

    # Test dynamic content loading
    assert "results-list-1" in html  # First results list
    assert "results-list-2" in html  # Second results list

    # Test JavaScript initialization
    assert "searchData(" in html
    assert "async function searchData(query, listIndex)" in html


def test_page1_search_results(auth_client):
    """Test that search results are properly structured when loaded."""
    # First make a search request to get some results
    response = auth_client.get("/diseasesearch?query=Hidradenitis")
    assert response.status_code == 200
    data = json.loads(response.data)

    # Verify we got some results
    assert len(data) > 0

    # Check first result has required structure
    first_result = data[0]
    assert "phecode" in first_result
    assert "description" in first_result
    assert "sex" in first_result
    assert "affected" in first_result
    assert "excluded" in first_result
    assert "phecode_exclude" in first_result
    assert "phecode_group" in first_result


def test_phecode_term_structure(auth_client):
    """Test the structure and content of a specific phecode term page."""
    phecode = "649.1"
    term = "GO:0035235"

    response = auth_client.get(f"/phecode/{phecode}/term/{term}")
    assert response.status_code == 200

    html = response.data.decode("utf-8")

    # Test phecode section
    assert "Diabetes or abnormal glucose tolerance complicating pregnancy" in html
    assert f'<span><a href="/phecode/{phecode}"' in html
    assert "Sex: Female" in html
    assert "<span>Affected: <strong>300</strong></span>" in html
    assert "Excluded: 138" in html
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


def test_phecode_term_variant_detail(auth_client):
    """Test the JSON response from the variant detail endpoint."""
    phecode = "649.1"
    term = "GO:0035235"

    response = auth_client.get(f"/phecode/{phecode}/term/{term}/tableVariantDetail")
    assert response.status_code == 200

    data = json.loads(response.data)

    # Check structure of response
    assert "data" in data
    assert "columns" in data
    assert "defaultColumns" in data
    assert "numColumns" in data

    # Check expected columns exist
    expected_columns = [
        "Variant",
        "Gene",
        "AA_Change",
        "HMM_Score",
        "Classification",
        "Drug_Program_Indication",
        "TDL",
        "TBIO",
        "vs00",
        "vs01",
        "vs11",
        "hmoz_alt",
        "hmoz_ref",
        "htrz",
        "GWAS_P",
        "GWAS_OR",
        "GWAS_RSID",
    ]
    assert all(col in data["columns"] for col in expected_columns)

    # Check data records if any exist
    if len(data["data"]) > 0:
        first_record = data["data"][0]
        # Check required fields in first record
        assert "Variant" in first_record
        assert "Gene" in first_record
        assert "HMM_Score" in first_record
        assert "GWAS_P" in first_record

        # Check numeric formatting
        assert float(first_record["HMM_Score"]) >= 0
        assert float(first_record["GWAS_P"]) >= 0


def test_phecode_page_structure(auth_client):
    """Test the structure and content of a specific phecode page."""
    phecode = "649.1"

    response = auth_client.get(f"/phecode/{phecode}")
    assert response.status_code == 200

    html = response.data.decode("utf-8")

    # Test main content
    assert "Diabetes or abnormal glucose tolerance complicating pregnancy" in html
    assert f"Phecode: {phecode} (Run v1)" in html
    assert "Sex: Female" in html
    assert "Affected: <strong>300</strong>" in html
    assert "Control: 263177" in html
    assert "Excluded: 138 (649-649.99)" in html
    assert "Disease: pregnancy complications" in html

    # Test Nomaly Results section
    assert "Per Term Nomaly Results" in html
    assert '<table id="resultsTable"' in html

    # Test JavaScript initialization
    assert f'const phecode = "{phecode}";' in html
    assert 'const runbatch = "Run v1";' in html


def test_phecode_page_with_gwas(auth_client):
    """Test the phecode page with GWAS functionality enabled."""
    phecode = "649.1"

    response = auth_client.get(f"/phecode/{phecode}?gwas=1")
    assert response.status_code == 200

    html = response.data.decode("utf-8")

    # Test GWAS button exists
    assert (
        '<button id="runTaskButton" class="btn btn-primary">GWAS Task</button>' in html
    )

    # Test GWAS results containers exist
    assert '<div id="taskResult"' in html
    assert '<div id="GWAStableContainer"' in html
    assert '<table id="variantTable"' in html

    # Test GWAS table headers
    assert "<th>Variant</th>" in html
    assert "<th>Gene</th>" in html
    assert "<th>RSID</th>" in html
    assert "<th>F_A</th>" in html
    assert "<th>F_U</th>" in html
    assert "<th>OR</th>" in html
    assert "<th>P</th>" in html


# def test_phecode_gwas_task_result(client):
#     """Test the GWAS task result endpoint."""
#     phecode = "649.1"

#     response = client.get(f"/task-result/{phecode}")
#     assert response.status_code == 200

#     data = json.loads(response.data)

#     # Check response structure
#     assert "result" in data
#     assert "associations" in data

#     # If GWAS has been run, check the results
#     if "GWAS identified" in data["result"]:
#         assert "2644 missense variants" in data["result"]
#         assert "1963 unique genes" in data["result"]

#         # Check associations data
#         associations = data["associations"]
#         assert len(associations) > 0
#         if len(associations) > 0:
#             first_assoc = associations[0]
#             assert "Variant" in first_assoc
#             assert "Gene" in first_assoc
#             assert "RSID" in first_assoc
#             assert "P" in first_assoc
#             assert "OR" in first_assoc


def test_phecode_nomaly_stats(auth_client):
    """Test the Nomaly stats endpoint."""
    phecode = "649.1"

    # Test both v1 and v2 endpoints
    for version in ["nomaly-stats", "nomaly-stats2"]:
        response = auth_client.post(f"/{version}/{phecode}")
        assert response.status_code == 200

        data = json.loads(response.data)

        # Check response structure
        assert "qqplot" in data
        assert "data" in data
        assert "columns" in data
        assert "columnNames" in data
        assert "defaultColumns" in data
        assert "numColumns" in data

        # Check expected columns exist
        expected_columns = [
            "minrank",
            "term",
            "name",
            "domain",
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
            assert "domain" in first_record
            assert "minrank" in first_record

            # Check p-value formatting
            for pval_field in ["mwu_pvalue", "metric1_pvalue", "mcc_pvalue"]:
                if pval_field in first_record:
                    pval = first_record[pval_field]
                    if pval != "nan":
                        # Should be in scientific notation
                        assert "e" in pval.lower()

            # Check data limits
            assert len(data["data"]) <= 1000  # Should be limited to 1000 entries


def test_variant_page_structure(auth_client):
    """Test the structure and content of a specific variant page."""
    variant = "17_80117714_G_A"

    response = auth_client.get(f"/variant/{variant}")
    assert response.status_code == 200

    html = response.data.decode("utf-8")

    # Test page title and headers
    assert "<h3>Variant Details</h3>" in html
    assert '<h5 class="card-title">Basic Information</h5>' in html

    # Test variant information
    assert "<strong>Variant ID:</strong> 17:80117714_G/A" in html
    assert "<strong>Chromosome:</strong> 17" in html
    assert "<strong>Position:</strong> 80117714" in html
    assert "<strong>Ref Allele:</strong> G" in html
    assert "<strong>Alt Allele:</strong> A" in html

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


def test_variant_phewas_results(auth_client):
    """Test the PheWAS results endpoint for a variant."""
    variant = "17_80117714_G_A"

    # First check initial state
    response = auth_client.get(f"/phewas-result/{variant}")
    assert response.status_code == 200
    initial_data = json.loads(response.data)
    assert initial_data["result"] == "Processing..."
    assert initial_data["associations"] == []

    # Trigger PheWAS analysis
    response = auth_client.post(f"/run-phewas/{variant}")
    assert response.status_code == 202
    assert json.loads(response.data)["status"] == "Task started"

    # Poll for results with timeout
    max_wait = 30  # Maximum wait time in seconds
    start_time = time.time()
    data = None

    while time.time() - start_time < max_wait:
        response = auth_client.get(f"/phewas-result/{variant}")
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


def test_variant_id_formats(auth_client):
    """Test different variant ID format handling."""
    variants = [
        "17_80117714_G_A",  # Underscore format
    ]

    for variant in variants:
        response = auth_client.get(f"/variant/{variant}")
        assert response.status_code == 200
        html = response.data.decode("utf-8")
        assert "<strong>Chromosome:</strong> 17" in html
        assert "<strong>Position:</strong> 80117714" in html


def test_phecode_gwas_pvalues(auth_client):
    """Test that GWAS P-values are present in phecode page for a specific case."""
    phecode = "561"

    # First check the page with GWAS enabled
    response = auth_client.get(f"/phecode/{phecode}?gwas=1")
    assert response.status_code == 200
    html = response.data.decode("utf-8")
    assert (
        '<button id="runTaskButton" class="btn btn-primary">GWAS Task</button>' in html
    )

    # Now simulate clicking the GWAS button by calling the run-task endpoint
    response = auth_client.post(f"/run-task/{phecode}")
    assert response.status_code == 200
    data = json.loads(response.data)

    # Verify GWAS completed successfully
    assert data["status"] == "completed"
    assert "associations" in data
    assert len(data["associations"]) > 0

    # Check P-values in GWAS results
    for assoc in data["associations"]:
        assert "P" in assoc
        assert assoc["P"] is not None
        p_value = float(assoc["P"])
        assert p_value >= 0
        assert p_value <= 1  # P-values should be between 0 and 1

    # Also check the nomaly stats which should include GWAS data
    response = auth_client.post(f"/nomaly-stats/{phecode}")
    assert response.status_code == 200
    data = json.loads(response.data)

    # Verify we have data
    assert "data" in data
    assert len(data["data"]) > 0

    # Check that P-values are present and not empty
    for record in data["data"]:
        if "mwu_pvalue" in record:  # This is a term with stats
            assert "metric1_pvalue" in record
            assert record["metric1_pvalue"] is not None
            assert record["metric1_pvalue"] != ""
            # Convert to float to ensure it's a valid number


def test_phecode_term_gwas_pvalues(auth_client):
    """Test that GWAS P-values are present in phecode term page for a specific case."""
    phecode = "561"
    term = "HP:0000789"

    # Get the variant detail data
    response = auth_client.get(f"/phecode/{phecode}/term/{term}/tableVariantDetail")
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
