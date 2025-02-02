import pytest
import pandas as pd

from db import (
    get_db_connection,
    get_all_phecodes,
    get_phecode_info,
    get_term_names,
    get_term_domains,
    get_term_genes,
    get_term_domain_genes_variant,
    get_term_domain_genes,
    get_term_variants,
    get_all_variants,
)
from errors import DataNotFoundError


def test_db_connection():
    """Test that we can connect to the database and execute a simple query."""
    conn = get_db_connection()
    with conn.cursor() as cur:
        cur.execute("SELECT 1")
        result = cur.fetchone()
        assert result is not None
        assert result[0] == 1


def test_get_all_phecodes():
    """Test fetching all phecodes and verify the dataframe structure."""
    df = get_all_phecodes()

    # Check dataframe structure
    assert isinstance(df, pd.DataFrame)
    assert not df.empty, "Phecode dataframe should not be empty"

    required_columns = ["phecode", "description", "sex", "phecode_group"]
    assert all(col in df.columns for col in required_columns), (
        f"Missing required columns. Expected {required_columns}"
    )

    # Verify data quality
    assert df["phecode"].is_unique, "Phecodes should be unique"
    assert not df["phecode"].isnull().any(), "No null phecodes should exist"
    assert not df["description"].isnull().any(), "No null descriptions should exist"


@pytest.mark.parametrize(
    "phecode,expected_exists",
    [
        ("250.2", True),  # Type 2 diabetes - should exist
        ("999.999", False),  # Invalid phecode - should not exist
    ],
)
def test_get_phecode_info(phecode, expected_exists):
    """Test fetching specific phecode info with different inputs."""
    if expected_exists:
        info = get_phecode_info(phecode)
        assert isinstance(info, dict)
        assert info["phecode"] == phecode
        assert "description" in info
        assert "sex" in info
    else:
        with pytest.raises(DataNotFoundError):
            get_phecode_info(phecode)


@pytest.mark.parametrize(
    "terms",
    [
        (["GO:0005524", "GO:0005525"]),  # ATP binding, GTP binding
        ([]),  # Empty list
        (["INVALID:123"]),  # Invalid term
    ],
)
def test_get_term_names(terms):
    """Test fetching term names with various inputs."""
    names = get_term_names(terms)
    assert isinstance(names, dict)

    if not terms:  # Empty input
        assert len(names) == 0
        return

    if "INVALID:123" in terms:  # Invalid input
        assert "INVALID:123" not in names
        return

    # Valid terms
    assert all(term in names for term in terms), (
        "All valid terms should have corresponding names"
    )
    assert all(isinstance(name, str) for name in names.values()), (
        "All names should be strings"
    )


def test_get_term_domains():
    """Test fetching term domains with known term."""
    terms = ["GO:0005524"]  # ATP binding
    domains = get_term_domains(terms)

    # Check structure
    assert isinstance(domains, dict)
    assert len(domains) > 0
    assert all(term in domains for term in terms)

    # Check domain sets
    assert all(isinstance(domain_set, set) for domain_set in domains.values())

    # Check content for known term
    atp_domains = domains["GO:0005524"]
    assert len(atp_domains) > 0, "ATP binding term should have associated domains"


def test_get_term_genes():
    """Test fetching term genes and verify dataframe structure."""
    terms = ["GO:0005524"]  # ATP binding
    df = get_term_genes(terms)

    # Check dataframe structure
    assert isinstance(df, pd.DataFrame)
    assert not df.empty, "Result should not be empty for ATP binding term"
    assert all(col in df.columns for col in ["term", "gene"])

    # Verify data quality
    assert not df["term"].isnull().any(), "No null terms should exist"
    assert not df["gene"].isnull().any(), "No null genes should exist"


def test_get_term_domain_genes_variant():
    """Test fetching term domain genes with variants."""
    term = "GO:0005524"  # ATP binding
    df = get_term_domain_genes_variant(term)

    # Check dataframe structure
    assert isinstance(df, pd.DataFrame)
    assert not df.empty, "Result should not be empty for ATP binding term"
    assert all(col in df.columns for col in ["variant_id", "gene"])

    # Verify data quality
    assert not df["variant_id"].isnull().any(), "No null variant IDs should exist"
    assert not df["gene"].isnull().any(), "No null genes should exist"


def test_get_term_domain_genes():
    """Test fetching term domain genes."""
    # with pytest.raises(DataNotFoundError):
    #    get_term_domain_genes(term)

    term = "DOES_NOT_EXIST"

    with pytest.raises(DataNotFoundError):
        get_term_domain_genes(term)

    term = "GO:0005524"  # ATP binding

    df = get_term_domain_genes(term)

    # Check dataframe structure
    assert isinstance(df, pd.DataFrame)
    assert not df.empty, "Result should not be empty for ATP binding term"
    assert all(col in df.columns for col in ["gene"])


def test_get_term_variants():
    """Test fetching term variants."""
    term = "GO:0005524"  # ATP binding
    df = get_term_variants(term)

    # Check dataframe structure
    assert isinstance(df, pd.DataFrame)
    assert not df.empty, "Result should not be empty for ATP binding term"
    assert all(col in df.columns for col in ["term", "variant_id"])


def test_get_all_variants():
    """Test fetching all variants."""
    df = get_all_variants()
    assert isinstance(df, pd.DataFrame)
    assert not df.empty, "Variants dataframe should not be empty"
    assert all(col in df.columns for col in ["variant_id"])


def test_edge_cases():
    """Test edge cases and error handling."""
    # Test with empty inputs
    assert len(get_term_names([])) == 0
    assert len(get_term_domains([])) == 0
    assert len(get_term_genes([])) == 0

    # Test with invalid inputs
    invalid_terms = ["INVALID:123"]
    assert len(get_term_names(invalid_terms)) == 0
    assert len(get_term_domains(invalid_terms)) == 0
    assert get_term_genes(invalid_terms).empty
    assert get_term_variants(invalid_terms[0]).empty

    # Test with large input
    many_terms = [
        f"GO:{i:07d}" for i in range(100)
    ]  # Reduced from 1000 to be more reasonable
    names = get_term_names(many_terms)
    assert isinstance(names, dict)
