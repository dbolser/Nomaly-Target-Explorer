import pytest
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
)
from errors import DataNotFoundError
import pandas as pd


def test_db_connection():
    """Test that we can connect to the database"""
    conn = get_db_connection()
    assert conn is not None
    cur = conn.cursor()
    cur.execute("SELECT 1")
    result = cur.fetchone()
    assert result is not None
    assert result[0] == 1
    cur.close()
    conn.close()


def test_get_all_phecodes():
    """Test fetching all phecodes"""
    df = get_all_phecodes()
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert all(
        col in df.columns for col in ["phecode", "description", "sex", "phecode_group"]
    )
    # Check a known phecode exists
    assert len(df[df["phecode"] == "250.2"]) > 0  # Type 2 diabetes


def test_get_phecode_info():
    """Test fetching specific phecode info"""
    # Test with known phecode
    info = get_phecode_info("250.2")  # Type 2 diabetes
    assert isinstance(info, dict)
    assert info["phecode"] == "250.2"
    assert "description" in info
    assert "sex" in info

    # Test with invalid phecode
    with pytest.raises(DataNotFoundError):
        get_phecode_info("999.999")


def test_get_term_names():
    """Test fetching term names"""
    # Test with known terms
    terms = ["GO:0005524", "GO:0005525"]  # ATP binding, GTP binding
    names = get_term_names(terms)
    assert isinstance(names, dict)
    assert len(names) > 0
    assert all(term in names for term in terms)


def test_get_term_domains():
    """Test fetching term domains"""
    terms = ["GO:0005524"]  # ATP binding
    domains = get_term_domains(terms)
    assert isinstance(domains, dict)
    assert len(domains) > 0
    assert all(term in domains for term in terms)
    # Check that domains are returned as sets
    assert all(isinstance(domain_set, set) for domain_set in domains.values())


def test_get_term_genes():
    """Test fetching term genes"""
    terms = ["GO:0005524"]  # ATP binding
    df = get_term_genes(terms)
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert all(col in df.columns for col in ["term", "gene"])


def test_get_term_domain_genes_variant():
    """Test fetching term domain genes with variants"""
    term = "GO:0005524"  # ATP binding
    df = get_term_domain_genes_variant(term)
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert all(col in df.columns for col in ["variant_id", "gene"])


def test_get_term_domain_genes():
    """Test fetching term domain genes"""
    term = "GO:0005524"  # ATP binding
    df = get_term_domain_genes(term)
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert "gene" in df.columns


def test_get_term_variants():
    """Test fetching term variants"""
    term = "GO:0005524"  # ATP binding
    df = get_term_variants(term)
    assert isinstance(df, pd.DataFrame)
    if not df.empty:  # Some terms might not have variants
        assert all(
            col in df.columns
            for col in ["term", "variant_id", "gene", "aa", "hmm_score"]
        )


def test_edge_cases():
    """Test edge cases and error conditions"""
    # Empty list of terms
    assert len(get_term_names([])) == 0
    assert len(get_term_domains([])) == 0
    assert len(get_term_genes([])) == 0

    # Invalid terms
    invalid_terms = ["INVALID:123"]
    assert len(get_term_names(invalid_terms)) == 0
    assert len(get_term_domains(invalid_terms)) == 0
    assert get_term_genes(invalid_terms).empty
    assert get_term_variants(invalid_terms[0]).empty

    # Large number of terms
    many_terms = [f"GO:{i:07d}" for i in range(1000)]
    names = get_term_names(many_terms)
    assert isinstance(names, dict)
    # assert len(names) == 1000
