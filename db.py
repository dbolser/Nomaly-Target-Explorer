import logging
from typing import Dict, List

import mysql.connector
import pandas as pd
from mysql.connector.abstracts import MySQLConnectionAbstract

from config import Config
from errors import DatabaseConnectionError, DataNotFoundError

logger = logging.getLogger(__name__)


def get_db_connection() -> MySQLConnectionAbstract:
    """Get a database connection using config values."""
    print(
        f"DEBUG DB_PY: MYSQL_HOST={Config.MYSQL_HOST}, MYSQL_PORT={Config.MYSQL_PORT}, MYSQL_DB={Config.MYSQL_DB}, MYSQL_USER={Config.MYSQL_USER}"
    )  # DEBUG LINE
    try:
        conn = mysql.connector.connect(
            host=Config.MYSQL_HOST,
            port=Config.MYSQL_PORT,
            user=Config.MYSQL_USER,
            # password=Config.MYSQL_PASSWORD,
            unix_socket="/var/run/mysqld/mysqld.sock",
            database=Config.MYSQL_DB,
        )
        assert isinstance(conn, MySQLConnectionAbstract)
        return conn
    except Exception as e:
        logger.error("Database connection failed", exc_info=True)
        raise DatabaseConnectionError(f"Could not connect to database: {str(e)}")


def get_all_phecodes() -> pd.DataFrame:
    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                cur.execute(
                    """
                    SELECT 
                      phecode, description, sex, phecode_group
                    FROM phecode_definition
                    """
                )
                results = cur.fetchall()

                if not results:
                    raise DataNotFoundError("No data found for all phecodes")

                assert cur.description is not None
                columns = [desc[0] for desc in cur.description]

                filtered_df = pd.DataFrame(results, columns=columns)

                return filtered_df
    except DatabaseConnectionError:
        raise
    except Exception as e:
        logger.error(f"Error fetching all phecodes: {str(e)}", exc_info=True)
        raise


def get_phecode_info(phecode: str) -> dict:
    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                cur.execute(
                    "SELECT * FROM phecode_definition WHERE phecode = %s", (phecode,)
                )
                results = cur.fetchone()

                if not results:
                    raise DataNotFoundError(f"No data found for phecode: {phecode}")

                assert isinstance(results, dict)
                return dict(results)
    except DatabaseConnectionError:
        raise
    except Exception as e:
        logger.error(f"Error fetching phecode info: {str(e)}", exc_info=True)
        raise


def get_term_names(terms_list: List[str]) -> Dict[str, str]:
    """Get the names of a list of terms."""

    if not terms_list:
        logger.warning("No terms provided")
        return {}

    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                placeholders = ", ".join(f"'{x}'" for x in terms_list)
                query = (
                    f"SELECT term, name FROM term_info WHERE term IN ({placeholders})"
                )

                cur.execute(query)
                results = cur.fetchall()

                if not results:
                    raise DataNotFoundError(f"No data found for terms: {terms_list}")

                # Result is a list of dictionaries
                term_name_dict = {result["term"]: result["name"] for result in results}

                return term_name_dict
    except DatabaseConnectionError:
        raise


def get_term_domains(terms_list: list[str]) -> dict[str, set[str]]:
    if not terms_list:
        logger.warning("No terms provided")
        return {}

    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                placeholders = ", ".join(f"'{x}'" for x in terms_list)
                query = f"SELECT term, domain, domaintype FROM term_domain WHERE term IN ({placeholders})"

                # filter the phecodes and the ICD10 table based on the query
                cur.execute(query)
                results = cur.fetchall()

                # turn into a dictionary, Note results is a list of dictionaries
                term_domain_dict = {}
                for result in results:
                    term = result["term"]
                    domain = result["domain"]
                    domaintype = result["domaintype"]

                    if term not in term_domain_dict:
                        term_domain_dict[term] = set()
                    # add domain with link
                    if domaintype == "sf" or domaintype == "fa":
                        # TODO: Don't do this in the DB module, separate concerns!
                        term_domain_dict[term].add(
                            f'<a href="https://supfam.org/SUPERFAMILY/cgi-bin/scop.cgi?sunid={domain}">{domain}</a>'
                        )
                    elif domain.startswith("PF"):
                        term_domain_dict[term].add(
                            f'<a href="https://www.ebi.ac.uk/interpro/entry/pfam/{domain}">{domain}</a>'
                        )
                    else:
                        term_domain_dict[term].add(domain)

        return term_domain_dict
    except DatabaseConnectionError:
        raise
    except Exception as e:
        logger.error(f"Error fetching term domains: {str(e)}", exc_info=True)
        raise


def get_term_genes(terms_list: list[str]) -> pd.DataFrame:
    if not terms_list:
        logger.warning("No terms provided")
        return pd.DataFrame()

    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                placeholders = ", ".join(f"'{x}'" for x in terms_list)

                query = f"""
                SELECT DISTINCT td.term, vc.gene
                FROM term_domain td
                JOIN variants_consequences vc ON td.domain = vc.sf
                WHERE td.domaintype = 'sf' AND td.term IN ({placeholders})
                """
                cur.execute(query)
                results = cur.fetchall()

                # dataframe
                term_gene_df = pd.DataFrame(results, columns=["term", "gene"])

                # get the FA domains
                query = f"""
                    SELECT DISTINCT td.term, vc.gene
                    FROM term_domain td
                    JOIN variants_consequences vc ON td.domain = vc.fa
                    WHERE td.domaintype = 'fa' AND td.term IN ({placeholders})
                """
                cur.execute(query)
                results = cur.fetchall()

                # dataframe
                term_gene_df = pd.concat(
                    [term_gene_df, pd.DataFrame(results, columns=["term", "gene"])],
                    axis=0,
                )

                # drop duplicates
                term_gene_df.drop_duplicates(inplace=True)

                return term_gene_df
    except DatabaseConnectionError:
        raise
    except Exception as e:
        logger.error(f"Error fetching term genes: {str(e)}", exc_info=True)
        raise


def get_term_domain_genes_variant(term: str) -> pd.DataFrame:
    if not term:
        logger.warning("No term provided")
        return pd.DataFrame()

    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                # vc.aa, vc.hmm, vc.hmm_pos, vc.sf, vc.fa
                query = f"""
                SELECT DISTINCT vc.variant_id, vc.gene
                FROM terms2snps ts
                JOIN variants_consequences vc ON ts.variant_id = vc.variant_id
                WHERE ts.term = '{term}'
                """
                cur.execute(query)
                results = cur.fetchall()

                if not results:
                    raise DataNotFoundError(f"No data found for term: {term}")

                assert cur.description is not None
                columns = [desc[0] for desc in cur.description]

                # dataframe
                term_df = pd.DataFrame(results, columns=columns)

                return term_df
    except DatabaseConnectionError:
        raise
    except Exception as e:
        logger.error(f"Error fetching term domain genes: {str(e)}", exc_info=True)
        raise


def get_term_domain_genes(term: str) -> pd.DataFrame:
    if not term:
        logger.warning("No term provided")
        return pd.DataFrame()

    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                # vc.aa, vc.hmm, vc.hmm_pos, vc.sf, vc.fa
                query = f"""
                SELECT DISTINCT vc.gene
                FROM terms2snps ts
                JOIN variants_consequences vc ON ts.variant_id = vc.variant_id
                WHERE ts.term = '{term}'
                """
                cur.execute(query)
                results = cur.fetchall()

                if not results:
                    raise DataNotFoundError(f"No data found for term: {term}")

                assert cur.description is not None
                columns = [desc[0] for desc in cur.description]

                # dataframe
                term_df = pd.DataFrame(results, columns=columns)

                return term_df
    except DatabaseConnectionError:
        raise
    except Exception as e:
        logger.error(f"Error fetching term domain genes: {str(e)}", exc_info=True)
        raise


def get_variant_gene(variant_id: str) -> str | None:
    """Get gene information for a variant.

    Args:
        variant_id (str): The variant ID in any format (will be normalized)

    Returns:
        str | None: The gene name if found, None otherwise
    """
    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                # Normalize variant ID format for database lookup
                parts = variant_id.replace(":", "_").split("_")
                if len(parts) == 4:  # Format: chr_pos_ref_alt
                    normalized_variant = "_".join(parts)
                else:  # Format: chr_pos_ref/alt
                    normalized_variant = "_".join(
                        [
                            parts[0],
                            parts[1],
                            parts[2].split("/")[0],
                            parts[2].split("/")[1],
                        ]
                    )

                # Updated query to match how we get gene info in other places
                query = """
                SELECT DISTINCT gene
                FROM variants_consequences
                WHERE variant_id = %s
                """
                cur.execute(query, (normalized_variant,))
                results = cur.fetchall()

                if not results:
                    return None

                # If we have multiple genes, join them with commas
                genes = [r["gene"] for r in results if r["gene"]]
                return ", ".join(genes) if genes else None

    except Exception as e:
        logger.error(
            f"Error fetching gene for variant {variant_id}: {str(e)}", exc_info=True
        )
        return None


def get_term_variants(term: str) -> pd.DataFrame:
    """
    Return a DataFrame of variant_id, gene, aa, hmm_score for the given term.
    """
    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                query = """
                    SELECT
                        term,
                        variant_id,
                        gene,
                        GROUP_CONCAT(DISTINCT aa) AS aa,
                        MAX(ABS(wild - mutant)) as hmm_score
                    FROM
                        variants_consequences
                    INNER JOIN
                        terms2snps
                    USING
                        (variant_id)
                    WHERE
                        term = %s
                    GROUP BY
                        term,
                        variant_id,
                        gene
                """
                logger.info(f"\nExecuting query for term: {term}")
                logger.debug(f"Query: {query}")

                cur.execute(query, (term,))
                results = cur.fetchall()
                logger.info(f"Number of results from DB: {len(results)}")

                if not results:
                    logger.warning(f"No data found for term: {term}")
                    raise DataNotFoundError(f"No data found for term: {term}")

                assert cur.description is not None
                columns = [desc[0] for desc in cur.description]

                df = pd.DataFrame(results, columns=columns)

                logger.debug(f"Created DataFrame with shape: {df.shape}")

                # Convert numpy types to Python native types
                if not df.empty:
                    df["hmm_score"] = df["hmm_score"].astype(float)
                    logger.debug("Sample of data:")
                    logger.debug(df.head())

                return df
    except DatabaseConnectionError:
        raise
    except Exception as e:
        logger.error(f"Error fetching term variants: {str(e)}", exc_info=True)
        raise


def get_all_variants() -> pd.DataFrame:
    try:
        with get_db_connection() as conn:
            with conn.cursor(dictionary=True) as cur:
                query = """
                SELECT DISTINCT variant_id FROM variants_consequences
                """
                cur.execute(query)
                results = cur.fetchall()
                return pd.DataFrame(results, columns=["variant_id"])
    except DatabaseConnectionError:
        raise
    except Exception as e:
        logger.error(f"Error fetching all variants: {str(e)}", exc_info=True)
        raise
