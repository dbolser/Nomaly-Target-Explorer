import mysql.connector
import pandas as pd
from errors import DatabaseConnectionError, DataNotFoundError
import logging
from flask import current_app
from utils import app_context

logger = logging.getLogger(__name__)


# Connect to the database
def get_db_connection():
    with app_context():
        try:
            conn = mysql.connector.connect(
                host=current_app.config['MYSQL_HOST'],
                user=current_app.config['MYSQL_USER'],
                password=current_app.config['MYSQL_PASSWORD'],
                database=current_app.config['MYSQL_DB']
            )
            return conn
        except Exception as e:
            logger.error("Database connection failed", exc_info=True)
            raise DatabaseConnectionError(f"Could not connect to database: {str(e)}")


def get_all_phecodes() -> pd.DataFrame:
    conn = get_db_connection()
    cur = conn.cursor()

    # filter the phecodes and the ICD10 table based on the query
    cur.execute(
        """
        SELECT 
          phecode, description, sex, phecode_group #, 
          #GROUP_CONCAT(icd10 ORDER BY icd10_count DESC),
          #GROUP_CONCAT(meaning ORDER BY icd10_count DESC), 
          #GROUP_CONCAT(icd10_count ORDER BY icd10_count DESC), 
          #SUM(icd10_count)
        FROM phecode_definition
        #INNER JOIN icd10_phecode USING(phecode)
        #INNER JOIN icd10_coding USING(icd10)
        #GROUP BY phecode
        ;
        """
    )
    results = cur.fetchall()

    # cur.close()
    # conn.close()

    # 4. Get the column names
    columns = [desc[0] for desc in cur.description]

    # 7. Convert the rows and columns into a Pandas DataFrame
    filtered_df = pd.DataFrame(results, columns=columns)

    return filtered_df


def get_phecode_info(phecode):
    try:
        conn = get_db_connection()
        cur = conn.cursor(dictionary=True)

        cur.execute("SELECT * FROM phecode_definition WHERE phecode = %s", (phecode,))
        result = cur.fetchone()

        if not result:
            raise DataNotFoundError(f"No data found for phecode: {phecode}")

        return result
    except DatabaseConnectionError:
        raise
    except Exception as e:
        logger.error(f"Error fetching phecode info: {str(e)}", exc_info=True)
        raise
    finally:
        if "cur" in locals():
            cur.close()
        if "conn" in locals():
            conn.close()


def get_term_names(terms_list):
    conn = get_db_connection()
    cur = conn.cursor()

    placeholders = ", ".join(f"'{x}'" for x in terms_list)
    query = f"SELECT term, name FROM term_info WHERE term IN ({placeholders})"

    # filter the phecodes and the ICD10 table based on the query
    cur.execute(query)
    results = cur.fetchall()

    cur.close()
    conn.close()

    # turn into a dictionary
    term_name_dict = {term: name for term, name in results}

    return term_name_dict

def get_term_domains(terms_list):
    conn = get_db_connection()
    cur = conn.cursor()

    placeholders = ", ".join(f"'{x}'" for x in terms_list)
    query = f"SELECT term, domain, domaintype FROM term_domain WHERE term IN ({placeholders})"

    # filter the phecodes and the ICD10 table based on the query
    cur.execute(query)
    results = cur.fetchall()

    cur.close()
    conn.close()

    # turn into a dictionary
    term_domain_dict = {}
    for term, domain, domaintype in results:
        if term not in term_domain_dict:
            term_domain_dict[term] = set()
        # add domain with link
        if domaintype == 'sf' or domaintype == 'fa':
            term_domain_dict[term].add(f'<a href="https://supfam.org/SUPERFAMILY/cgi-bin/scop.cgi?sunid={domain}">{domain}</a>')
        elif domain.startswith('PF'):
            term_domain_dict[term].add(f'<a href="https://www.ebi.ac.uk/interpro/entry/pfam/{domain}">{domain}</a>')
        else:
            term_domain_dict[term].add(domain)

    return term_domain_dict

def get_term_genes(terms_list):
    conn = get_db_connection()
    cur = conn.cursor()

    placeholders = ", ".join(f"'{x}'" for x in terms_list)

    query = f"""
    SELECT DISTINCT td.term, vc.gene
    FROM term_domain td
    JOIN variants_consequences vc ON td.domain = vc.sf
    WHERE td.domaintype = 'sf' AND td.term IN ({placeholders})
    """
    # filter the phecodes and the ICD10 table based on the query
    cur.execute(query)
    results = cur.fetchall()

    # dataframe
    term_gene_df = pd.DataFrame(results, columns=['term', 'gene'])
    
    query = f"""
    SELECT DISTINCT td.term, vc.gene
    FROM term_domain td
    JOIN variants_consequences vc ON td.domain = vc.fa
    WHERE td.domaintype = 'fa' AND td.term IN ({placeholders})
    """
    cur.execute(query)
    results = cur.fetchall()

    cur.close()
    conn.close()

    # dataframe
    term_gene_df = pd.concat([term_gene_df, pd.DataFrame(results, columns=['term', 'gene'])], axis=0)

    # drop duplicates
    term_gene_df.drop_duplicates(inplace=True)

    return term_gene_df

def get_term_domain_genes_variant(term):
    conn = get_db_connection()
    cur = conn.cursor()

    # vc.aa, vc.hmm, vc.hmm_pos, vc.sf, vc.fa
    query = f"""
    SELECT DISTINCT vc.variant_id, vc.gene
    FROM terms2snps ts
    JOIN variants_consequences vc ON ts.variant_id = vc.variant_id
    WHERE ts.term = '{term}'
    """
    # filter the phecodes and the ICD10 table based on the query
    cur.execute(query)
    results = cur.fetchall()

    columns = [desc[0] for desc in cur.description]

    cur.close()

    # dataframe
    term_df = pd.DataFrame(results, columns=columns)

    return term_df

def get_term_domain_genes(term):
    conn = get_db_connection()
    cur = conn.cursor()

    # vc.aa, vc.hmm, vc.hmm_pos, vc.sf, vc.fa
    query = f"""
    SELECT DISTINCT vc.gene
    FROM terms2snps ts
    JOIN variants_consequences vc ON ts.variant_id = vc.variant_id
    WHERE ts.term = '{term}'
    """
    # filter the phecodes and the ICD10 table based on the query
    cur.execute(query)
    results = cur.fetchall()

    columns = [desc[0] for desc in cur.description]

    cur.close()

    # dataframe
    term_df = pd.DataFrame(results, columns=columns)

    return term_df
