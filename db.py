import mysql.connector
from mysql.connector import Error
import pandas as pd

# Connect to the database
def get_db_connection():
    # Database connection configuration
    db_config = {
        'host': 'localhost',
        'database': 'ukbb',
        'user': 'clu',
        'password': 'mysqlOutse3',
        'raise_on_warnings': False
    }
    try:
        conn = mysql.connector.connect(**db_config)
        return conn
    except Exception as e:
        print(f"Error connecting to database: {e}")
        return None



def get_all_phecodes() -> pd.DataFrame:

    conn = get_db_connection()
    cur = conn.cursor()

    # filter the phecodes and the ICD10 table based on the query
    cur.execute(
        f"""
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

    #cur.close()
    #conn.close()

    # 4. Get the column names
    columns = [desc[0] for desc in cur.description]

    # 7. Convert the rows and columns into a Pandas DataFrame
    filtered_df = pd.DataFrame(results, columns=columns)

    # Print the first 5 rows
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(filtered_df.head())

    return filtered_df


def get_phecode_info(phecode):

    conn = get_db_connection()
    cur = conn.cursor()

    # filter the phecodes and the ICD10 table based on the query
    cur.execute(
        f"""
        SELECT p.phecode, p.description, p.sex, p.phecode_group, p.phecode_exclude, p.affected, p.excluded, icd10.icd10, icd10.meaning, icd10.icd10_count
        FROM phecode_definition p
        JOIN icd10_phecode ip ON p.phecode = ip.phecode
        JOIN icd10_coding icd10 ON ip.icd10 = icd10.icd10
        WHERE p.phecode = '{phecode}';
        """
    )
    results = cur.fetchall()

    cur.close()
    conn.close()

    # 4. Get the column names
    columns = [desc[0] for desc in cur.description]

    # 7. Convert the rows and columns into a Pandas DataFrame
    filtered_df = pd.DataFrame(results, columns=columns)

    # Add icd10_count to meaning
    filtered_df['meaning'] = filtered_df['meaning'] + ' | ' + filtered_df['icd10_count'].astype(str)

    # Group by 'description' and aggregate meanings into a list
    grouped = filtered_df.groupby(['description', 'phecode', 'sex', 'affected', 'excluded', 'phecode_exclude', 'phecode_group'
    ]).agg({
        'meaning': lambda x: list(x)  # Group meanings into a list
    }).reset_index()

    # Convert results to a list of dictionaries
    results = grouped.to_dict(orient='records')

    # results only has one row
    data = results[0]

    # # Example data (replace with actual data retrieval logic)
    # data = {
    #     "phecode": phecode,
    #     "sex": "Both",
    #     "affected": 150,
    #     "excluded": 10,
    #     "phecode_exclude": "None",
    #     "phecode_group": "Cardiovascular"
    # }

    return data

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