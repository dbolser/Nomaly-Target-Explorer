import mysql.connector
from mysql.connector import Error


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

# def get_db_connection():
#     connection = psycopg2.connect(
#         host="localhost",
#         database="ukbb_eur",
#     )
#     return connection