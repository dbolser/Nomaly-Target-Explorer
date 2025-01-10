from flask import request, jsonify, Blueprint
from db import get_db_connection
import pandas as pd

main_bp = Blueprint("main", __name__)

# # Route to serve the HTML page
# @main_bp.route('/', methods=['GET', 'POST'])
# def index():
#     return render_template('index.html')


# API route for search functionality in index.html
@main_bp.route("/diseasesearch")
def search():
    query = request.args.get("query", "").lower()

    results = []

    # Filter based on the query
    if query:
        # if query length is less than 3, return empty results
        if len(query) < 3:
            return results

        conn = get_db_connection()
        cur = conn.cursor()

        # Step 1: Search for matching Phecodes and the ICD10 table based on the query
        cur.execute(
            f"""
            SELECT p.phecode
            FROM phecode_definition p
            JOIN icd10_phecode ip ON p.phecode = ip.phecode
            JOIN icd10_coding icd10 ON ip.icd10 = icd10.icd10
            WHERE p.description LIKE '%{query}%'
                OR icd10.meaning LIKE '%{query}%';
            """
        )
        results = cur.fetchall()
        phecodes = [row[0] for row in results]

        # If no phecodes found, return an empty result
        if not phecodes:
            print("No matching Phecodes found.")
            return []

        # Step 2: Retrieve detailed information for the matching Phecodes
        placeholders = ", ".join(["%s"] * len(phecodes))
        details_query = f"""
        SELECT p.phecode, p.description, p.sex, p.phecode_group, p.phecode_exclude,
                p.affected, p.excluded, icd10.icd10, icd10.meaning, icd10.icd10_count
        FROM phecode_definition p
        JOIN icd10_phecode ip ON p.phecode = ip.phecode
        JOIN icd10_coding icd10 ON ip.icd10 = icd10.icd10
        WHERE p.phecode IN ({placeholders});
        """
        cur.execute(details_query, tuple(phecodes))
        results = cur.fetchall()

        cur.close()
        conn.close()

        # 4. Get the column names
        columns = [desc[0] for desc in cur.description]

        # 7. Convert the rows and columns into a Pandas DataFrame
        filtered_df = pd.DataFrame(results, columns=columns)

        # Add icd10_count to meaning
        filtered_df["meaning"] = (
            filtered_df["meaning"] + " | " + filtered_df["icd10_count"].astype(str)
        )

        # Group by 'description' and aggregate meanings into a list
        grouped = (
            filtered_df.groupby(
                [
                    "description",
                    "phecode",
                    "sex",
                    "affected",
                    "excluded",
                    "phecode_exclude",
                    "phecode_group",
                ]
            )
            .agg(
                {
                    "meaning": lambda x: list(x)  # Group meanings into a list
                }
            )
            .reset_index()
        )

        # Convert results to a list of dictionaries
        results = grouped.to_dict(orient="records")

    return jsonify(results)
