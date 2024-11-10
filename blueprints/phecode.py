from flask import Blueprint, render_template, request, jsonify, session
import threading
import time
from db import get_db_connection
import pandas as pd

# Create the blueprint
phecode_bp = Blueprint('phecode', __name__, template_folder='../templates')

# Dictionary to store background task results (in-memory for simplicity)
task_results = {}

# Route to render the Phecode page
@phecode_bp.route('/phecode/<string:phecode>', methods=['GET'])
def show_phecode(phecode):

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
    return render_template('phecode.html', data=data)

# Background task function
def background_task(phecode):
    # Simulate a time-consuming process
    time.sleep(5)
    result = f"Processed data for Phecode {phecode}"
    
    # Store the result in the task_results dictionary
    task_results[phecode] = result

# Endpoint to trigger the background task
@phecode_bp.route('/run-task/<string:phecode>', methods=['POST'])
def run_task(phecode):
    # Start the background task using threading
    task_thread = threading.Thread(target=background_task, args=(phecode,))
    task_thread.start()
    return jsonify({"status": "Task started"}), 202

# Endpoint to get the result of the background task
@phecode_bp.route('/task-result/<string:phecode>', methods=['GET'])
def get_task_result(phecode):
    result = task_results.get(phecode, "Processing...")
    return jsonify({"result": result})