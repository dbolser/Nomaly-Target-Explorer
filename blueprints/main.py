from flask import render_template, request, jsonify, Blueprint
# from db import get_db_connection
# import pandas as pd

main_bp = Blueprint('main', __name__)

# Route to serve the HTML page
@main_bp.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')