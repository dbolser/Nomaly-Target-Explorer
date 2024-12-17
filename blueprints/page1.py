from flask import render_template, request, jsonify, Blueprint
from flask import session, redirect, url_for
from db import get_db_connection
import pandas as pd

# ----------------------------------------------------- #
# Phecode Blueprint
# ----------------------------------------------------- #

# Create the blueprint
# page1_bp = Blueprint('page1', __name__)
disease_select = Blueprint('disease', __name__, template_folder='../templates')


# Route for page1 with tables
@disease_select.route('/page1')
def render_page1():
    return render_template('page1.html')
