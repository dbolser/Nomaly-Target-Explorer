from flask import Flask, render_template, request, redirect, url_for, session, flash, jsonify, Response
from flask_mysqldb import MySQL
from flask_session import Session
import MySQLdb.cursors
from werkzeug.exceptions import HTTPException
from logging_config import setup_logging
from errors import *
from config import config
import os

# Import blueprints after creating the app
from blueprints.main import main_bp
from blueprints.phecode import phecode_bp
from blueprints.phecode_term import phecode_term_bp
from blueprints.page1 import disease_select
from blueprints.variant import variant_bp

app = Flask(__name__)

# Load config based on environment
env = os.getenv('FLASK_ENV', 'default')
app.config.from_object(config[env])

# Set up logging
logger = setup_logging(app)

# Configure MySQL
mysql = MySQL(app)

# Configure session for access control
Session(app)

# Home route
@app.route('/')
def index():
    if 'loggedin' in session:
        return render_template('index.html', username=session['username'])
    return redirect(url_for('login'))

# Login route
@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']

        try:
            cursor = mysql.connection.cursor(MySQLdb.cursors.DictCursor)
            if cursor is None:
                raise DatabaseConnectionError("Could not establish database connection")
                
            cursor.execute('SELECT * FROM users WHERE username = %s AND password = %s', (username, password))
            user = cursor.fetchone()
            
            if user:
                session['loggedin'] = True
                session['id'] = user['id']
                session['username'] = user['username']
                return redirect(url_for('index'))
            else:
                flash('Incorrect username or password!')
        except MySQLdb.Error as e:
            app.logger.error(f"Database error: {str(e)}")
            flash('A database error occurred. Please try again later.')
        finally:
            if 'cursor' in locals():
                cursor.close()
                
    return render_template('login.html')

# Logout route
@app.route('/logout')
def logout():
    session.pop('loggedin', None)
    session.pop('id', None)
    session.pop('username', None)
    return redirect(url_for('login'))

# Protected routes
# register main blueprint: /diseasesearch
app.register_blueprint(main_bp)

# register phecode blueprint: /phecode/<string:phecode>
app.register_blueprint(phecode_bp)

# register phecode_term_bp blueprint: /phecode/<string:phecode>/term/<string:term>  
app.register_blueprint(phecode_term_bp)

# register disease blueprint: /page1
app.register_blueprint(disease_select)

# register variant blueprint: /variant/<string:variant>
app.register_blueprint(variant_bp)

@app.route('/page2')
def page2():
    if 'loggedin' in session:
        return render_template('page2.html', username=session['username'])
    return redirect(url_for('login'))

# Error handlers
@app.errorhandler(Exception)
def handle_exception(e: Exception) -> tuple[str | Response, int]:
    """Handle all unhandled exceptions"""
    if isinstance(e, HTTPException):
        response = {
            "error": True,
            "message": e.description,
            "status_code": e.code
        }
        status_code = e.code
    else:
        app.logger.error(f"Unhandled exception: {str(e)}", exc_info=True)
        response = {
            "error": True,
            "message": "An unexpected error occurred",
            "status_code": 500
        }
        status_code = 500

    if request.headers.get('Accept') == 'application/json':
        return jsonify(response), status_code
    return render_template('error.html', error=response["message"]), status_code

@app.errorhandler(DataNotFoundError)
def handle_not_found(e):
    app.logger.warning(f"Data not found: {str(e)}")
    return render_template('error.html', error=str(e)), 404

@app.errorhandler(DatabaseConnectionError)
def handle_db_error(e):
    app.logger.error(f"Database error: {str(e)}")
    return render_template('error.html', error="Database connection error"), 503

# Run app
if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8756)
