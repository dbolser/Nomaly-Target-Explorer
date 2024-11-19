from flask import Flask, render_template, request, redirect, url_for, session, flash
from flask_mysqldb import MySQL
from flask_session import Session
import MySQLdb.cursors
import re

# Import blueprints after creating the app
from blueprints.main import main_bp
from blueprints.phecode import phecode_bp
from blueprints.phecode_term import phecode_term_bp
# from blueprints.page1 import page1_bp

app = Flask(__name__)

# Secret key for session management
app.secret_key = 'your_secret_key'

# Configure MySQL
app.config['MYSQL_HOST'] = 'localhost'
app.config['MYSQL_USER'] = 'clu'
app.config['MYSQL_PASSWORD'] = 'mysqlOutse3'
app.config['MYSQL_DB'] = 'ukbb'

mysql = MySQL(app)

# Configure session for access control
app.config['SESSION_TYPE'] = 'filesystem'
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

        cursor = mysql.connection.cursor(MySQLdb.cursors.DictCursor)
        cursor.execute('SELECT * FROM users WHERE username = %s AND password = %s', (username, password))
        user = cursor.fetchone()
        
        if user:
            session['loggedin'] = True
            session['id'] = user['id']
            session['username'] = user['username']
            return redirect(url_for('index'))
        else:
            flash('Incorrect username or password!')
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

# app.register_blueprint(page1_bp)
@app.route('/page1')
def page1():
    if 'loggedin' in session:
        return render_template('page1.html', username=session['username'])
    return redirect(url_for('login'))

@app.route('/page2')
def page2():
    if 'loggedin' in session:
        return render_template('page2.html', username=session['username'])
    return redirect(url_for('login'))

# Search routes for each search section
@app.route('/search1', methods=['GET'])
def search1():
    if 'loggedin' in session:
        # Process search query here
        return f"Search results for Section 1: {request.args.get('query1')}"
    return redirect(url_for('login'))

# Run app
if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8756)
