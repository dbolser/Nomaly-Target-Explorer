from flask import (
    Flask,
    render_template,
    request,
    redirect,
    url_for,
    session,
    flash,
    jsonify,
    Response,
)
from flask_mysqldb import MySQL
from flask_session import Session
import MySQLdb.cursors
from werkzeug.exceptions import HTTPException
from logging_config import setup_logging
from errors import DataNotFoundError, DatabaseConnectionError
from config import config
from auth import check_page_permission
import os

# Import blueprints after creating the app
from blueprints.main import main_bp
from blueprints.phecode import phecode_bp
from blueprints.phecode_term import phecode_term_bp
from blueprints.page1 import disease_select
from blueprints.variant import variant_bp

app = Flask(__name__)

# Load config based on environment
env = os.getenv("FLASK_ENV", "default")
app.config.from_object(config[env])

# Set up logging
logger = setup_logging(app)

# Configure MySQL
mysql = MySQL(app)

# Configure session for access control
Session(app)


# Home route
@app.route("/")
def index():
    """Home page that shows different content based on auth status."""
    if "loggedin" in session:
        return render_template(
            "index.html", username=session["username"], authenticated=True
        )
    return render_template("index.html", authenticated=False)


# Login route
@app.route("/login", methods=["GET", "POST"])
def login():
    # If user is already logged in, redirect to index
    if "loggedin" in session:
        return redirect(url_for("index"))

    if request.method == "POST":
        username = request.form["username"]
        password = request.form["password"]

        cursor = None
        try:
            cursor = mysql.connection.cursor(MySQLdb.cursors.DictCursor)  # type: ignore
            if cursor is None:
                raise DatabaseConnectionError("Could not establish database connection")

            cursor.execute(
                "SELECT * FROM users2 WHERE username = %s AND password = %s",
                (username, password),
            )
            user = cursor.fetchone()

            if user:
                session["loggedin"] = True
                session["id"] = user["id"]
                session["username"] = user["username"]
                return redirect(url_for("index"))
            else:
                flash("Incorrect username or password!")
        except MySQLdb.Error as e:
            app.logger.error(f"Database error: {str(e)}")
            flash("A database error occurred. Please try again later.")
        finally:
            if cursor:
                cursor.close()

    # GET request or failed POST request
    return render_template("login.html")


# Logout route
@app.route("/logout")
def logout():
    session.pop("loggedin", None)
    session.pop("id", None)
    session.pop("username", None)
    flash("You have been logged out successfully")
    return redirect(url_for("index"))


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


@app.before_request
def require_login():
    # Public routes that don't require authentication
    public_routes = ["login", "logout", "static"]
    public_paths = ["/", "/favicon.ico"]  # Root path should be public

    # Check if accessing public route or path
    if request.endpoint in public_routes or request.path in public_paths:
        return

    # Require login for everything else
    if "loggedin" not in session:
        flash("Please log in to access this page")
        return redirect(url_for("login"))

    # Check permissions using existing auth system
    user_id = session.get("id")
    if not check_page_permission(user_id, request.path):
        flash("Access denied: Insufficient permissions")
        return redirect(url_for("index"))


@app.route("/page2")
def page2():
    if "loggedin" in session:
        return render_template("page2.html", username=session["username"])
    return redirect(url_for("login"))


# Error handlers
@app.errorhandler(Exception)
def handle_exception(e: Exception) -> tuple[str | Response, int]:
    """Handle all unhandled exceptions"""
    if isinstance(e, HTTPException):
        response = {"error": True, "message": e.description, "status_code": e.code}
        status_code = e.code or 500
    else:
        app.logger.error(f"Unhandled exception: {str(e)}", exc_info=True)
        response = {
            "error": True,
            "message": "An unexpected error occurred",
            "status_code": 500,
        }
        status_code = 500

    if request.headers.get("Accept") == "application/json":
        return jsonify(response), status_code
    return render_template("error.html", error=response["message"]), status_code


@app.errorhandler(DataNotFoundError)
def handle_not_found(e):
    app.logger.warning(f"Data not found: {str(e)}")
    return render_template("error.html", error=str(e)), 404


@app.errorhandler(DatabaseConnectionError)
def handle_db_error(e):
    app.logger.error(f"Database error: {str(e)}")
    return render_template("error.html", error="Database connection error"), 503


@app.route("/search")
def search():
    """Disease search page."""
    if "loggedin" not in session:
        flash("Please log in to access the search page")
        return redirect(url_for("login"))
    return render_template("search.html", username=session["username"])


# Run app
if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8756)
