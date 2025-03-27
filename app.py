import os

from flask import (
    Flask,
    Response,
    current_app,
    flash,
    jsonify,
    redirect,
    render_template,
    request,
    session,
    url_for,
)
from flask_cors import CORS
from flask_login import (
    LoginManager,
    UserMixin,
    current_user,
    login_user,
    logout_user,
)
from flask_mysqldb import MySQL
from werkzeug.exceptions import HTTPException
from werkzeug.security import check_password_hash

from auth import check_page_permission
from blueprints.admin import admin_bp
from blueprints.disease_sets import disease_sets_bp
from blueprints.phecode import phecode_bp
from blueprints.phecode_term import phecode_term_bp
from blueprints.prioritisation_by_nomaly_scores import prioritisation_bp

# Import blueprints after creating the app
from blueprints.search import search_bp
from blueprints.variant import variant_bp
from config import config
from db import get_db_connection
from errors import DatabaseConnectionError, DataNotFoundError
from flask_session import Session  # Server-side session extension
from dotenv import load_dotenv

# Initialize extensions
mysql = MySQL()
flask_session = Session()  # Renamed to avoid conflict with Flask's session
login_manager = LoginManager()
cors = CORS()


def create_app(config_name=None):
    """Create and configure the Flask application.

    Args:
        config_name: The name of the configuration to use (default, testing, production)
    """
    app = Flask(__name__)

    # Use FLASK_ENV if config_name not provided
    if config_name is None:
        config_name = os.getenv("FLASK_ENV", "development")

    print(f"Loading config: {config_name}")
    app.config.from_object(config[config_name])

    # Optionally load environment-specific .env file
    env_file = f".env.{config_name}"
    if os.path.exists(env_file):
        print(f"Loading environment file: {env_file}")
        load_dotenv(env_file)

    # Initialize extensions
    mysql.init_app(app)
    flask_session.init_app(app)  # Use the renamed variable
    login_manager.init_app(app)
    cors.init_app(app)

    from data_services.registry import ServiceRegistry

    ServiceRegistry(app)

    # Register blueprints
    register_blueprints(app)

    # Register routes and error handlers
    register_routes(app)
    register_error_handlers(app)

    return app


def register_blueprints(app):
    """Register Flask blueprints."""
    app.register_blueprint(search_bp)
    app.register_blueprint(phecode_bp)
    app.register_blueprint(phecode_term_bp)
    app.register_blueprint(disease_sets_bp)
    app.register_blueprint(variant_bp)
    app.register_blueprint(admin_bp)
    app.register_blueprint(prioritisation_bp)

    from api import stats_bp

    app.register_blueprint(stats_bp)


def register_routes(app):
    """Register all routes and related functions."""

    # User Class
    class User(UserMixin):
        def __init__(self, id, username, allowed_paths, is_admin=False):
            self.id = id
            self.username = username
            self.allowed_paths = allowed_paths.split(",") if allowed_paths else []
            self.is_admin = is_admin

    # User Loader Callback
    @login_manager.user_loader
    def load_user(user_id):
        conn = get_db_connection()
        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT
                u.username,
                IF(up.allowed_paths='*', 1, 0) AS is_admin,
                up.allowed_paths
            FROM users2 u
            LEFT JOIN user_permissions up ON u.id = up.user_id 
            WHERE u.id = %s AND u.is_active = TRUE
        """,
            (user_id,),
        )

        user = cursor.fetchone()
        cursor.close()
        conn.close()

        if user:
            return User(
                id=user_id,
                username=user[0],
                allowed_paths=user[2] or "",
                is_admin=bool(user[1]),
            )
        return None

    @app.before_request
    def require_permissions():
        """Check if user has permission to access the page."""
        if not check_page_permission(current_user, request.path):
            # If user is not logged in, redirect to login with the requested URL as 'next' parameter
            if not current_user.is_authenticated:
                return redirect(url_for("login", next=request.url))
            # If user is logged in but has insufficient permissions, go back to previous page
            flash("Access denied: Insufficient permissions")
            return redirect(request.referrer or url_for("index"))
        return None

    @app.route("/login", methods=["GET", "POST"])
    def login():
        if request.method == "POST":
            username = request.form.get("username")
            password = request.form.get("password")

            if not username or not password:
                flash("Please enter both username and password")
                # Preserve the 'next' parameter during form submission
                next_url = request.args.get("next")
                return redirect(
                    url_for("login", next=next_url) if next_url else url_for("login")
                )

            conn = get_db_connection()
            cursor = conn.cursor()
            cursor.execute(
                "SELECT id, password FROM users2 WHERE username = %s AND is_active = TRUE",
                (username,),
            )
            user = cursor.fetchone()
            cursor.close()
            conn.close()

            if user and check_password_hash(user[1], password):
                user_obj = load_user(user[0])
                if user_obj:
                    login_user(user_obj)
                    current_app.logger.info(f"User {user_obj.username} logged in.")
                    flash("Logged in successfully.")

                    # Get the next URL from the query parameters
                    next_page = request.args.get("next")

                    if next_page:
                        # Handle full URLs by extracting just the path
                        from urllib.parse import urlparse

                        parsed_url = urlparse(next_page)
                        # Use the path component if it's a complete URL
                        if parsed_url.netloc:
                            next_page = parsed_url.path

                        # Check if user has permission for the target page
                        if not check_page_permission(current_user, next_page):
                            flash("Access denied: Insufficient permissions")
                            return redirect(url_for("index"))
                    else:
                        next_page = url_for("index")
                    return redirect(next_page)
            flash("Invalid username or password.")

        # Pass next parameter to template
        next_url = request.args.get("next")
        return render_template("login.html", next=next_url)

    @app.route("/logout")
    def logout():
        logout_user()
        flash("You have been logged out successfully.")
        return redirect(url_for("index"))

    @app.route("/page2")
    def page2():
        return render_template("page2.html")

    @app.route("/search")
    def search():
        """Disease search page."""
        return render_template("search.html")

    @app.route("/")
    def index():
        """Home page route."""
        return render_template("index.html")


def register_error_handlers(app):
    """Register error handlers."""

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



# Only create app if running directly
if __name__ == "__main__":
    app = create_app(os.getenv("FLASK_ENV", "default"))
    app.run()
