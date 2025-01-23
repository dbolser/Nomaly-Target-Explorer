from functools import wraps
from flask import abort
from flask_login import current_user


def admin_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if not current_user.is_authenticated or not current_user.is_admin:
            abort(403)  # Forbidden
        return f(*args, **kwargs)

    return decorated_function


def check_page_permission(current_user, path):
    """Check if user has permission to access the page."""

    # Define public endpoints that don't require authentication
    public_paths = ["/", "/login", "/static", "/search", "/favicon.ico"]

    # Skip authentication for public routes
    if path in public_paths or path.startswith("/static/"):
        return True

    # All other routes require authentication
    if not current_user.is_authenticated:
        return False

    # Admin can access everything
    if current_user.is_admin:
        return True

    # Next we define what a non-admin user can access

    # Define basic routes that are always allowed once authenticated
    basic_routes = ["/logout", "/disease-sets", "/diseasesearch"]
    variant_routes = ["/variant", "/run-phewas", "/phewas-result"]

    # Check if the path starts with any of the basic or variant routes
    if any(path.startswith(route) for route in basic_routes + variant_routes):
        return True

    # Check phecode-specific routes
    phecode_routes = [
        "/phecode/",
        "/nomaly-stats/",
        "/phecode2/",
        "/nomaly-stats2/",
        "/run-task/",
        "/variant_scores/",
        "/stream_progress/",
    ]

    # Check if the path starts with any of the phecode-specific routes for the
    # current user
    for route in phecode_routes:
        if path.startswith(route):
            phecode = path[len(route) :].split("/")[0].split(".")[0]
            if phecode in current_user.allowed_paths:
                return True

    return False
