from functools import wraps
from flask import session, redirect, url_for, request
from db import get_db_connection
import MySQLdb.cursors


def check_page_permission(user_id, path):
    """Check if user has permission to access the page."""
    print(f"Checking permission for user {user_id} to access path: {path}")
    conn = get_db_connection()
    cursor = conn.cursor()

    # Get user's allowed phecodes and routes
    cursor.execute(
        """
        SELECT allowed_paths FROM user_permissions 
        WHERE user_id = %s
        """,
        (user_id,),
    )
    result = cursor.fetchone()
    cursor.close()
    conn.close()

    if not result:
        return False

    allowed_paths = result[0]

    # Admin check - full access
    if allowed_paths == "*":
        print("Admin access granted")
        return True

    allowed_paths = allowed_paths.split(",")

    # Check basic routes first
    if path in ["/page1", "/diseasesearch"]:
        print("Allowed basic routes")
        return True
    elif (
        False
        or path.startswith("/variant")
        or path.startswith("/run-phewas")
        or path.startswith("/phewas-result")
    ):
        print("Allowed variant and phewas routes")
        return True

    # For phecode routes, extract the phecode number
    for route in [
        "/phecode/",
        "/nomaly-stats/",
        "/phecode2/",
        "/nomaly-stats2/",
        "/run-task/",
    ]:
        if path.startswith(route):
            phecode = (
                path[len(route) :].split("/")[0].split(".")[0]
            )  # Get base phecode number
            if phecode in allowed_paths:
                return True

    return False
