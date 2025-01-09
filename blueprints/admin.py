from flask import Blueprint, render_template, request, redirect, url_for
from db import get_db_connection
import MySQLdb

admin_bp = Blueprint("admin", __name__)


@admin_bp.route("/admin/users", methods=["GET"])
def manage_users():
    conn = get_db_connection()
    cursor = conn.cursor(MySQLdb.cursors.DictCursor)  # type: ignore

    cursor.execute("SELECT * FROM users")
    users = cursor.fetchall()

    cursor.execute("SELECT * FROM pages")
    pages = cursor.fetchall()

    cursor.close()
    conn.close()

    return render_template("admin/users.html", users=users, pages=pages)


@admin_bp.route("/admin/permissions", methods=["POST"])
def update_permissions():
    user_id = request.form.get("user_id")
    page_ids = request.form.getlist("pages")
    wildcards = request.form.get("wildcards", "").split("\n")

    conn = get_db_connection()
    cursor = conn.cursor()

    # Clear existing permissions
    cursor.execute("DELETE FROM user_permissions WHERE user_id = %s", (user_id,))

    # Add new page permissions
    for page_id in page_ids:
        cursor.execute(
            """
            INSERT INTO user_permissions (user_id, page_id)
            VALUES (%s, %s)
        """,
            (user_id, page_id),
        )

    # Add wildcard permissions
    for wildcard in wildcards:
        if wildcard.strip():
            cursor.execute(
                """
                INSERT INTO user_permissions (user_id, wildcard_path)
                VALUES (%s, %s)
            """,
                (user_id, wildcard.strip()),
            )

    conn.commit()
    cursor.close()
    conn.close()

    return redirect(url_for("admin.manage_users"))
