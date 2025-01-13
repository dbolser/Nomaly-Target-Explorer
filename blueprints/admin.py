from flask import Blueprint, render_template, request, redirect, url_for, flash
from db import get_db_connection
from auth import admin_required
from werkzeug.security import generate_password_hash

admin_bp = Blueprint("admin", __name__)


@admin_bp.route("/admin/users")
@admin_required
def manage_users():
    conn = get_db_connection()
    cursor = conn.cursor(dictionary=True)

    # Get all users
    cursor.execute("""
        SELECT u.id, u.username, u.email, u.is_active, up.allowed_paths
        FROM users2 u
        LEFT JOIN user_permissions up ON u.id = up.user_id
    """)
    users = cursor.fetchall()

    # Get all available phecodes (you might want to store these in a separate table)
    cursor.execute("""
        SELECT DISTINCT phecode, description 
        FROM phecode_definition
        ORDER BY phecode
    """)
    available_phecodes = cursor.fetchall()

    cursor.close()
    conn.close()

    return render_template(
        "admin/users.html", users=users, available_phecodes=available_phecodes
    )


@admin_bp.route("/admin/permissions", methods=["POST"])
@admin_required
def update_permissions():
    user_id = request.form.get("user_id")
    phecodes = request.form.getlist("phecodes")

    # Join phecodes with commas for storage
    allowed_paths = ",".join(phecodes)

    conn = get_db_connection()
    cursor = conn.cursor()

    try:
        # Update or insert permissions
        cursor.execute(
            """
            INSERT INTO user_permissions (user_id, allowed_paths)
            VALUES (%s, %s)
            ON DUPLICATE KEY UPDATE allowed_paths = VALUES(allowed_paths)
        """,
            (user_id, allowed_paths),
        )

        conn.commit()
        flash("Permissions updated successfully!", "success")
    except Exception as e:
        conn.rollback()
        flash(f"Error updating permissions: {str(e)}", "error")
    finally:
        cursor.close()
        conn.close()

    return redirect(url_for("admin.manage_users"))


@admin_bp.route("/admin/users/add", methods=["POST"])
@admin_required
def add_user():
    username = request.form.get("username")
    email = request.form.get("email")
    password = request.form.get("password")

    if not all([username, email, password]):
        flash("All fields are required", "error")
        return redirect(url_for("admin.manage_users"))

    assert password is not None

    conn = get_db_connection()
    cursor = conn.cursor()

    try:
        # Check if username already exists
        cursor.execute("SELECT id FROM users2 WHERE username = %s", (username,))
        if cursor.fetchone():
            flash("Username already exists", "error")
            return redirect(url_for("admin.manage_users"))

        # Insert new user
        cursor.execute(
            """
            INSERT INTO users2 (username, email, password, is_active)
            VALUES (%s, %s, %s, TRUE)
            """,
            (username, email, generate_password_hash(password)),
        )
        conn.commit()
        flash(f"User {username} created successfully!", "success")

    except Exception as e:
        conn.rollback()
        flash(f"Error creating user: {str(e)}", "error")

    finally:
        cursor.close()
        conn.close()

    return redirect(url_for("admin.manage_users"))
