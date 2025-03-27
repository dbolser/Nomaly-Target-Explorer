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
        SELECT u.id, u.username, u.email, u.is_active, u.is_admin, up.allowed_paths
        FROM users2 u
        LEFT JOIN user_permissions up ON u.id = up.user_id
    """)
    users = cursor.fetchall()

    # Get all available phecodes
    cursor.execute("""
        SELECT DISTINCT phecode, description 
        FROM phecode_definition
        ORDER BY phecode
    """)
    all_phecodes = cursor.fetchall()

    # Organize phecodes into a hierarchical structure
    hierarchical_phecodes = {}
    for phecode_data in all_phecodes:
        phecode = phecode_data["phecode"]
        # Check if it's a base code (no decimal) or a sub-code
        if "." not in phecode:
            # It's a base code
            if phecode not in hierarchical_phecodes:
                hierarchical_phecodes[phecode] = {
                    "description": phecode_data["description"],
                    "subcodes": [],
                }
        else:
            # It's a sub-code, get the base code
            base_code = phecode.split(".")[0]
            if base_code not in hierarchical_phecodes:
                # Create base code entry if it doesn't exist
                hierarchical_phecodes[base_code] = {
                    "description": f"Base: {base_code}",
                    "subcodes": [],
                }
            hierarchical_phecodes[base_code]["subcodes"].append(phecode_data)

    cursor.close()
    conn.close()

    return render_template(
        "admin/users.html",
        users=users,
        hierarchical_phecodes=hierarchical_phecodes,
        all_phecodes=all_phecodes,
    )


@admin_bp.route("/admin/permissions", methods=["POST"])
@admin_required
def update_permissions():
    user_id = request.form.get("user_id")
    phecodes = request.form.getlist("phecodes")
    is_admin = request.form.get("is_admin") == "on"

    # Join phecodes with commas for storage
    allowed_paths = ",".join(phecodes)

    conn = get_db_connection()
    cursor = conn.cursor()

    try:
        # Update admin status
        cursor.execute(
            "UPDATE users2 SET is_admin = %s WHERE id = %s", (is_admin, user_id)
        )

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
    is_admin = request.form.get("is_admin_new") == "on"

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
            INSERT INTO users2 (username, email, password, is_active, is_admin)
            VALUES (%s, %s, %s, TRUE, %s)
            """,
            (username, email, generate_password_hash(password), is_admin),
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


@admin_bp.route("/admin/toggle_user_status", methods=["POST"])
@admin_required
def toggle_user_status():
    user_id = request.form.get("user_id")

    if not user_id:
        flash("User ID is required", "error")
        return redirect(url_for("admin.manage_users"))

    conn = get_db_connection()

    try:
        # Toggle status directly in the SQL
        with conn.cursor() as cursor:
            # First check if user exists
            cursor.execute(
                "SELECT EXISTS(SELECT 1 FROM users2 WHERE id = %s)", (user_id,)
            )
            result = cursor.fetchone()
            user_exists = result[0] if result else False

            if not user_exists:
                flash("User not found", "error")
                return redirect(url_for("admin.manage_users"))

            # Toggle status
            cursor.execute(
                "UPDATE users2 SET is_active = NOT is_active WHERE id = %s", (user_id,)
            )

            # Check new status
            cursor.execute("SELECT is_active FROM users2 WHERE id = %s", (user_id,))
            result = cursor.fetchone()
            is_active = result[0] if result else False

        conn.commit()
        status_text = "activated" if is_active else "deactivated"
        flash(f"User {status_text} successfully!", "success")

    except Exception as e:
        conn.rollback()
        flash(f"Error updating user status: {str(e)}", "error")

    finally:
        conn.close()

    return redirect(url_for("admin.manage_users"))
